import asyncio
from functools import partial, reduce
from os.path import isfile
from random import randint
import re

from astropy import units as u
from astroquery.splatalogue import Splatalogue
from bs4 import BeautifulSoup
import numpy as np
from requests_html import AsyncHTMLSession
from tqdm import tqdm

from .spectraldata import convert_units
from .utility import doctuple

Splatalogue.QUERY_URL = 'https://splatalogue.online/c_export.php'


@asyncio.coroutine
async def _bs4_tool(url):
    session = AsyncHTMLSession()
    page = await session.get(url)
    await session.close()
    # Status code 200 means it was successful.
    if page.status_code == 200:
        #print(f"Finished {url}")
        return BeautifulSoup(page.html.raw_html, "html.parser")


@asyncio.coroutine
async def _get_jpl_data(url, tag):
    return str(await _bs4_tool(f"{url}/c{tag}.cat"))


def _td_finder(t, tag):
    return t.name == "td" and tag in t.text


def _a_finder(t):
    return t.name == "a" and "HTML" in t.text


@asyncio.coroutine
async def _get_cdms_data(url, tag):
    soup = await _bs4_tool(f"{url}/classic/entries/")
    td_find = partial(_td_finder, tag=tag)
    link = reduce(lambda x, cond: x.find(cond), ["table", td_find, _a_finder], soup)[
        "href"
    ]
    linepage = await _bs4_tool(f"{url}{link}")
    second_link = linepage.find("a")["href"]
    return str((await _bs4_tool(f"{url}{second_link}")).find("pre"))


def _regex_helper(s):
    match = re.match("^(?:[^\..]*[\..]){2}", s)
    try:
        return float(match.group(0)[:-1]) if match else float(s)
    except ValueError:
        # verify with an adult
        return float(re.sub("\D", "", s[: s.find("-")]))


SpectralQuery = doctuple(
    """documentation goes here""",
    "SpectralQuery",
    (
        "molecule",
        "data",
        "frequency",
        "uncertainty",
        "intensity",
        "deg",
        "elo",
        "units",
    ),
)


@asyncio.coroutine
async def get_molecule_data(db, tag, molecule=None):
    db_dict = {
        "JPL": {
            "url": "https://spec.jpl.nasa.gov/ftp/pub/catalog",
            "fun": _get_jpl_data,
        },
        "CDMS": {"url": "https://cdms.astro.uni-koeln.de", "fun": _get_cdms_data},
    }
    db = db.upper().strip()
    lines = await db_dict[db]["fun"](db_dict[db]["url"], tag)
    # regex genius
    lines = re.sub("(?i)[\n]?[\n]?</?pre[^>]*>[\n]?[\s]?", "", lines).split("\n")
    data = (
        np.array([list(map(_regex_helper, x.split()[:5])) for x in lines[:-1]]) * u.MHz
    )
    try:
        result = SpectralQuery(molecule, data.value, *data.value.T, "MHz")
    except TypeError:
        result = None
    return result


# def get_molecule_csv():
#    # put the csv online
#    with open("all_molecules.csv",'r') as f:
#        reader = csv.DictReader(f)
#        out = {}
#        for row in reader:
#            for k, v in row.items():
#                out.setdefault(k,[]).append(v)
#    out['Molecule'] = out.pop('')
#    return out


splat_query = doctuple(
    """placeholder""", "Splatalogue_Query", ("db", "tag", "molecule")
)


def query_splatalogue(spikes):
    # from a spike to the splatalogue (linelist, tag, moleculename)
    # verify this is right
    spacing = (spikes.data[1:, 0] - spikes.data[:-1, 0]).mean()
    bounds = (
        np.column_stack([spikes.frequency - spacing, spikes.frequency + spacing])
        * u.GHz
    )
    out = []
    for i in range(spikes.frequency.shape[0]):
        res = Splatalogue.query_lines(
            *bounds[i, :].T,
            show_molecule_tag=True,
            top20="planet",
            line_lists=["CDMS", "JPL"],
            line_strengths="ls1",
        )
        if "Molecule<br>Tag" in res.keys():
            tag = list(
                map(lambda x: str(x).zfill(6).replace("-", "0"), res["Molecule<br>Tag"])
            )
            db = list(res["Linelist"])
            molecule = list(res["Chemical Name"])
            if res:
                out.append(splat_query(db, tag, molecule))
    # tuple transpositions!
    return set(zip(*[sum(o, []) for o in zip(*out)]))


def _get_molecules_from_spikes(spikes):
    print("Querying SPlatalogue...")
    molecule_set = query_splatalogue(spikes)
    # Do the async thing
    jobs = [get_molecule_data(*m) for m in molecule_set]
    print("done!")
    loop = asyncio.get_event_loop()
    print("Query spectroscopy databases...")
    results = []
    for i in tqdm(range(0, len(jobs) + 5, 5)):
        results += loop.run_until_complete(asyncio.gather(*jobs[i : i + 5]))
    out = {m[-1]: result for m, result in zip(molecule_set, results) if result is not None}
    print("Cleaning up...")
    return {k: convert_units(v, "GHz") for k, v in out.items()}

def get_molecules_from_spikes(spikes, save_path=None):
    if save_path is not None:
        if isfile(save_path):
            return load_molecules(save_path)
        else:
            result = _get_molecules_from_spikes(spikes)
            save_molecules(result, save_path)
            return result
    else:
        return _get_molecules_from_spikes(spikes)


def save_molecules(mols, path):
    import joblib

    mols2 = {k: v._asdict() for k, v in mols.items()}
    joblib.dump(mols2, path)


def load_molecules(path):
    import joblib

    return {k: SpectralQuery(**v) for k, v in joblib.load(path).items()}
