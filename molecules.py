from spectraldata import convert_units
import numpy as np
from collections import namedtuple
from pprint import pprint
from bs4 import BeautifulSoup
import requests
from multiprocessing import pool
import re
import csv
from astropy import units as u
from astroquery.splatalogue import Splatalogue
from utility import doctuple
from functools import reduce, partial


def _bs4_tool(url):
    return BeautifulSoup(requests.get(url).content, 'html.parser')

def _get_jpl_data(url, tag):
    return str(_bs4_tool(f"{url}/c{tag}.cat"))

def _td_finder(t, tag):
    return t.name == 'td' and tag in t.text

def _a_finder(t):
    return t.name == 'a' and 'HTML' in t.text


def _get_cdms_data(url, tag):
    soup = _bs4_tool(f"{url}/classic/entries/")
    td_find = partial(_td_finder, tag=tag)
    link = reduce(lambda x, cond: x.find(cond), ["table" ,td_find, _a_finder], soup)['href']
    linepage = _bs4_tool(f"{url}{link}")
    second_link = linepage.find("a")['href']
    return str((_bs4_tool(f"{url}{second_link}")).find("pre"))

def _regex_helper(s):
    match = re.match("^(?:[^\..]*[\..]){2}",s)
    try:
        return float(match.group(0)[:-1]) if match else float(s)
    except ValueError:
        # verify with an adult
        return float(s[:s.find("-")])

SpectralQuery = doctuple(
    """documentation goes here""",
    "SpectralQuery",
    (
        "molecule",
        "data",
        "frequency",
        "uncertainty",
        "intensity",
        "units"
    )
)


def get_molecule_data(db, tag, molecule = None):
    db_dict = {
        "JPL":{
            "url":"https://spec.jpl.nasa.gov/ftp/pub/catalog",
            "fun":_get_jpl_data
        },
        "CDMS":{
            "url":"https://cdms.astro.uni-koeln.de",
            "fun":_get_cdms_data
        }
    }
    db = db.upper().strip()
    lines = db_dict[db]['fun'](db_dict[db]['url'], tag)
    # regex genius
    lines = re.sub("(?i)[\n]?[\n]?</?pre[^>]*>[\n]?[\s]?","",lines).split("\n")
    data = np.array([list(map(_regex_helper, x.split()[:3])) for x in lines[:-1]]) * u.MHz
    return SpectralQuery(molecule, data.value, *data.value.T, 'MHz')
#def get_molecule_csv():
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
    """placeholder""","Splatalogue_Query",('db','tag','molecule')
)
def query_splatalogue(spikes):
    # from a spike to the splatalogue (linelist, tag, moleculename)
    # verify this is right
    spacing = (spikes.data[1:, 0] - spikes.data[:-1, 0]).mean() / 2
    bounds = np.column_stack([spikes.frequency - spacing, spikes.frequency + spacing]) * u.GHz
    out = []
    for i in range(spikes.frequency.shape[0]):
        res = Splatalogue.query_lines(
            *bounds[i,:].T,
            show_molecule_tag = True,
            top20='planet',
            line_lists = ['CDMS', 'JPL'],
            line_strengths = 'ls1'
        )
        tag = list(map(lambda x: str(x).zfill(6).replace('-','0'), res['Molecule<br>Tag']))
        db = list(res['Linelist'])
        molecule = list(res['Chemical Name'])
        if res:
            out.append(splat_query(db, tag, molecule))
    # tuple transpositions!
    return set(zip(*[sum(o, []) for o in zip(*out)]))

def get_molecules_from_spikes(spikes):
    molecule_set = query_splatalogue(spikes)
    out = {}
    for m in molecule_set:
        out[m[-1]] = get_molecule_data(*m)
    return {k: convert_units(v, 'GHz') for k, v in out.items()}
