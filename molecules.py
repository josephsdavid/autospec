import numpy as np
from collections import namedtuple
from pprint import pprint
from bs4 import BeautifulSoup
import requests
from multiprocessing import pool
import re
import csv
from functools import reduce, partial

def doctuple(doc, *args, **kwargs):
    # so we can use docstrings
    nt = namedtuple(*args, **kwargs)
    nt.__doc__ = doc
    return nt


SpectralQuery = doctuple('''
                         documentation goes here
                         ''',
                         'SpectralQuery', ('molecule_name', 'data', 'frequency','uncertainty','intensity', 'units'))

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
    data = np.array([list(map(_regex_helper, x.split()[:3])) for x in lines[:-1]])
    return SpectralQuery(molecule, data, *data.T, 'MHz')


def get_molecule_csv():
    # put the csv online
    with open("all_molecules.csv",'r') as f:
        reader = csv.DictReader(f)
        out = {}
        for row in reader:
            for k, v in row.items():
                out.setdefault(k,[]).append(v)
    out['Molecule'] = out.pop('')
    return out

def construct_mega_query():
    # talk to partee
    molecules = get_molecule_csv()
    out = []
    for i in range(len(molecules['Molecule'])):
        print(i)
        info = [v[i] for v in list(molecules.values())]
        if info[1] in ["JPL","CDMS"]:
            data = get_molecule_data(info[1], info[0], info[-1])
            if data:
                out.append(data)
    return out
    #return SpectralQuery(*zip(*out))

