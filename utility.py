from collections import namedtuple

def doctuple(doc, *args, **kwargs):
    # so we can use docstrings
    nt = namedtuple(*args, **kwargs)
    nt.__doc__ = doc
    return nt



