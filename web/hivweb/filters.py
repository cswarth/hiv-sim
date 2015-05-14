import itertools

CLASSES = {
        'A': 'adenine',
        'C': 'cytosine',
        'G': 'guanine',
        'T': 'thymine',
        '-': 'gap',
        '.': 'match'}

def group_by_base(sequence):
    """
    Split a string into continuous identical sections, e.g.
    >>> group_by_base('ACCTTTGA')
    ['A', 'CC', 'TTT', 'G', 'A']
    """
    sequence = str(sequence)
    return [(CLASSES.get(k, ''), ''.join(v)) for k, v in itertools.groupby(sequence, lambda x: x)]


def percent_filter(f, default=""):
    if f is None:
        return default
    return "{0:0.2f}%".format(f*100)

def register(app):
    """
    Register all filters with an application
    """
    app.jinja_env.filters['percent'] = percent_filter
    app.jinja_env.filters['group_by_base'] = group_by_base
