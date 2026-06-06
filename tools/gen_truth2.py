"""Ground truth for the additional pure core functions (breadth parity)."""
import json
import os, sys
_HERE = os.path.dirname(os.path.abspath(__file__))
_SUB = os.path.join(_HERE, '..', 'repos', 'pypgx')
os.chdir(_SUB)
import math
import pypgx

def fnum(x):
    # Represent floats/NaN as strings so the Rust side can compare exactly.
    if isinstance(x, float) and math.isnan(x):
        return "nan"
    return str(x)

out = {}

# predict_score
ps = [
    ('CYP2D6', '*1'), ('CYP2D6', '*1x2'), ('CYP2D6', '*1x4'), ('CYP2D6', '*4'),
    ('CYP2D6', '*4x2'), ('CYP2D6', '*22'), ('CYP2D6', '*22x2'),
    ('CYP2D6', '*36+*10'), ('CYP2D6', '*1x2+*4x2+*10'),
    ('DPYD', 'Reference'), ('DPYD', 'c.1905+1G>A (*2A)'),
    ('CYP2B6', '*1'), ('CYP2B6', '*2'), ('CYP2C9', '*1'), ('CYP2C9', '*2'), ('CYP2C9', '*3'),
]
out['predict_score'] = {f'{g}|{a}': fnum(pypgx.predict_score(g, a)) for g, a in ps}

# predict_phenotype
pp = [
    ('CYP2D6', '*4', '*5'), ('CYP2D6', '*5', '*4'), ('CYP2D6', '*1', '*22'),
    ('CYP2D6', '*1', '*1x2'), ('CYP2B6', '*1', '*4'), ('CYP2D6', '*1', '*1'),
    ('CYP2C9', '*1', '*1'), ('CYP2C9', '*2', '*3'), ('CYP2C19', '*1', '*1'),
    ('DPYD', 'Reference', 'Reference'),
    ('CACNA1S', 'Reference', 'Reference'), ('CACNA1S', 'Reference', 'c.520C>T'),
]
out['predict_phenotype'] = {f'{g}|{a}|{b}': pypgx.predict_phenotype(g, a, b) for g, a, b in pp}

# get_priority
gp = [
    ('CYP2D6', 'Normal Metabolizer'), ('CYP2D6', 'Ultrarapid Metabolizer'),
    ('CYP3A5', 'Normal Metabolizer'), ('CYP3A5', 'Poor Metabolizer'),
]
out['get_priority'] = {f'{g}|{p}': pypgx.get_priority(g, p) for g, p in gp}

# get_region / exons / strand / paralog
genes = ['CYP2D6', 'CYP4F2', 'CYP2B6', 'CFTR']
out['get_region'] = {f'{g}|{asm}': pypgx.get_region(g, assembly=asm)
                     for g in genes for asm in ['GRCh37', 'GRCh38']}
out['get_exon_starts'] = {f'{g}|{asm}': pypgx.get_exon_starts(g, assembly=asm)
                          for g in genes for asm in ['GRCh37', 'GRCh38']}
out['get_exon_ends'] = {f'{g}|{asm}': pypgx.get_exon_ends(g, assembly=asm)
                        for g in genes for asm in ['GRCh37', 'GRCh38']}
out['get_strand'] = {g: pypgx.get_strand(g) for g in genes}
out['get_paralog'] = {g: pypgx.get_paralog(g) for g in genes + ['CYP2E1']}

# list_functions (NaN -> None -> "nan")
def lf(v):
    return ["nan" if (isinstance(x, float) and math.isnan(x)) else x for x in v]
out['list_functions'] = {g: lf(pypgx.list_functions(gene=g)) for g in ['CYP2D6', 'CYP4F2']}
out['list_functions_all'] = lf(pypgx.list_functions())

# list_phenotypes
out['list_phenotypes'] = {g: pypgx.list_phenotypes(gene=g) for g in ['CYP2D6', 'CYP2C9']}
out['list_phenotypes_all'] = pypgx.list_phenotypes()

# is_legit_allele / has_sv(allele) / has_score / get_score / get_function
out['is_legit_allele'] = {f'CYP2D6|{a}': bool(pypgx.is_legit_allele('CYP2D6', a))
                          for a in ['*1', '*4', '*999']}
out['has_score'] = {g: bool(pypgx.has_score(g)) for g in ['CYP2D6', 'CYP2B6', 'CYP4F2']}
out['has_sv'] = {g: bool(pypgx.has_sv(g)) for g in ['CYP2D6', 'CYP3A5', 'CYP4F2']}
out['get_function'] = {f'CYP2D6|{a}': ("nan" if (isinstance(pypgx.get_function('CYP2D6', a), float)) else pypgx.get_function('CYP2D6', a))
                       for a in ['*1', '*4', '*22']}
out['get_score'] = {f'CYP2D6|{a}': fnum(pypgx.get_score('CYP2D6', a)) for a in ['*1', '*4', '*22']}

# get_recommendation
import warnings
gr_cases = [
    ('codeine', 'CYP2D6', 'Normal Metabolizer', None, None),
    ('codeine', 'CYP2D6', 'Ultrarapid Metabolizer', None, None),
    ('codeine', 'CYP2D6', 'Poor Metabolizer', None, None),
    ('codeine', 'CYP2D6', 'Indeterminate', None, None),
    ('tacrolimus', 'CYP3A5', 'Normal Metabolizer', None, None),
    ('fluvastatin', 'CYP2C9', 'Normal Metabolizer', 'SLCO1B1', 'Normal Function'),
    ('fluvastatin', 'SLCO1B1', 'Normal Function', 'CYP2C9', 'Normal Metabolizer'),
]
gr = {}
for c in gr_cases:
    args = [x for x in c if x is not None]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        gr['|'.join(args)] = pypgx.get_recommendation(*args)
out['get_recommendation'] = gr

print(json.dumps(out, indent=1))
