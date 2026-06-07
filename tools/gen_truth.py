"""Dump exact ground-truth values from the Python reference so the Rust port
can assert byte-for-byte parity, including the 3 data-driven test failures."""
import json
import os, sys
_HERE = os.path.dirname(os.path.abspath(__file__))
_SUB = os.path.join(_HERE, '..', 'repos', 'pypgx')
os.chdir(_SUB)
import pypgx
import pandas as pd
import numpy as np
from fuc import pyvcf, common

out = {}

# ---- list_genes ----
out['genes_target'] = pypgx.list_genes(mode='target')
out['genes_control'] = pypgx.list_genes(mode='control')
out['genes_all'] = pypgx.list_genes(mode='all')

# ---- test_allele_table ----
at = pypgx.load_allele_table()
allele_dups = {}
allele_unsorted = {}
for assembly in ['GRCh37', 'GRCh38']:
    i = at[[f'{assembly}Core', 'SV']].dropna().duplicated(keep=False)
    l = sorted(at[[f'{assembly}Core', 'SV']].dropna()[i][f'{assembly}Core'].to_list())
    allele_dups[assembly] = l
    unsorted = []
    for _, r in at.iterrows():
        if pd.isna(r[f'{assembly}Core']):
            continue
        ordered = ','.join(sorted(r[f'{assembly}Core'].split(','), key=lambda x: common.parse_variant(x)[1]))
        if r[f'{assembly}Core'] != ordered:
            unsorted.append([r[f'{assembly}Core'], ordered])
    allele_unsorted[assembly] = unsorted
out['allele_dups'] = allele_dups
out['allele_unsorted'] = allele_unsorted
out['allele_gene_unique'] = list(at.Gene.unique())
out['list_genes_default'] = pypgx.list_genes()

# ---- test_diplotype_table ----
d1 = pypgx.load_diplotype_table()
gt = pypgx.load_gene_table()
out['diplotype_gene_unique_count'] = int(len(d1.Gene.unique()))
out['gene_phenotypemethod_diplotype_count'] = int(gt.PhenotypeMethod.value_counts()['Diplotype'])

# ---- test_equation_table ----
e1 = pypgx.load_equation_table()
out['equation_gene_unique_count'] = int(len(e1.Gene.unique()))
out['gene_phenotypemethod_score_count'] = int(gt.PhenotypeMethod.value_counts()['Score'])

# ---- test_priority_table ----
ph = pypgx.load_phenotype_table()
out['priority_a'] = [x for x in pypgx.list_genes() if pypgx.has_phenotype(x)]
out['priority_b'] = list(ph.Gene.unique())

# ---- test_definition_table ----
vt = pypgx.load_variant_table()
# part 1: variant-table self consistency (collect violations, don't raise)
def one_row(r):
    bad = []
    for assembly in ['GRCh37', 'GRCh38']:
        other = 'GRCh38' if assembly == 'GRCh37' else 'GRCh37'
        variant = r[f'{assembly}Name']
        if pd.isna(variant):
            continue
        chrom, pos, ref, alt = common.parse_variant(variant)
        if not (chrom == r.Chromosome and pos == r[f'{assembly}Position'] and
                ref == r[f'{assembly}Allele'] and (alt == r.Variant or alt == r[f'{other}Allele'])):
            bad.append(variant)
    return bad
part1_bad = []
for _, r in vt.iterrows():
    part1_bad += one_row(r)
out['definition_part1_bad'] = part1_bad
# part 2: allele-table vs variant-table symmetric diff per gene/assembly
diffs = []
for gene in pypgx.list_genes():
    t1 = at[at.Gene == gene]
    t2 = vt[vt.Gene == gene]
    for assembly in ['GRCh37', 'GRCh38']:
        variants = []
        for _, r in t1.iterrows():
            if not pd.isna(r[f'{assembly}Core']):
                for v in r[f'{assembly}Core'].split(','):
                    if v not in variants:
                        variants.append(v)
            if not pd.isna(r[f'{assembly}Tag']):
                for v in r[f'{assembly}Tag'].split(','):
                    if v not in variants:
                        variants.append(v)
        s = t2[f'{assembly}Name'].unique()
        diff = set(variants) ^ set(s[~pd.isna(s)])
        if diff:
            diffs.append([gene, assembly, sorted(diff)])
out['definition_diffs'] = diffs

# ---- test_predict_alleles ----
for tag in ['GRCh37', 'GRCh38']:
    a = pypgx.predict_alleles(f'test-data/CYP4F2-{tag}.zip')
    out[f'predict_{tag}'] = {idx: a.data.loc[idx].to_list() for idx in a.data.index}

# ---- supporting: build_definition_table + helpers for CYP4F2 ----
for tag in ['GRCh37', 'GRCh38']:
    vf = pypgx.build_definition_table('CYP4F2', tag)
    out[f'deftable_CYP4F2_{tag}'] = {
        'samples': list(vf.samples),
        'df': vf.df.astype(str).to_dict(orient='records'),
    }
    out[f'list_variants_CYP4F2_{tag}'] = pypgx.list_variants('CYP4F2', assembly=tag)
out['ref_allele_CYP4F2'] = pypgx.get_ref_allele('CYP4F2')
out['default_allele_CYP4F2_GRCh37'] = pypgx.get_default_allele('CYP4F2', 'GRCh37')
out['default_allele_CYP4F2_GRCh38'] = pypgx.get_default_allele('CYP4F2', 'GRCh38')
out['synonyms_CYP4F2_GRCh37'] = pypgx.get_variant_synonyms('CYP4F2', 'GRCh37')

# ---- parse_variant / sort_variants examples ----
pv = ['19-16008388-A-C', '2-234668879-C-CAT', '22-42127941', '22-42127941-G',
      '2-234668879---AT', '1:100:A>C']
out['parse_variant'] = {x: list(common.parse_variant(x)) for x in pv}
sv = ['5-200-G-T', '5:100:T:C', '1:100:A>C', '10-100-G-C', '19-16008388-A-C', '19-15990431-C-T']
out['sort_variants'] = common.sort_variants(set(sv))

print(json.dumps(out, indent=1, default=str))
