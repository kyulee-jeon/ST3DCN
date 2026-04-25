"""Per-essential-tag heterogeneity analysis.

For each of 26 essential CT tags:
  1. Bin values (categorical passthrough or numeric cut)
  2. Group lesions by bin
  3. Compute within-bin metric (sensitivity for HCC+ cohort, specificity for
     HCC− cohort) and Wilson 95% CI
  4. Kruskal-Wallis test on p(HCC) distributions across bins with n ≥ 5
  5. Classify tag: fully_missing / invariant / insufficient / testable
  6. BH-FDR adjustment over testable tags (α = 0.05)

Outputs:
  - classification.csv     — one row per essential tag
  - per_bin.csv             — all bins × all tags (wide for auditing)
  - report.md               — readable writeup

Supports direct-only and derivation-augmented modes via substitute columns.
"""
import re
from typing import Optional

import numpy as np
import pandas as pd
from scipy.stats import kruskal

from .binners import BINNERS, normalize_cat
from .essential_tags import ESSENTIAL_26


DERIVATION_MAP = {
    # direct column → derived column to substitute when direct is mostly missing
    'SpacingBetweenSlices':  'SpacingBetweenSlices_derived',
    'Exposure_uAs':          'Exposure_uAs_derived',
    'Manufacturer':          'Manufacturer_derived',
    'ManufacturerModelName': 'ManufacturerModelName_derived',
}


def wilson(p: float, n: int, z: float = 1.96):
    if n == 0: return (0.0, 0.0)
    d = 1 + z**2 / n
    c = (p + z**2 / (2*n)) / d
    h = z * ((p*(1-p)/n + z**2/(4*n**2))**0.5) / d
    return max(0.0, c-h), min(1.0, c+h)


def bh_fdr(pvals):
    """Benjamini-Hochberg FDR-adjusted p-values."""
    pvals = np.asarray(pvals, dtype=float)
    mask = ~np.isnan(pvals)
    p = pvals[mask]
    if len(p) == 0: return np.full_like(pvals, np.nan)
    order = np.argsort(p); ranked = p[order]; n = len(p)
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out = np.full(n, np.nan); out[order] = adj
    fdr = np.full_like(pvals, np.nan); fdr[mask] = out
    return fdr


def classify(tag_rows, fill_rate, n_testable_bins):
    if fill_rate < 0.05:
        return 'fully_missing'
    present = [r for r in tag_rows if r['bin'] not in ('(missing)','(column missing)')]
    if len({r['bin'] for r in present}) <= 1:
        return 'invariant'
    if n_testable_bins < 2:
        return 'insufficient'
    return 'testable'


def _apply_derivations(df: pd.DataFrame):
    """In-place: replace direct column with derived when direct fill <5%."""
    applied = {}
    for direct, derived in DERIVATION_MAP.items():
        if derived not in df.columns:
            continue
        direct_vals = df.get(direct)
        direct_fill = 0
        if direct_vals is not None:
            non_missing = (direct_vals.notna()
                           & (direct_vals.astype(str).str.strip() != '')
                           & (direct_vals.astype(str).str.lower() != 'nan'))
            direct_fill = int(non_missing.sum())
        N = len(df)
        derived_fill = int(df[derived].notna().sum())
        if derived_fill > direct_fill:
            if direct_fill / max(1, N) < 0.05 and derived_fill > 0:
                df[direct] = df[derived]
                applied[direct] = f'substituted from {derived} (direct fill {direct_fill}/{N})'
            else:
                df[direct] = df[direct].where(
                    df[direct].notna() & (df[direct].astype(str).str.strip() != ''),
                    df[derived])
                applied[direct] = f'filled-in from {derived} where direct missing'
    return applied


def analyze(
    in_csv: str,
    out_cls_csv: str,
    out_bin_csv: str,
    out_md: str,
    threshold: float = 0.8,
    all_positive_label: bool = True,
    apply_derivations: bool = True,
    min_bin_n: int = 5,
    fdr_alpha: float = 0.05,
):
    """all_positive_label=True  → cohort is HCC+ (sensitivity primary).
    all_positive_label=False → cohort is HCC− (specificity primary).
    """
    df = pd.read_csv(in_csv)
    work = df[df['status'] == 'ok'].copy()
    work['prob'] = work['prob'].astype(float)
    work['pred'] = (work['prob'] >= threshold).astype(int)
    N = len(work)

    applied_note = _apply_derivations(work) if apply_derivations else {}

    rows_all = []; cls_rows = []; tag_tag_rows = {}

    for col, hex_tag, cls_code, cls_name, bin_type, mentioned in ESSENTIAL_26:
        binner = BINNERS.get(bin_type, normalize_cat)
        vals = work[col].apply(binner) if col in work.columns \
               else pd.Series(['(column missing)']*N)
        g = work.assign(_bin=vals.values)
        grouped = g.groupby('_bin').agg(
            n=('pred', 'size'), hit=('pred', 'sum'),
            mean_prob=('prob', 'mean'), median_prob=('prob', 'median')).reset_index()

        def sort_key(s):
            m = re.match(r'\[([-\d.]+),', str(s))
            if m: return (0, float(m.group(1)))
            try: return (1, float(s))
            except (ValueError, TypeError): return (2, s)
        grouped = grouped.assign(_sk=grouped['_bin'].map(sort_key))\
                           .sort_values('_sk').drop(columns='_sk')

        tag_rows = []
        for _, r in grouped.iterrows():
            n = int(r['n']); hit = int(r['hit'])
            if all_positive_label:
                sens = hit / n if n else 0.0
                lo, hi = wilson(sens, n)
                rec = {'bin': r['_bin'], 'n': n, 'tp': hit, 'sens': sens,
                       'ci_lo': lo, 'ci_hi': hi,
                       'mean_prob': float(r['mean_prob']),
                       'median_prob': float(r['median_prob'])}
            else:
                tn = n - hit; spec = tn / n if n else 0.0
                lo, hi = wilson(spec, n)
                rec = {'bin': r['_bin'], 'n': n, 'fp': hit, 'tn': tn, 'spec': spec,
                       'ci_lo': lo, 'ci_hi': hi,
                       'mean_prob': float(r['mean_prob']),
                       'median_prob': float(r['median_prob'])}
            tag_rows.append(rec)
            rows_all.append({'class': cls_name, 'tag': col, 'hex': hex_tag,
                              'mentioned': mentioned,
                              'derived': applied_note.get(col, ''), **rec})
        tag_tag_rows[col] = tag_rows

        present = [r for r in tag_rows if r['bin'] not in ('(missing)','(column missing)')]
        fill = sum(r['n'] for r in present) / N if N else 0.0
        testable_bins = [r for r in present if r['n'] >= min_bin_n]
        n_testable = len(testable_bins)

        kw_H = kw_p = np.nan
        if n_testable >= 2:
            groups = [g[g['_bin'] == r['bin']]['prob'].values for r in testable_bins]
            try: kw_H, kw_p = kruskal(*groups)
            except Exception: pass

        cls_rows.append({
            'tag': col, 'hex': hex_tag, 'class': cls_name,
            'mentioned_in_paper': mentioned,
            'derived_applied': applied_note.get(col, ''),
            'fill_rate': fill,
            'n_unique_bins': len({r['bin'] for r in present}),
            'n_testable_bins': n_testable,
            'status': classify(tag_rows, fill, n_testable),
            'kw_H': kw_H, 'kw_p': kw_p,
            'mean_min': min((r['mean_prob'] for r in testable_bins), default=np.nan),
            'mean_max': max((r['mean_prob'] for r in testable_bins), default=np.nan),
        })

    cls_df = pd.DataFrame(cls_rows)
    cls_df['kw_p_fdr'] = bh_fdr(cls_df['kw_p'].values)
    cls_df['significant_fdr'] = (cls_df['kw_p_fdr'] < fdr_alpha).fillna(False)

    cls_df.to_csv(out_cls_csv, index=False)
    pd.DataFrame(rows_all).to_csv(out_bin_csv, index=False)

    # Minimal markdown report
    label = 'sensitivity' if all_positive_label else 'specificity'
    N_total = len(work)
    hit_total = int(work['pred'].sum())
    base = hit_total / N_total if all_positive_label else (N_total - hit_total) / N_total
    blo, bhi = wilson(base, N_total)
    with open(out_md, 'w') as f:
        f.write(f'# Per-tag {label} + probability heterogeneity\n\n')
        f.write(f'- **N lesions**: {N_total}\n')
        f.write(f'- **Threshold on p(HCC)**: {threshold}\n')
        f.write(f'- **Baseline {label}**: {base:.4f} (Wilson 95% CI {blo:.4f}, {bhi:.4f})\n')
        f.write(f'- **BH-FDR α**: {fdr_alpha}\n\n')
        if applied_note:
            f.write('## Derivation substitutions applied\n\n')
            for t, note in applied_note.items():
                f.write(f'- **{t}**: {note}\n')
            f.write('\n')
        f.write('## Classification summary\n\n')
        counts = cls_df['status'].value_counts().to_dict()
        f.write(f'| status | n | notes |\n|---|---|---|\n')
        for s in ['testable','insufficient','invariant','fully_missing']:
            f.write(f'| {s} | {counts.get(s,0)} | |\n')
        f.write(f'\nFDR-significant (α={fdr_alpha}): {int(cls_df["significant_fdr"].sum())}\n\n')
        sig = cls_df[cls_df['significant_fdr']].sort_values('kw_p')
        if len(sig):
            f.write('## FDR-significant tags\n\n')
            f.write('| tag | mentioned | derived | raw p | FDR p |\n|---|---|---|---|---|\n')
            for _, r in sig.iterrows():
                f.write(f"| {r['tag']} | "
                        f"{'✓' if r['mentioned_in_paper'] else '✗'} | "
                        f"{'◆' if str(r['derived_applied']).strip() else ''} | "
                        f"{r['kw_p']:.4f} | **{r['kw_p_fdr']:.4f}** |\n")
    print(f'classification → {out_cls_csv}')
    print(f'per-bin         → {out_bin_csv}')
    print(f'markdown        → {out_md}')
    return cls_df
