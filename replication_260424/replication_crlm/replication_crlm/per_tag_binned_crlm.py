"""Per-essential-tag SPECIFICITY + probability analysis on CRLM cohort.

Mirrors /home/ubuntu/hcc_workspace/replication/per_tag_binned.py but adapted for
CRLM's 100% negative cohort:
  - Primary metric per bin: specificity (= TN/n = fraction with prob < threshold)
  - Secondary: mean prob, Wilson 95% CI on specificity
  - Also performs Kruskal-Wallis + BH-FDR on probability distributions across bins

Outputs:
  - per_tag_binned.csv
  - per_tag_binned.md
  - per_tag_classification.csv (fully_missing / invariant / insufficient / testable / significant)
"""
import re

import numpy as np
import pandas as pd
from scipy.stats import kruskal

IN_CSV = '/home/ubuntu/hcc_workspace/replication_crlm/batch_results_with_tags.csv'
OUT_MD = '/home/ubuntu/hcc_workspace/replication_crlm/per_tag_binned.md'
OUT_CSV = '/home/ubuntu/hcc_workspace/replication_crlm/per_tag_binned.csv'
OUT_CLS = '/home/ubuntu/hcc_workspace/replication_crlm/per_tag_classification.csv'
THRESHOLD = 0.8

ESSENTIAL = [
    ('KVP',                     '180060', 'i', 'Intensity/exposure', 'cat', True),
    ('XRayTubeCurrent',         '181151', 'i', 'Intensity/exposure', 'num_tube_current', True),
    ('ConvolutionKernel',       '181210', 'i', 'Intensity/exposure', 'kernel_family', False),
    ('ContrastBolusAgent',      '180010', 'i', 'Intensity/exposure', 'cat', True),
    ('CTDIvol',                 '189345', 'i', 'Intensity/exposure', 'cat', False),
    ('ContrastFlowRate',        '181046', 'i', 'Intensity/exposure', 'cat', True),
    ('Exposure_uAs',            '181153', 'i', 'Intensity/exposure', 'num_exposure_uas', False),
    ('RevolutionTime',          '189305', 'i', 'Intensity/exposure', 'num_revtime', False),
    ('ReconstructionAlgorithm', '189315', 'i', 'Intensity/exposure', 'cat', False),
    ('ContrastBolusAgentPhase', '189344', 'i', 'Intensity/exposure', 'cat', True),
    ('SliceThickness',          '180050', 'g', 'Geometry',           'cat', True),
    ('Rows',                    '280010', 'g', 'Geometry',           'cat', False),
    ('Columns',                 '280011', 'g', 'Geometry',           'cat', False),
    ('PixelSpacing',            '280030', 'g', 'Geometry',           'num_pixspacing', False),
    ('ReconstructionDiameter',  '181100', 'g', 'Geometry',           'num_recondiam', False),
    ('WindowCenter',            '281050', 'g', 'Geometry',           'cat_multival', False),
    ('WindowWidth',             '281051', 'g', 'Geometry',           'cat_multival', False),
    ('TotalCollimationWidth',   '189307', 'g', 'Geometry',           'num_tcw', False),
    ('SpiralPitchFactor',       '189311', 'g', 'Geometry',           'num_pitch', True),
    ('TableSpeed',              '189309', 'g', 'Geometry',           'num_tablespeed', False),
    ('SpacingBetweenSlices',    '180088', 'g', 'Geometry',           'num_spacing', False),
    ('ReformattingThickness',   '720512', 'g', 'Geometry',           'cat', True),
    ('ManufacturerModelName',   '081090', 'd', 'Device',             'cat', True),
    ('Manufacturer',            '080070', 'd', 'Device',             'cat', True),
    ('BodyPartExamined',        '180015', 'a', 'Anatomy/position',   'cat', True),
    ('PatientPosition',         '185100', 'a', 'Anatomy/position',   'cat', False),
]


def wilson(p, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    denom = 1 + z**2 / n
    centre = (p + z**2 / (2*n)) / denom
    half = z * ((p*(1-p)/n + z**2/(4*n**2))**0.5) / denom
    return max(0.0, centre - half), min(1.0, centre + half)


def try_float(v):
    if v is None:
        return None
    try:
        if pd.isna(v): return None
    except (TypeError, ValueError):
        pass
    try:
        x = float(str(v).strip("'\""))
        if np.isnan(x): return None
        return x
    except (ValueError, TypeError, AttributeError):
        return None


def parse_pixspacing_mean(v):
    if pd.isna(v): return None
    s = str(v)
    nums = [float(x) for x in re.findall(r'\d+\.?\d*', s)]
    if not nums: return None
    return sum(nums) / len(nums)


def _num_canon(s):
    try:
        x = float(s)
        if x == int(x): return str(int(x))
        return f'{x:g}'
    except (ValueError, TypeError):
        return s


def normalize_multival_cat(v):
    if pd.isna(v): return '(missing)'
    s = str(v).strip().strip("'\"")
    if not s or s.lower() == 'nan': return '(missing)'
    m = re.match(r'\[\s*([^,\]]+)\s*(?:,\s*([^,\]]+))?\s*\]', s)
    if m:
        parts = [p.strip() for p in m.groups() if p]
        if parts and all(p == parts[0] for p in parts):
            return _num_canon(parts[0])
    return _num_canon(s)


def kernel_family(v):
    if pd.isna(v): return '(missing)'
    s = str(v).upper().strip().strip("'\"")
    if not s: return '(missing)'
    for fam in ('STANDARD', 'SOFT', 'B30', 'B40', 'B50', 'FC', 'BONE', 'LUNG'):
        if fam in s: return fam
    return s[:8]


def _bin_num(v, edges):
    x = try_float(v)
    if x is None: return '(missing)'
    for lo, hi in edges:
        if lo <= x < hi:
            return f'[{lo}, {hi})'
    return f'[{edges[-1][1]}+)'


def bin_tube_current(v): return _bin_num(v, [(0,200),(200,300),(300,400),(400,500),(500,10_000)])
def bin_exposure_uas(v):
    """Bin μAs (= ms×mA) into coarse ranges to avoid per-value bins."""
    return _bin_num(v, [(0,150_000),(150_000,250_000),(250_000,350_000),
                         (350_000,450_000),(450_000,1_000_000),(1_000_000,10_000_000)])
def bin_spacing(v):
    """Bin SpacingBetweenSlices in mm."""
    x = try_float(v)
    if x is None: return '(missing)'
    edges = [(0, 1.5), (1.5, 3.0), (3.0, 5.0), (5.0, 8.0), (8.0, 20.0)]
    for lo, hi in edges:
        if lo <= x < hi:
            return f'[{lo:.1f}, {hi:.1f})'
    return '[20+)'
def bin_recondiam(v):    return _bin_num(v, [(0,320),(320,360),(360,400),(400,500)])
def bin_tcw(v):          return _bin_num(v, [(0,20),(20,40),(40,80),(80,160)])
def bin_tablespeed(v):   return _bin_num(v, [(0,20),(20,40),(40,60),(60,100)])


def bin_pixspacing(v):
    x = parse_pixspacing_mean(v)
    if x is None: return '(missing)'
    edges = [(0.50, 0.65), (0.65, 0.75), (0.75, 0.85), (0.85, 1.00)]
    for lo, hi in edges:
        if lo <= x < hi:
            return f'[{lo:.2f}, {hi:.2f})'
    return '(out-of-range)'


def bin_pitch(v):
    x = try_float(v)
    if x is None: return '(missing)'
    return f'{round(x, 3)}'


def bin_revtime(v):
    x = try_float(v)
    if x is None: return '(missing)'
    return f'{round(x, 2)}'


def normalize_cat(v):
    if pd.isna(v): return '(missing)'
    s = str(v).strip().strip("'\"")
    return s if s and s.lower() != 'nan' else '(missing)'


BINNERS = {
    'cat':              normalize_cat,
    'cat_multival':     normalize_multival_cat,
    'kernel_family':    kernel_family,
    'num_tube_current': bin_tube_current,
    'num_recondiam':    bin_recondiam,
    'num_pixspacing':   bin_pixspacing,
    'num_pitch':        bin_pitch,
    'num_tcw':          bin_tcw,
    'num_tablespeed':   bin_tablespeed,
    'num_revtime':      bin_revtime,
    'num_exposure_uas': bin_exposure_uas,
    'num_spacing':      bin_spacing,
}


def bh_fdr(pvals):
    """Benjamini-Hochberg FDR-adjusted p-values."""
    pvals = np.asarray(pvals, dtype=float)
    mask = ~np.isnan(pvals)
    p = pvals[mask]
    if len(p) == 0:
        return np.full_like(pvals, np.nan)
    order = np.argsort(p)
    ranked = p[order]
    n = len(p)
    adj = ranked * n / (np.arange(n) + 1)
    # monotone from right
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out = np.full(n, np.nan)
    out[order] = adj
    fdr = np.full_like(pvals, np.nan)
    fdr[mask] = out
    return fdr


def classify(rows, fill_rate, n_testable_bins):
    if fill_rate < 0.05:
        return 'fully_missing'
    # count unique non-missing bins
    present = [r for r in rows if r['bin'] != '(missing)' and r['bin'] != '(column missing)']
    unique_bins = len({r['bin'] for r in present})
    if unique_bins <= 1:
        return 'invariant'
    if n_testable_bins < 2:
        return 'insufficient'
    return 'testable'


def main():
    df = pd.read_csv(IN_CSV)
    pu = df[df['status'] == 'ok'].copy()
    pu['prob'] = pu['prob'].astype(float)
    pu['pred'] = (pu['prob'] >= THRESHOLD).astype(int)
    N = len(pu)
    FP = int((pu['pred'] == 1).sum())
    TN = N - FP
    base_spec = TN / N
    base_lo, base_hi = wilson(base_spec, N)

    rows_all = []
    sections = {}
    cls_rows = []

    for col, hex_tag, cls, cls_name, bin_type, mentioned in ESSENTIAL:
        binner = BINNERS.get(bin_type, normalize_cat)
        vals = pu[col].apply(binner) if col in pu.columns else pd.Series(['(column missing)']*N)
        pu_work = pu.assign(_bin=vals.values)

        # Per-bin aggregation (specificity = TN/n where TN = prob<THRESHOLD)
        grouped = (
            pu_work.groupby('_bin')
                   .agg(n=('pred', 'size'), fp=('pred', 'sum'),
                        mean_prob=('prob', 'mean'))
                   .reset_index()
        )
        grouped['tn'] = grouped['n'] - grouped['fp']
        grouped['spec'] = grouped['tn'] / grouped['n']

        # Sort bins
        def sort_key(s):
            m = re.match(r'\[([-\d.]+),', s)
            if m: return (0, float(m.group(1)))
            try: return (1, float(s))
            except (ValueError, TypeError): return (2, s)
        grouped = grouped.assign(_sk=grouped['_bin'].map(sort_key)).sort_values('_sk').drop(columns='_sk')

        # Collect rows
        tag_rows = []
        for _, r in grouped.iterrows():
            n = int(r['n']); fp = int(r['fp']); tn = int(r['tn'])
            spec = tn / n if n else 0.0
            lo, hi = wilson(spec, n)
            tag_rows.append({
                'bin': r['_bin'], 'n': n, 'tn': tn, 'fp': fp,
                'spec': spec, 'ci_lo': lo, 'ci_hi': hi,
                'mean_prob': float(r['mean_prob']),
            })
            rows_all.append({
                'class': cls_name, 'tag': col, 'hex': hex_tag,
                'mentioned': mentioned, 'bin': r['_bin'],
                'n': n, 'tn': tn, 'fp': fp, 'spec': spec,
                'ci_lo': lo, 'ci_hi': hi, 'mean_prob': float(r['mean_prob']),
            })
        sections.setdefault(cls_name, []).append((col, hex_tag, bin_type, mentioned, tag_rows))

        # Fill rate, testable-bin count
        present = [r for r in tag_rows if r['bin'] not in ('(missing)', '(column missing)')]
        present_n = sum(r['n'] for r in present)
        fill = present_n / N if N else 0.0
        testable_bins = [r for r in present if r['n'] >= 5]
        n_testable = len(testable_bins)

        # Kruskal-Wallis on probs across testable bins
        kw_H = kw_p = np.nan
        if n_testable >= 2:
            groups = []
            for r in testable_bins:
                groups.append(pu_work[pu_work['_bin'] == r['bin']]['prob'].values)
            try:
                kw_H, kw_p = kruskal(*groups)
            except Exception:
                kw_H, kw_p = np.nan, np.nan

        cls_rows.append({
            'tag': col, 'hex': hex_tag, 'class': cls_name,
            'mentioned_in_paper': mentioned,
            'fill_rate': fill,
            'n_unique_bins': len({r['bin'] for r in present}),
            'n_testable_bins': n_testable,
            'status': classify(tag_rows, fill, n_testable),
            'kw_H': kw_H, 'kw_p': kw_p,
            'mean_min': min((r['mean_prob'] for r in testable_bins), default=np.nan),
            'mean_max': max((r['mean_prob'] for r in testable_bins), default=np.nan),
            'spec_min': min((r['spec'] for r in testable_bins), default=np.nan),
            'spec_max': max((r['spec'] for r in testable_bins), default=np.nan),
        })

    cls_df = pd.DataFrame(cls_rows)
    # BH-FDR across testable tags' kw_p
    cls_df['kw_p_fdr'] = bh_fdr(cls_df['kw_p'].values)
    cls_df['significant_fdr'] = (cls_df['kw_p_fdr'] < 0.05).fillna(False)
    cls_df['mean_delta'] = cls_df['mean_max'] - cls_df['mean_min']

    cls_df.to_csv(OUT_CLS, index=False)
    pd.DataFrame(rows_all).to_csv(OUT_CSV, index=False)

    with open(OUT_MD, 'w') as f:
        f.write('# Per-essential-tag SPECIFICITY + probability heterogeneity (CRLM, ST3DCN)\n\n')
        f.write(f'**Cohort**: {N} CRLM patients, 100% HCC-negative (specificity-only, no HCC ground truth).\n\n')
        f.write(f'**Threshold**: {THRESHOLD}  ·  **Binary pred**: `prob ≥ {THRESHOLD}` → pred HCC (FP here).\n\n')
        f.write(f'**Baseline pu**: n={N}, TN={TN}, FP={FP}, '
                f'specificity = **{base_spec:.4f}**, Wilson 95% CI = ({base_lo:.4f}, {base_hi:.4f}).\n\n')
        f.write(f'Primary analysis: Kruskal-Wallis on probabilities across bins, BH-FDR across testable tags '
                f'(same framework as HCC-TACE-Seg RESULTS_v5).\n\n')

        # Summary of classification
        f.write('## Essential-26 classification\n\n')
        counts = cls_df['status'].value_counts().to_dict()
        n_sig = int(cls_df['significant_fdr'].sum())
        f.write(f'- fully_missing: {counts.get("fully_missing", 0)}\n')
        f.write(f'- invariant:     {counts.get("invariant", 0)}\n')
        f.write(f'- insufficient:  {counts.get("insufficient", 0)}\n')
        f.write(f'- testable:      {counts.get("testable", 0)}\n')
        f.write(f'  - significant (FDR<0.05): {n_sig}\n')
        f.write(f'  - not significant:                {counts.get("testable", 0) - n_sig}\n\n')

        # Significant tags table
        sig_df = cls_df[cls_df['significant_fdr']].sort_values('kw_p')
        if len(sig_df):
            f.write('## Significant tags (FDR<0.05)\n\n')
            f.write('| Tag | Hex | Class | Mentioned | KW H | raw p | FDR p | mean range | spec range |\n')
            f.write('|---|---|---|---|---|---|---|---|---|\n')
            for _, r in sig_df.iterrows():
                f.write(f"| **{r['tag']}** | {r['hex']} | {r['class']} | "
                        f"{'✓' if r['mentioned_in_paper'] else '✗'} | "
                        f"{r['kw_H']:.2f} | {r['kw_p']:.4f} | **{r['kw_p_fdr']:.4f}** | "
                        f"{r['mean_min']:.3f}–{r['mean_max']:.3f} | "
                        f"{r['spec_min']:.3f}–{r['spec_max']:.3f} |\n")
            f.write('\n')

        # Tag availability overview
        f.write('## Tag availability overview\n\n')
        f.write('| Tag | Hex | Class | Fill% | # bins | Spec range | Mean prob range | Status |\n')
        f.write('|---|---|---|---|---|---|---|---|\n')
        for _, r in cls_df.iterrows():
            sr = (f"{r['spec_min']:.3f}–{r['spec_max']:.3f}"
                  if not pd.isna(r['spec_min']) else '—')
            mr = (f"{r['mean_min']:.3f}–{r['mean_max']:.3f}"
                  if not pd.isna(r['mean_min']) else '—')
            f.write(f"| {r['tag']} | {r['hex']} | {r['class'][0]} | "
                    f"{r['fill_rate']*100:.1f} | {r['n_unique_bins']} | "
                    f"{sr} | {mr} | {r['status']} |\n")
        f.write('\n*Class*: **i** = intensity/exposure, **g** = geometry, **d** = device, **a** = anatomy/position.\n\n')

        # Per-class detailed tables
        for cls_order in ['Intensity/exposure', 'Geometry', 'Device', 'Anatomy/position']:
            if cls_order not in sections:
                continue
            f.write(f'## {cls_order}\n\n')
            for col, hex_tag, bin_type, mentioned, tag_rows in sections[cls_order]:
                present = [r for r in tag_rows if r['bin'] not in ('(missing)', '(column missing)')]
                if not present:
                    f.write(f'- **{col}** (`{hex_tag}`): fully missing in CRLM (fill=0%).\n\n')
                    continue
                if len(present) == 1:
                    r = present[0]
                    missing_r = next((rr for rr in tag_rows if rr['bin'] in ('(missing)', '(column missing)')), None)
                    mn = missing_r['n'] if missing_r else 0
                    f.write(f'- **{col}** (`{hex_tag}`): single value `{r["bin"]}` '
                            f'(n={r["n"]}, spec={r["spec"]:.3f}, mean p={r["mean_prob"]:.3f})'
                            f"{f'; missing in {mn}/{N}' if mn else ' — invariant'}.\n\n")
                    continue
                f.write(f'### {col} (`{hex_tag}`)\n\n')
                f.write('| Value / bin | n | TN | FP | Spec | Wilson 95% CI | mean p |\n')
                f.write('|---|---|---|---|---|---|---|\n')
                for r in tag_rows:
                    if r['bin'] in ('(missing)', '(column missing)') and r['n'] == 0:
                        continue
                    f.write(f'| `{r["bin"]}` | {r["n"]} | {r["tn"]} | {r["fp"]} | '
                            f'{r["spec"]:.4f} | ({r["ci_lo"]:.4f}, {r["ci_hi"]:.4f}) | '
                            f'{r["mean_prob"]:.3f} |\n')
                f.write('\n')

    print(f'markdown    → {OUT_MD}')
    print(f'per-bin csv → {OUT_CSV}')
    print(f'classif csv → {OUT_CLS}')
    print(f'\npu baseline: n={N} TN={TN} FP={FP} spec={base_spec:.4f}')


if __name__ == '__main__':
    main()
