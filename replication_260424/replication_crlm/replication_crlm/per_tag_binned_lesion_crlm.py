"""Lesion-level per-essential-tag SPECIFICITY + probability heterogeneity on CRLM.

Same framework as per_tag_binned_crlm.py but unit is lesion (multi per patient).
"""
import importlib.util
spec = importlib.util.spec_from_file_location(
    'per_tag_binned_crlm',
    '/home/ubuntu/hcc_workspace/replication_crlm/per_tag_binned_crlm.py')
ptb = importlib.util.module_from_spec(spec); spec.loader.exec_module(ptb)

import re
import numpy as np
import pandas as pd
from scipy.stats import kruskal

IN_CSV = '/home/ubuntu/hcc_workspace/replication_crlm/lesion_results_with_tags.csv'
OUT_MD = '/home/ubuntu/hcc_workspace/replication_crlm/per_tag_binned_lesion.md'
OUT_CSV = '/home/ubuntu/hcc_workspace/replication_crlm/per_tag_binned_lesion.csv'
OUT_CLS = '/home/ubuntu/hcc_workspace/replication_crlm/per_tag_classification_lesion.csv'
THRESHOLD = 0.8


def main():
    df = pd.read_csv(IN_CSV)
    les = df[df['status'] == 'ok'].copy()
    les['prob'] = les['prob'].astype(float)
    les['pred'] = (les['prob'] >= THRESHOLD).astype(int)
    N = len(les)
    FP = int(les['pred'].sum())
    TN = N - FP
    base_spec = TN / N
    base_lo, base_hi = ptb.wilson(base_spec, N)

    rows_all = []
    sections = {}
    cls_rows = []

    for col, hex_tag, cls, cls_name, bin_type, mentioned in ptb.ESSENTIAL:
        binner = ptb.BINNERS.get(bin_type, ptb.normalize_cat)
        vals = les[col].apply(binner) if col in les.columns else pd.Series(['(column missing)']*N)
        les_work = les.assign(_bin=vals.values)

        grouped = (
            les_work.groupby('_bin')
                    .agg(n=('pred', 'size'), fp=('pred', 'sum'),
                         mean_prob=('prob', 'mean'), median_prob=('prob', 'median'))
                    .reset_index()
        )
        grouped['tn'] = grouped['n'] - grouped['fp']
        grouped['spec'] = grouped['tn'] / grouped['n']

        def sort_key(s):
            m = re.match(r'\[([-\d.]+),', s)
            if m: return (0, float(m.group(1)))
            try: return (1, float(s))
            except (ValueError, TypeError): return (2, s)
        grouped = grouped.assign(_sk=grouped['_bin'].map(sort_key)).sort_values('_sk').drop(columns='_sk')

        tag_rows = []
        for _, r in grouped.iterrows():
            n = int(r['n']); fp = int(r['fp']); tn = int(r['tn'])
            spec = tn / n if n else 0.0
            lo, hi = ptb.wilson(spec, n)
            tag_rows.append({
                'bin': r['_bin'], 'n': n, 'tn': tn, 'fp': fp,
                'spec': spec, 'ci_lo': lo, 'ci_hi': hi,
                'mean_prob': float(r['mean_prob']),
                'median_prob': float(r['median_prob']),
            })
            rows_all.append({
                'class': cls_name, 'tag': col, 'hex': hex_tag,
                'mentioned': mentioned, 'bin': r['_bin'],
                'n': n, 'tn': tn, 'fp': fp, 'spec': spec,
                'ci_lo': lo, 'ci_hi': hi, 'mean_prob': float(r['mean_prob']),
            })
        sections.setdefault(cls_name, []).append((col, hex_tag, bin_type, mentioned, tag_rows))

        present = [r for r in tag_rows if r['bin'] not in ('(missing)', '(column missing)')]
        present_n = sum(r['n'] for r in present)
        fill = present_n / N if N else 0.0
        testable_bins = [r for r in present if r['n'] >= 5]
        n_testable = len(testable_bins)

        kw_H = kw_p = np.nan
        if n_testable >= 2:
            groups = []
            for r in testable_bins:
                groups.append(les_work[les_work['_bin'] == r['bin']]['prob'].values)
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
            'status': ptb.classify(tag_rows, fill, n_testable),
            'kw_H': kw_H, 'kw_p': kw_p,
            'mean_min': min((r['mean_prob'] for r in testable_bins), default=np.nan),
            'mean_max': max((r['mean_prob'] for r in testable_bins), default=np.nan),
            'spec_min': min((r['spec'] for r in testable_bins), default=np.nan),
            'spec_max': max((r['spec'] for r in testable_bins), default=np.nan),
        })

    cls_df = pd.DataFrame(cls_rows)
    cls_df['kw_p_fdr'] = ptb.bh_fdr(cls_df['kw_p'].values)
    cls_df['significant_fdr'] = (cls_df['kw_p_fdr'] < 0.05).fillna(False)
    cls_df['mean_delta'] = cls_df['mean_max'] - cls_df['mean_min']

    cls_df.to_csv(OUT_CLS, index=False)
    pd.DataFrame(rows_all).to_csv(OUT_CSV, index=False)

    with open(OUT_MD, 'w') as f:
        f.write(f'# Per-essential-tag SPECIFICITY + probability heterogeneity '
                f'(CRLM lesion-level, ST3DCN)\n\n')
        n_pat = les['patient_id'].nunique()
        f.write(f'**Unit**: lesion (Tumor_N segments).  '
                f'**N lesions**: {N}  ·  **N patients**: {n_pat}.\n\n')
        f.write(f'**Threshold**: {THRESHOLD}  ·  **Binary**: `prob ≥ {THRESHOLD}` → pred HCC (FP here).\n\n')
        f.write(f'**Baseline pu**: n={N}, TN={TN}, FP={FP}, '
                f'specificity = **{base_spec:.4f}**, Wilson 95% CI = ({base_lo:.4f}, {base_hi:.4f}).\n\n')

        f.write('## Essential-26 classification (lesion)\n\n')
        counts = cls_df['status'].value_counts().to_dict()
        n_sig = int(cls_df['significant_fdr'].sum())
        f.write(f'- fully_missing: {counts.get("fully_missing", 0)}\n')
        f.write(f'- invariant:     {counts.get("invariant", 0)}\n')
        f.write(f'- insufficient:  {counts.get("insufficient", 0)}\n')
        f.write(f'- testable:      {counts.get("testable", 0)}\n')
        f.write(f'  - significant (FDR<0.05): {n_sig}\n')
        f.write(f'  - not significant:                {counts.get("testable", 0) - n_sig}\n\n')

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
        f.write('\n')

        for cls_order in ['Intensity/exposure', 'Geometry', 'Device', 'Anatomy/position']:
            if cls_order not in sections:
                continue
            f.write(f'## {cls_order}\n\n')
            for col, hex_tag, bin_type, mentioned, tag_rows in sections[cls_order]:
                present = [r for r in tag_rows if r['bin'] not in ('(missing)', '(column missing)')]
                if not present:
                    f.write(f'- **{col}** (`{hex_tag}`): fully missing (fill=0%).\n\n')
                    continue
                if len(present) == 1:
                    r = present[0]
                    f.write(f'- **{col}** (`{hex_tag}`): single value `{r["bin"]}` '
                            f'(n={r["n"]}, spec={r["spec"]:.3f}, mean p={r["mean_prob"]:.3f}).\n\n')
                    continue
                f.write(f'### {col} (`{hex_tag}`)\n\n')
                f.write('| Value / bin | n | TN | FP | Spec | Wilson 95% CI | mean p | median p |\n')
                f.write('|---|---|---|---|---|---|---|---|\n')
                for r in tag_rows:
                    if r['bin'] in ('(missing)', '(column missing)') and r['n'] == 0:
                        continue
                    f.write(f'| `{r["bin"]}` | {r["n"]} | {r["tn"]} | {r["fp"]} | '
                            f'{r["spec"]:.4f} | ({r["ci_lo"]:.4f}, {r["ci_hi"]:.4f}) | '
                            f'{r["mean_prob"]:.3f} | {r["median_prob"]:.3f} |\n')
                f.write('\n')

    print(f'markdown    → {OUT_MD}')
    print(f'per-bin csv → {OUT_CSV}')
    print(f'classif csv → {OUT_CLS}')
    print(f'\nlesion baseline: n={N} TN={TN} FP={FP} spec={base_spec:.4f}  '
          f'n_patients={les["patient_id"].nunique()}')


if __name__ == '__main__':
    main()
