"""Augmented per-tag binned analysis for both cohorts.

Uses `lesion_results_with_derived.csv`. Where a direct field is fully-missing
but a `_derived` counterpart is populated, substitute the derived value and
LABEL the tag's status as `derivable → <new_status>` in the classification output.

Outputs (per cohort): per_tag_classification_augmented.csv, per_tag_binned_augmented.md
"""
import sys
import re
import importlib.util
import numpy as np
import pandas as pd
from scipy.stats import kruskal

# Reuse CRLM binners, wilson, bh_fdr (same logic as HCC)
spec = importlib.util.spec_from_file_location(
    'ptb_crlm', '/home/ubuntu/hcc_workspace/replication_crlm/per_tag_binned_crlm.py')
ptb = importlib.util.module_from_spec(spec); spec.loader.exec_module(ptb)


def run(cohort_key, in_csv, out_md, out_csv, out_cls, all_positive_label=True):
    """all_positive_label: HCC → sensitivity; CRLM → specificity."""
    df = pd.read_csv(in_csv)
    pu = df[df['status'] == 'ok'].copy()
    pu['prob'] = pu['prob'].astype(float)
    THR = 0.8
    pu['pred'] = (pu['prob'] >= THR).astype(int)
    N = len(pu)

    # Derivation substitution: for these tags, use derived column if present
    # and collapse non-missing direct values (both should agree; direct wins if both exist)
    DERIVATION_COLUMNS = {
        'SpacingBetweenSlices': 'SpacingBetweenSlices_derived',
        'Exposure_uAs':         'Exposure_uAs_derived',
        'Manufacturer':         'Manufacturer_derived',           # HCC only meaningful
        'ManufacturerModelName':'ManufacturerModelName_derived',
    }
    derived_applied = {}
    for direct, derived in DERIVATION_COLUMNS.items():
        if derived in pu.columns:
            direct_vals = pu.get(direct)
            direct_fill = 0
            if direct_vals is not None:
                direct_fill = (direct_vals.notna() & (direct_vals.astype(str).str.strip() != '')
                                & (direct_vals.astype(str).str.lower() != 'nan')).sum()
            derived_fill = pu[derived].notna().sum()
            # Use derived when direct is missing for that row, or use derived entirely
            # when direct is <5% filled (per our "derivable" rule)
            if derived_fill > direct_fill:
                pu[direct] = pu[direct].where(pu[direct].notna() & (pu[direct].astype(str).str.strip() != ''),
                                               pu[derived])
                # But if direct was mostly missing, use derived outright
                if direct_fill / max(1, N) < 0.05 and derived_fill > 0:
                    pu[direct] = pu[derived]
                    derived_applied[direct] = f'substituted from {derived} (direct fill {direct_fill}/{N})'
                else:
                    derived_applied[direct] = f'filled-in from {derived} where direct missing'

    # ===== per-tag classification =====
    rows_all = []; sections = {}; cls_rows = []
    for col, hex_tag, cls, cls_name, bin_type, mentioned in ptb.ESSENTIAL:
        binner = ptb.BINNERS.get(bin_type, ptb.normalize_cat)
        vals = pu[col].apply(binner) if col in pu.columns else pd.Series(['(column missing)']*N)
        work = pu.assign(_bin=vals.values)

        grouped = (
            work.groupby('_bin')
                .agg(n=('pred', 'size'), fp_or_tp=('pred', 'sum'),
                     mean_prob=('prob', 'mean'), median_prob=('prob', 'median'))
                .reset_index()
        )

        def sort_key(s):
            m = re.match(r'\[([-\d.]+),', s)
            if m: return (0, float(m.group(1)))
            try: return (1, float(s))
            except (ValueError, TypeError): return (2, s)
        grouped = grouped.assign(_sk=grouped['_bin'].map(sort_key)).sort_values('_sk').drop(columns='_sk')

        tag_rows = []
        for _, r in grouped.iterrows():
            n = int(r['n']); hits = int(r['fp_or_tp'])
            if all_positive_label:
                sens = hits / n if n else 0.0
                lo, hi = ptb.wilson(sens, n)
                rec = {'bin': r['_bin'], 'n': n, 'tp': hits, 'sens': sens,
                       'mean_prob': float(r['mean_prob']), 'median_prob': float(r['median_prob']),
                       'ci_lo': lo, 'ci_hi': hi}
            else:
                fp = hits; tn = n - fp
                spec_ = tn / n if n else 0.0
                lo, hi = ptb.wilson(spec_, n)
                rec = {'bin': r['_bin'], 'n': n, 'tn': tn, 'fp': fp, 'spec': spec_,
                       'mean_prob': float(r['mean_prob']), 'median_prob': float(r['median_prob']),
                       'ci_lo': lo, 'ci_hi': hi}
            tag_rows.append(rec)
            rows_all.append({'class': cls_name, 'tag': col, 'hex': hex_tag,
                             'mentioned': mentioned,
                             'derived_applied': derived_applied.get(col, ''),
                             **rec})
        sections.setdefault(cls_name, []).append((col, hex_tag, bin_type, mentioned, tag_rows))

        present = [r for r in tag_rows if r['bin'] not in ('(missing)', '(column missing)')]
        present_n = sum(r['n'] for r in present)
        fill = present_n / N if N else 0.0
        testable_bins = [r for r in present if r['n'] >= 5]
        n_testable = len(testable_bins)

        kw_H = kw_p = np.nan
        if n_testable >= 2:
            groups = [work[work['_bin'] == r['bin']]['prob'].values for r in testable_bins]
            try:
                kw_H, kw_p = kruskal(*groups)
            except Exception:
                pass

        status = ptb.classify(tag_rows, fill, n_testable)
        cls_rows.append({
            'tag': col, 'hex': hex_tag, 'class': cls_name,
            'mentioned_in_paper': mentioned,
            'derived_applied': derived_applied.get(col, ''),
            'fill_rate': fill,
            'n_unique_bins': len({r['bin'] for r in present}),
            'n_testable_bins': n_testable,
            'status': status,
            'kw_H': kw_H, 'kw_p': kw_p,
            'mean_min': min((r['mean_prob'] for r in testable_bins), default=np.nan),
            'mean_max': max((r['mean_prob'] for r in testable_bins), default=np.nan),
        })

    cls_df = pd.DataFrame(cls_rows)
    cls_df['kw_p_fdr'] = ptb.bh_fdr(cls_df['kw_p'].values)
    cls_df['significant_fdr'] = (cls_df['kw_p_fdr'] < 0.05).fillna(False)
    cls_df.to_csv(out_cls, index=False)
    pd.DataFrame(rows_all).to_csv(out_csv, index=False)

    # ===== write markdown =====
    with open(out_md, 'w') as f:
        label = 'sensitivity' if all_positive_label else 'specificity'
        f.write(f'# Augmented per-tag {label} + probability heterogeneity '
                f'(cohort={cohort_key.upper()}, ST3DCN)\n\n')
        f.write('Derivation substitutions applied:\n')
        for t, note in derived_applied.items():
            f.write(f'- **{t}**: {note}\n')
        if not derived_applied:
            f.write('- (no direct fields were <5% AND had non-zero derived counterpart)\n')
        f.write('\n')

        f.write('## Essential-26 classification (augmented)\n\n')
        counts = cls_df['status'].value_counts().to_dict()
        n_sig = int(cls_df['significant_fdr'].sum())
        for k in ['fully_missing','invariant','insufficient','testable']:
            f.write(f'- {k}: {counts.get(k, 0)}\n')
        f.write(f'  - significant (FDR<0.05): {n_sig}\n\n')

        sig_df = cls_df[cls_df['significant_fdr']].sort_values('kw_p')
        if len(sig_df):
            f.write('## Significant tags (FDR<0.05)\n\n')
            f.write('| Tag | Hex | Class | Mentioned | Derived? | KW H | raw p | FDR p | mean range |\n')
            f.write('|---|---|---|---|---|---|---|---|---|\n')
            for _, r in sig_df.iterrows():
                dr = '✓' if r['derived_applied'] else ''
                f.write(f"| **{r['tag']}** | {r['hex']} | {r['class']} | "
                        f"{'✓' if r['mentioned_in_paper'] else '✗'} | {dr} | "
                        f"{r['kw_H']:.2f} | {r['kw_p']:.4f} | **{r['kw_p_fdr']:.4f}** | "
                        f"{r['mean_min']:.3f}–{r['mean_max']:.3f} |\n")
            f.write('\n')

        f.write('## Tag classification after derivation\n\n')
        f.write('| Tag | Hex | Mentioned | Derived? | Fill% | status | KW FDR |\n')
        f.write('|---|---|---|---|---|---|---|\n')
        for _, r in cls_df.iterrows():
            dr = '✓' if r['derived_applied'] else ''
            fdr = f"{r['kw_p_fdr']:.4f}" if not pd.isna(r['kw_p_fdr']) else '—'
            f.write(f"| {r['tag']} | {r['hex']} | "
                    f"{'✓' if r['mentioned_in_paper'] else '✗'} | {dr} | "
                    f"{r['fill_rate']*100:.1f} | {r['status']} | {fdr} |\n")

    print(f'[{cohort_key}] → {out_md}')
    print(f'[{cohort_key}] → {out_cls}')


def main():
    # HCC: sensitivity
    run('hcc',
        '/home/ubuntu/hcc_workspace/replication/lesion_results_with_derived.csv',
        '/home/ubuntu/hcc_workspace/replication/per_tag_binned_augmented.md',
        '/home/ubuntu/hcc_workspace/replication/per_tag_binned_augmented.csv',
        '/home/ubuntu/hcc_workspace/replication/per_tag_classification_augmented.csv',
        all_positive_label=True)
    # CRLM: specificity
    run('crlm',
        '/home/ubuntu/hcc_workspace/replication_crlm/lesion_results_with_derived.csv',
        '/home/ubuntu/hcc_workspace/replication_crlm/per_tag_binned_augmented.md',
        '/home/ubuntu/hcc_workspace/replication_crlm/per_tag_binned_augmented.csv',
        '/home/ubuntu/hcc_workspace/replication_crlm/per_tag_classification_augmented.csv',
        all_positive_label=False)


if __name__ == '__main__':
    main()
