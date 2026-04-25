"""Combine multi-cohort per-tag analysis into a cross-cohort thesis report
with observations 1/2/3 and per-tag bin-level side-by-side tables.

Call `build_report(cohort_configs, out_md)` where cohort_configs is a dict:
  {
    'hcc_tace_seg': {
        'cls_csv': '.../per_tag_classification_augmented.csv',
        'bin_csv': '.../per_tag_binned_augmented.csv',
        'label':   'HCC-TACE-Seg (HCC+)',
        'metric':  'sens',                 # or 'spec' for HCC−
        'n_lesions': 146,
    },
    'crlm': {...},
    ...
  }
"""
import pandas as pd

from .essential_tags import ESSENTIAL_26


def _fmt(s):  return ', '.join(sorted(s)) or '∅'


def build_report(cohort_configs: dict, out_md: str, thesis_title: str = None):
    """Produce cross-cohort essential-vs-mentioned markdown report.

    cohort_configs keys become cohort identifiers in the report.
    """
    cohorts = {}
    for key, cfg in cohort_configs.items():
        cls = pd.read_csv(cfg['cls_csv']).set_index('tag')
        bins = pd.read_csv(cfg['bin_csv'])
        cohorts[key] = {'cls': cls, 'bins': bins, **cfg}

    mentioned = set(t for t,_,_,_,_,m in ESSENTIAL_26 if m)

    # Global metrics
    sig_sets = {k: set(v['cls'][v['cls']['significant_fdr']].index) for k, v in cohorts.items()}
    union_sig = set().union(*sig_sets.values())
    inter_sig = set.intersection(*sig_sets.values()) if len(sig_sets) > 1 else set()

    with open(out_md, 'w') as f:
        f.write(f'# {thesis_title or "Essential vs Mentioned DICOM Tags — cross-cohort"}\n\n')
        f.write('Cohorts:\n')
        for k, v in cohorts.items():
            f.write(f'- **{k}** — {v["label"]}, n_lesions={v.get("n_lesions", "?")}, metric={v["metric"]}\n')
        f.write(f'\nMentioned-in-paper subset: 11 / 26 essential tags.\n\n')

        # Section 1 — Observation 1
        f.write('## 1. Observation 1 — Mentioned-11 captures only a fraction of sig tags\n\n')
        f.write('| Cohort | FDR-sig total | mentioned | not-mentioned |\n')
        f.write('|---|---|---|---|\n')
        for k, v in cohorts.items():
            sig = sig_sets[k]; m = sig & mentioned; n = sig - mentioned
            f.write(f'| {k} | {len(sig)} | {_fmt(m)} ({len(m)}/{len(sig)}) | {_fmt(n)} |\n')
        f.write('\n')

        # Section 2 — Observation 2
        f.write('## 2. Observation 2 — Cross-cohort intersection\n\n')
        f.write(f'Union of FDR-sig tags: {{ {_fmt(union_sig)} }} (n={len(union_sig)})\n\n')
        f.write(f'Intersection: {{ {_fmt(inter_sig)} }}\n\n')
        for t in sorted(inter_sig):
            is_m = 'mentioned' if t in mentioned else '**not mentioned**'
            f.write(f'- `{t}`: {is_m}\n')

        # Section 3 — Observation 3
        f.write('\n## 3. Observation 3 — Mentioned tags often unusable\n\n')
        f.write('| Cohort | Mentioned testable | Mentioned sig | Mentioned unusable |\n')
        f.write('|---|---|---|---|\n')
        for k, v in cohorts.items():
            d = v['cls']
            mt = ((d['status']=='testable') & (d['mentioned_in_paper'])).sum()
            ms = ((d['significant_fdr']) & (d['mentioned_in_paper'])).sum()
            mu = ((d['status'].isin(['fully_missing','invariant','insufficient'])) & (d['mentioned_in_paper'])).sum()
            f.write(f'| {k} | {mt} | {ms} | {mu} |\n')

        # Section 4 — Full table
        f.write('\n## 4. Essential-26 classification — all cohorts side-by-side\n\n')
        hdr = '| Tag | Hex | Class | Mentioned |'
        sep = '|---|---|---|---|'
        for k in cohorts:
            hdr += f' {k} fill% | {k} status | {k} FDR |'
            sep += '---|---|---|'
        f.write(hdr + '\n' + sep + '\n')
        for col, hexa, ccode, cname, _, mentioned_flag in ESSENTIAL_26:
            row = f'| {col} | {hexa} | {cname[0]} | {"✓" if mentioned_flag else "✗"} |'
            for k, v in cohorts.items():
                cls = v['cls']
                if col in cls.index:
                    r = cls.loc[col]
                    fdr = f"{r['kw_p_fdr']:.4f}" if pd.notna(r['kw_p_fdr']) else '—'
                    sig = '★' if r['significant_fdr'] else ''
                    dr  = '◆' if str(r.get('derived_applied','')).strip() not in ('','nan') else ''
                    row += f' {r["fill_rate"]*100:.1f} | {r["status"]} | {fdr}{sig}{dr} |'
                else:
                    row += ' — | — | — |'
            f.write(row + '\n')
        f.write('\n★ = FDR-sig; ◆ = derivation-augmented.\n')

        # Section 5 — per-tag bin tables
        f.write('\n## 5. Per-tag bin tables (all cohorts side-by-side)\n\n')
        cur_cls_name = None
        for col, hexa, ccode, cname, _, mentioned_flag in ESSENTIAL_26:
            if cname != cur_cls_name:
                f.write(f'\n### ▸ {cname}\n')
                cur_cls_name = cname
            f.write(f'\n#### {col} `({hexa})` — {"✓ mentioned" if mentioned_flag else "✗ not mentioned"}\n\n')
            # status line
            for k, v in cohorts.items():
                cls = v['cls']
                if col in cls.index:
                    r = cls.loc[col]
                    fdr = f"FDR {r['kw_p_fdr']:.4f}" if pd.notna(r['kw_p_fdr']) else 'FDR —'
                    marker = '**★ sig**' if r['significant_fdr'] else r['status']
                    dr = ' ◆' if str(r.get('derived_applied','')).strip() not in ('','nan') else ''
                    f.write(f'- **{k}**: fill {r["fill_rate"]*100:.1f}% · {marker} · {fdr}{dr}\n')
            # combined bin table
            all_bins = []
            per_cohort = {}
            for k, v in cohorts.items():
                sub = v['bins'][v['bins']['tag'] == col].copy()
                if len(sub): all_bins.extend(sub['bin'].tolist())
                per_cohort[k] = {row['bin']: row for _, row in sub.iterrows()}
            seen = set(); all_bins = [b for b in all_bins if not (b in seen or seen.add(b))]
            if all_bins:
                hdr = '| bin |'
                sep = '|---|'
                for k, v in cohorts.items():
                    metric = v['metric']
                    hdr += f' {k} n | {k} {metric} | {k} mean p |'
                    sep += '---|---|---|'
                f.write(hdr + '\n' + sep + '\n')
                for b in all_bins:
                    row = f'| `{b}` |'
                    for k, v in cohorts.items():
                        r = per_cohort[k].get(b)
                        if r is None:
                            row += ' — | — | — |'
                        else:
                            met_col = v['metric']
                            metric_val = r.get(met_col)
                            row += (f' {int(r["n"])} | '
                                    f'{f"{metric_val:.3f}" if pd.notna(metric_val) else "—"} | '
                                    f'{f"{r["mean_prob"]:.3f}" if pd.notna(r["mean_prob"]) else "—"} |')
                    f.write(row + '\n')
    print(f'→ {out_md}')
