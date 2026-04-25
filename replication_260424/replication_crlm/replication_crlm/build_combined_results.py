"""Build RESULTS_ESSENTIAL_vs_MENTIONED.md — unified report on essential-tag
coverage across HCC-TACE-Seg (sens) and CRLM (spec), direct + derivation.

Sections:
  0. Headline thesis numbers
  1. Observation 1 — Mentioned coverage
  2. Observation 2 — Cross-cohort intersection emerges after derivation
  3. Observation 3 — Mentioned-tag unusability
  4. Per-tag side-by-side bin tables (for every essential tag)
"""
import pandas as pd, numpy as np

OUT_MD = '/home/ubuntu/hcc_workspace/replication_crlm/RESULTS_ESSENTIAL_vs_MENTIONED.md'

HCC_CLS  = '/home/ubuntu/hcc_workspace/replication/per_tag_classification_augmented.csv'
CRL_CLS  = '/home/ubuntu/hcc_workspace/replication_crlm/per_tag_classification_augmented.csv'
HCC_BIN  = '/home/ubuntu/hcc_workspace/replication/per_tag_binned_augmented.csv'
CRL_BIN  = '/home/ubuntu/hcc_workspace/replication_crlm/per_tag_binned_augmented.csv'

ESSENTIAL_ORDER = [
    # (tag, hex, class, mentioned)
    ('KVP',                     '180060', 'Intensity/exposure', True),
    ('XRayTubeCurrent',         '181151', 'Intensity/exposure', True),
    ('ConvolutionKernel',       '181210', 'Intensity/exposure', False),
    ('ContrastBolusAgent',      '180010', 'Intensity/exposure', True),
    ('CTDIvol',                 '189345', 'Intensity/exposure', False),
    ('ContrastFlowRate',        '181046', 'Intensity/exposure', True),
    ('Exposure_uAs',            '181153', 'Intensity/exposure', False),
    ('RevolutionTime',          '189305', 'Intensity/exposure', False),
    ('ReconstructionAlgorithm', '189315', 'Intensity/exposure', False),
    ('ContrastBolusAgentPhase', '189344', 'Intensity/exposure', True),
    ('SliceThickness',          '180050', 'Geometry',           True),
    ('Rows',                    '280010', 'Geometry',           False),
    ('Columns',                 '280011', 'Geometry',           False),
    ('PixelSpacing',            '280030', 'Geometry',           False),
    ('ReconstructionDiameter',  '181100', 'Geometry',           False),
    ('WindowCenter',            '281050', 'Geometry',           False),
    ('WindowWidth',             '281051', 'Geometry',           False),
    ('TotalCollimationWidth',   '189307', 'Geometry',           False),
    ('SpiralPitchFactor',       '189311', 'Geometry',           True),
    ('TableSpeed',              '189309', 'Geometry',           False),
    ('SpacingBetweenSlices',    '180088', 'Geometry',           False),
    ('ReformattingThickness',   '720512', 'Geometry',           True),
    ('ManufacturerModelName',   '081090', 'Device',             True),
    ('Manufacturer',            '080070', 'Device',             True),
    ('BodyPartExamined',        '180015', 'Anatomy/position',   True),
    ('PatientPosition',         '185100', 'Anatomy/position',   False),
]


def load_all():
    hcc_cls = pd.read_csv(HCC_CLS).set_index('tag')
    crl_cls = pd.read_csv(CRL_CLS).set_index('tag')
    hcc_bin = pd.read_csv(HCC_BIN)
    crl_bin = pd.read_csv(CRL_BIN)
    return hcc_cls, crl_cls, hcc_bin, crl_bin


def cls_cell(cls_df, tag):
    if tag not in cls_df.index:
        return ('—', '—', '—', False, '')
    r = cls_df.loc[tag]
    fill = f"{r['fill_rate']*100:.1f}"
    status = r['status']
    fdr = f"{r['kw_p_fdr']:.4f}" if not pd.isna(r['kw_p_fdr']) else '—'
    sig = bool(r['significant_fdr'])
    derived = str(r.get('derived_applied','')).strip()
    derived_flag = '◆' if derived and derived.lower() != 'nan' else ''
    return (fill, status, fdr, sig, derived_flag)


def build_overview_table(hcc_cls, crl_cls):
    lines = []
    lines.append('| Tag | Hex | Class | Mentioned | HCC fill% | HCC status | HCC FDR | HCC | CRLM fill% | CRLM status | CRLM FDR | CRLM |')
    lines.append('|---|---|---|---|---|---|---|---|---|---|---|---|')
    for tag, hexa, cls, ment in ESSENTIAL_ORDER:
        m = '✓' if ment else '✗'
        hfill, hst, hfdr, hsig, hd = cls_cell(hcc_cls, tag)
        cfill, cst, cfdr, csig, cd = cls_cell(crl_cls, tag)
        hmark = ('**★**' if hsig else '') + hd
        cmark = ('**★**' if csig else '') + cd
        lines.append(f"| {tag} | {hexa} | {cls[0]} | {m} | {hfill} | {hst} | {hfdr} | {hmark} | {cfill} | {cst} | {cfdr} | {cmark} |")
    lines.append('')
    lines.append('★ = FDR-significant (<0.05) ; ◆ = derivation-augmented (SpacingBetweenSlices ← IPP; Exposure_uAs ← ExposureTime × XRayTubeCurrent)')
    return '\n'.join(lines)


def bin_rows(df, tag):
    """Filter per-bin rows for a tag, stable order."""
    sub = df[df['tag'] == tag].copy()
    if len(sub) == 0:
        return sub
    # Sort numeric bins first, then categorical
    import re
    def k(s):
        s = str(s)
        m = re.match(r'\[([-\d.]+),', s)
        if m: return (0, float(m.group(1)))
        try: return (1, float(s))
        except (ValueError, TypeError): return (2, s)
    sub['_sk'] = sub['bin'].map(k)
    return sub.sort_values('_sk')


def per_tag_bin_section(tag, hexa, cls, ment, hcc_bin, crl_bin, hcc_cls, crl_cls):
    lines = []
    m = '✓ mentioned' if ment else '✗ not mentioned'
    lines.append(f'\n### {tag} `({hexa})` — {cls} · {m}\n')
    # status line for each cohort
    hfill, hst, hfdr, hsig, hd = cls_cell(hcc_cls, tag)
    cfill, cst, cfdr, csig, cd = cls_cell(crl_cls, tag)
    hmarker = '**★ FDR-sig**' if hsig else ('testable' if hst=='testable' else hst)
    cmarker = '**★ FDR-sig**' if csig else ('testable' if cst=='testable' else cst)
    der_note = ''
    if hd or cd:
        der_note = ' — derivation-augmented (ExposureTime×XRayTubeCurrent / IPP z-diff)'
    lines.append(f'- **HCC**  · fill {hfill}% · status {hmarker} · FDR {hfdr}{der_note}')
    lines.append(f'- **CRLM** · fill {cfill}% · status {cmarker} · FDR {cfdr}')
    lines.append('')

    # combined table
    hsub = bin_rows(hcc_bin, tag)
    csub = bin_rows(crl_bin, tag)
    all_bins = list(dict.fromkeys(list(hsub['bin']) + list(csub['bin'])))
    if not all_bins:
        lines.append('(no bin data)')
        return '\n'.join(lines)

    lines.append('| bin | HCC n | HCC sens | HCC mean p | CRLM n | CRLM spec | CRLM mean p |')
    lines.append('|---|---|---|---|---|---|---|')
    hmap = {r['bin']: r for _, r in hsub.iterrows()}
    cmap = {r['bin']: r for _, r in csub.iterrows()}
    for b in all_bins:
        h = hmap.get(b); c = cmap.get(b)
        hn  = int(h['n']) if h is not None else 0
        hs  = f"{h['sens']:.3f}" if h is not None and 'sens' in h and not pd.isna(h['sens']) else '—'
        hm  = f"{h['mean_prob']:.3f}" if h is not None and not pd.isna(h['mean_prob']) else '—'
        cn  = int(c['n']) if c is not None else 0
        cs  = f"{c['spec']:.3f}" if c is not None and 'spec' in c and not pd.isna(c['spec']) else '—'
        cm  = f"{c['mean_prob']:.3f}" if c is not None and not pd.isna(c['mean_prob']) else '—'
        hn_s = str(hn) if hn > 0 else '—'
        cn_s = str(cn) if cn > 0 else '—'
        lines.append(f'| `{b}` | {hn_s} | {hs} | {hm} | {cn_s} | {cs} | {cm} |')
    return '\n'.join(lines)


def mentioned_sets(hcc_cls, crl_cls):
    """Compute observation-1/2/3 numbers."""
    hcc_sig = set(hcc_cls[hcc_cls['significant_fdr']].index)
    crl_sig = set(crl_cls[crl_cls['significant_fdr']].index)
    intersection = hcc_sig & crl_sig
    union = hcc_sig | crl_sig
    mentioned = set(t for t,_,_,m in ESSENTIAL_ORDER if m)
    return {
        'hcc_sig': hcc_sig,
        'crl_sig': crl_sig,
        'intersection': intersection,
        'union': union,
        'mentioned': mentioned,
    }


def main():
    hcc_cls, crl_cls, hcc_bin, crl_bin = load_all()
    sets = mentioned_sets(hcc_cls, crl_cls)

    hcc_sig = sets['hcc_sig']; crl_sig = sets['crl_sig']
    mentioned = sets['mentioned']
    hcc_m = hcc_sig & mentioned; hcc_n = hcc_sig - mentioned
    crl_m = crl_sig & mentioned; crl_n = crl_sig - mentioned
    union = sets['union']; inter = sets['intersection']
    union_m = union & mentioned; union_n = union - mentioned

    # Direct-only counts (from pre-derivation classification files)
    hcc_direct = pd.read_csv('/home/ubuntu/hcc_workspace/replication/essential_26_classification.csv').set_index('tag')
    crl_direct = pd.read_csv('/home/ubuntu/hcc_workspace/replication_crlm/per_tag_classification_lesion.csv').set_index('tag')
    hcc_direct_sig = set(hcc_direct[hcc_direct['significant_fdr']==True].index)
    crl_direct_sig = set(crl_direct[crl_direct['significant_fdr']==True].index)
    direct_union = hcc_direct_sig | crl_direct_sig
    direct_inter = hcc_direct_sig & crl_direct_sig
    direct_union_m = direct_union & mentioned

    def fmt(s): return ', '.join(sorted(s)) or '∅'

    with open(OUT_MD, 'w') as f:
        # Header
        f.write('# Essential vs Mentioned DICOM Tags — ST3DCN replication on HCC-TACE-Seg + CRLM\n\n')
        f.write('**Cohorts**:\n')
        f.write('- **HCC-TACE-Seg** — 105 HCC+ patients / 146 lesions (TCIA). ST3DCN sensitivity path.\n')
        f.write('- **CRLM** — 197 HCC− patients / 464 lesions (TCIA, Colorectal-Liver-Metastases). ST3DCN specificity path.\n\n')
        f.write('Same pipeline, same preprocessing (70³ crop, HU window [-160,240]), same threshold 0.8 on p(HCC). Lesion-level analysis in both. Per-essential-tag heterogeneity tested by Kruskal-Wallis on p(HCC), with **BH-FDR α = 0.05** across testable tags.\n\n')
        f.write('**Essential set**: 26 CT tags (protocol2dicom). **Mentioned set**: 11 of those 26 (in ST3DCN paper Methods / Table S1). Table S1 also mentions 3 non-CT attributes (AcquisitionSequence / InjectionDelay / AgentSequence) which are not value-typed DICOM fields; excluded from per-tag testing.\n\n')
        f.write('**Two DICOM-standard derivations applied** to rescue fields that are direct-populated <5%:\n')
        f.write('- `SpacingBetweenSlices (0018,0088)` ← `median(diff(sorted ImagePositionPatient[2]))` — DICOM PS3.3 §C.7.6.2.1.1\n')
        f.write('- `Exposure_uAs (0018,1153)` ← `ExposureTime × XRayTubeCurrent` — DICOM PS3.3 §C.8.7.2.1 (standard-listed example calc; chosen over `Exposure×1000` because of 0018,1152 vendor-inconsistent semantics, intra-scanner CV 0.74-1.31 vs 0.13-0.38)\n\n')
        f.write('Non-recoverable despite effort: Manufacturer / ModelName on HCC-TACE-Seg (NBIA digest itself has only 111/677 rows populated).\n\n')

        # Headline — all numbers dynamic
        f.write('## 0. Headline — thesis numbers\n\n')
        hcc_direct_m = hcc_direct_sig & mentioned; hcc_direct_n = hcc_direct_sig - mentioned
        crl_direct_m = crl_direct_sig & mentioned; crl_direct_n = crl_direct_sig - mentioned
        f.write('| | Direct-only | + DICOM-standard derivation |\n')
        f.write('|---|---|---|\n')
        f.write(f'| HCC sig tags | {len(hcc_direct_sig)} ({len(hcc_direct_m)} mentioned + {len(hcc_direct_n)} not) | **{len(hcc_sig)}** ({len(hcc_m)} mentioned + {len(hcc_n)} not) |\n')
        f.write(f'| CRLM sig tags | {len(crl_direct_sig)} ({len(crl_direct_m)} mentioned + {len(crl_direct_n)} not) | **{len(crl_sig)}** ({len(crl_m)} mentioned + {len(crl_n)} not) |\n')
        f.write(f'| Union sig tags (both cohorts) | {len(direct_union)} | **{len(union)}** |\n')
        f.write(f'| HCC ∩ CRLM sig tags | **{fmt(direct_inter)}** | **{{ {fmt(inter)} }}** |\n')
        f.write(f'| Mentioned-11 coverage of union | {len(direct_union_m)}/{len(direct_union)} = **{len(direct_union_m)/max(1,len(direct_union))*100:.0f}%** | {len(union_m)}/{len(union)} = **{len(union_m)/max(1,len(union))*100:.0f}%** |\n')
        f.write(f'| Essential-26 coverage of union | {len(direct_union)}/{len(direct_union)} = **100%** | {len(union)}/{len(union)} = **100%** |\n\n')
        f.write('Derivation widens the essential-vs-mentioned gap — newly-testable derived tags are all *not mentioned* in the paper.\n\n')

        # Observation 1 — dynamic
        f.write(f'## 1. Observation 1 — Mentioned-11 captures only {len(hcc_m)}/{len(hcc_sig)} (HCC) and {len(crl_m)}/{len(crl_sig)} (CRLM) of FDR-sig tags\n\n')
        f.write('| Cohort | FDR-sig total | mentioned in paper | NOT mentioned in paper |\n')
        f.write('|---|---|---|---|\n')
        f.write(f'| HCC-TACE-Seg (all HCC+) | {len(hcc_sig)} | {fmt(hcc_m)} (**{len(hcc_m)}/{len(hcc_sig)}**) | {fmt(hcc_n)} |\n')
        f.write(f'| CRLM (all HCC−)         | {len(crl_sig)} | {fmt(crl_m)} (**{len(crl_m)}/{len(crl_sig)}**) | {fmt(crl_n)} |\n\n')
        f.write(f'Across both cohorts, {len(union_m)} of {len(union)} FDR-sig tags are mentioned; the remaining {len(union_n)} ({fmt(union_n)}) are essential-set-only discoveries that the paper did not cite.\n\n')

        # Observation 2 — dynamic
        if len(inter):
            inter_txt = '{' + fmt(inter) + '}'
        else:
            inter_txt = '∅'
        f.write(f'## 2. Observation 2 — Cross-cohort intersection: {fmt(direct_inter)} (direct) → {inter_txt} (derived)\n\n')
        f.write(f'- **Direct-only**: HCC ∩ CRLM = {fmt(direct_inter)}\n')
        f.write(f'- **After derivation**: HCC ∩ CRLM = {inter_txt}\n')
        for t in sorted(inter):
            is_m = 'mentioned' if t in mentioned else '**not mentioned**'
            f.write(f'  - `{t}`: {is_m} in paper; essential-26 ∈\n')
        f.write('\n')
        f.write('Which tags express variability is deployment-specific. A single-cohort study (like the ST3DCN paper) cannot enumerate them a priori — it can only find the tags visible in its own distribution. Consequences:\n')
        f.write('1. The safety net needs to be wider than any single cohort\'s sig set → essential-26 exists for this reason.\n')
        f.write(f'2. Using just the paper\'s mentioned-11 on HCC-TACE-Seg would miss {len(hcc_n)} tags; on CRLM, {len(crl_n)} tags — all in essential-26 but not in Methods.\n')
        f.write('3. Derivation reveals tags that would be invisible in naive fill-only analysis. In particular, Exposure_uAs is direct-0% in both cohorts but derivable via the DICOM-sanctioned formula — turning a "fully missing" row into a testable (and in HCC, sig) one.\n\n')

        # Observation 3
        f.write('## 3. Observation 3 — Most mentioned tags are simply unusable\n\n')
        f.write('Not every mentioned tag is populated. The paper mentioning a tag in Methods does not guarantee the tag is observable on deployment DICOMs. *Even after derivation attempts, most mentioned tags remain non-testable.*\n\n')
        n_hcc_test = ((hcc_cls['status']=='testable') & (hcc_cls['mentioned_in_paper'])).sum()
        n_hcc_unus = ((hcc_cls['status'].isin(['fully_missing','invariant','insufficient'])) & (hcc_cls['mentioned_in_paper'])).sum()
        n_crl_test = ((crl_cls['status']=='testable') & (crl_cls['mentioned_in_paper'])).sum()
        n_crl_unus = ((crl_cls['status'].isin(['fully_missing','invariant','insufficient'])) & (crl_cls['mentioned_in_paper'])).sum()
        f.write(f'| Cohort | Mentioned-11 testable | Mentioned-11 FDR-sig | Mentioned-11 unusable |\n')
        f.write(f'|---|---|---|---|\n')
        f.write(f'| HCC-TACE-Seg | {n_hcc_test} | {len(hcc_m)} | **{n_hcc_unus}** |\n')
        f.write(f'| CRLM         | {n_crl_test} | {len(crl_m)} | **{n_crl_unus}** |\n\n')
        f.write('Unusable mentioned tags per cohort:\n\n')
        for lab, df in [('HCC', hcc_cls), ('CRLM', crl_cls)]:
            unusable = df[(df['status'].isin(['fully_missing','invariant','insufficient'])) & (df['mentioned_in_paper'])]
            items = [f"{t} ({r['status']})" for t, r in unusable.iterrows()]
            f.write(f'- **{lab}**: {", ".join(items)}\n')
        f.write('\nMentioned-in-BOTH-cohorts-unusable (structural paper-mentioned-vs-data gap):\n\n')
        hcc_unusable = set(hcc_cls[(hcc_cls['status'].isin(['fully_missing','invariant','insufficient'])) & (hcc_cls['mentioned_in_paper'])].index)
        crl_unusable = set(crl_cls[(crl_cls['status'].isin(['fully_missing','invariant','insufficient'])) & (crl_cls['mentioned_in_paper'])].index)
        both_unusable = hcc_unusable & crl_unusable
        for t in sorted(both_unusable):
            hst = hcc_cls.loc[t, 'status']; cst = crl_cls.loc[t, 'status']
            f.write(f'- `{t}` — HCC {hst}, CRLM {cst}\n')
        f.write('\nEssential-26 as a whole nearly doubles observable variability:\n\n')
        n_hcc_ess_test = (hcc_cls['status']=='testable').sum()
        n_crl_ess_test = (crl_cls['status']=='testable').sum()
        n_hcc_ess_sig  = int(hcc_cls['significant_fdr'].sum())
        n_crl_ess_sig  = int(crl_cls['significant_fdr'].sum())
        f.write(f'| Cohort | Essential-26 testable | Essential-26 FDR-sig | Mentioned-11 testable | Mentioned-11 sig |\n')
        f.write(f'|---|---|---|---|---|\n')
        f.write(f'| HCC-TACE-Seg | {n_hcc_ess_test} | {n_hcc_ess_sig} | {n_hcc_test} | {len(hcc_m)} |\n')
        f.write(f'| CRLM         | {n_crl_ess_test} | {n_crl_ess_sig} | {n_crl_test} | {len(crl_m)} |\n\n')

        # Overview table
        f.write('## 4. Full classification — all 26 essential tags (HCC ‖ CRLM)\n\n')
        f.write(build_overview_table(hcc_cls, crl_cls))
        f.write('\n\n')

        # Per-tag bin tables
        f.write('## 5. Per-tag bin tables — HCC sens / mean-p  vs  CRLM spec / mean-p\n\n')
        f.write('Tables grouped by class, in order of essential spec.\n\n')
        current_class = None
        for tag, hexa, cls, ment in ESSENTIAL_ORDER:
            if cls != current_class:
                f.write(f'\n### ▸ {cls}\n')
                current_class = cls
            f.write(per_tag_bin_section(tag, hexa, cls, ment, hcc_bin, crl_bin, hcc_cls, crl_cls))
            f.write('\n')

    print(f'→ {OUT_MD}')


if __name__ == '__main__':
    main()
