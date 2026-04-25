"""CRLM lesion-level subset SPECIFICITY analysis (pu / pc / pce).

Same matching logic as replication/subset_analysis.py. Adds space-insensitive
ManufacturerModelName normalization (fixes "LightSpeed VCT" vs "Light Speed VCT"
style mismatches seen in CRLM — this is a genuine normalization fix, not a
relaxation of the spec).

Runs both STRICT (original matcher + space fix only) and RELAXED (additional
tube-current tolerance 5%) to show how pc is distribution-sensitive.

Unit = lesion; patient-level pc/pce flag broadcasts to all of that patient's lesions.
"""
import sys
import importlib.util
import numpy as np
import pandas as pd

spec = importlib.util.spec_from_file_location(
    'subset_analysis', '/home/ubuntu/hcc_workspace/replication/subset_analysis.py')
sa = importlib.util.module_from_spec(spec); spec.loader.exec_module(sa)

IN_LESION = '/home/ubuntu/hcc_workspace/replication_crlm/lesion_results_with_tags.csv'
IN_PATIENT = '/home/ubuntu/hcc_workspace/replication_crlm/batch_results_with_tags.csv'
OUT_SUMMARY = '/home/ubuntu/hcc_workspace/replication_crlm/subset_summary_lesion.txt'
OUT_DETAIL = '/home/ubuntu/hcc_workspace/replication_crlm/subset_detail_lesion.csv'
THRESHOLD = 0.8


def wilson(p, n, z=1.96):
    if n == 0: return (0.0, 0.0)
    d = 1 + z**2 / n
    c = (p + z**2 / (2*n)) / d
    h = z * ((p*(1-p)/n + z**2/(4*n**2))**0.5) / d
    return max(0.0, c-h), min(1.0, c+h)


_orig_value_matches = sa.value_matches


def patched_value_matches(tag_hex, patient_val, protocol_val,
                            tube_current_tol_pct=None, agent_stem_match=False):
    """Patched matcher for CRLM normalization quirks.

    - ModelName: space-insensitive substring (fixes "LightSpeed VCT" vs "Light Speed VCT")
    - [optional] Tube current: relative tolerance (default strict ±0.5)
    - [optional] Agent: brand-stem match — MSKCC uses "OMNI" for "Omnipaque",
      "IOPAM" for "Iopamidol", etc. Consider a 4+ char prefix overlap a match.
    """
    pv = sa.norm_value(patient_val)
    ov = sa.norm_value(protocol_val)
    if pv is None or ov is None:
        return None

    if tag_hex == '00081090':  # ManufacturerModelName
        pv_n = ''.join(pv.lower().split())
        ov_n = ''.join(ov.lower().split())
        if '(' in ov_n and 'only' in ov_n: return None
        return pv_n in ov_n or ov_n in pv_n

    if tag_hex == '00181151' and tube_current_tol_pct is not None:
        pfv, ofv = sa.try_float(pv), sa.try_float(ov)
        if pfv is None or ofv is None: return None
        if isinstance(ov, str) and ov.lower() == 'automatic': return None
        return abs(pfv - ofv) / max(ofv, 1.0) <= tube_current_tol_pct

    if tag_hex == '00180010' and agent_stem_match:
        import re as _re
        def stems(s):
            out = set()
            for t in _re.split(r'[/\s&,]+', s.lower()):
                t = _re.sub(r'\d+', '', t).strip()
                if len(t) >= 4:
                    out.add(t[:4])  # 4-char brand stem (omni = omnipaque[:4])
            return out
        return bool(stems(pv) & stems(ov))

    return _orig_value_matches(tag_hex, patient_val, protocol_val)


def check_pc_patched(row, protocols, essential_only=False, essential_set=None,
                      tube_current_tol_pct=None, agent_stem_match=False):
    """Copy of sa.check_pc using patched_value_matches (with optional relaxations)."""
    best_match = False; best_p = None; per_tag_final = {}
    for p_idx in range(15):
        all_present_match = True; per_tag = {}
        for tag_hex, name, values in protocols:
            full_hex_8 = tag_hex.zfill(8).lower()
            if full_hex_8 in sa.EXCLUDE_TAGS or tag_hex in sa.EXCLUDE_TAGS:
                continue
            if essential_only and essential_set is not None:
                forms = {full_hex_8, full_hex_8.lstrip('0').zfill(8),
                         tag_hex.lstrip('0').zfill(8)}
                if not (essential_set & forms):
                    continue
            col_name = sa.PATIENT_TAG_FROM_CSV.get(full_hex_8) or sa.PATIENT_TAG_FROM_CSV.get(tag_hex)
            if col_name is None:
                continue
            patient_val = row.get(col_name)
            protocol_val = values[p_idx]
            m = patched_value_matches(full_hex_8, patient_val, protocol_val,
                                       tube_current_tol_pct=tube_current_tol_pct,
                                       agent_stem_match=agent_stem_match)
            if m is None:
                per_tag[name] = 'skip'
            elif m:
                per_tag[name] = 'match'
            else:
                per_tag[name] = 'mismatch'; all_present_match = False
        checked = sum(1 for v in per_tag.values() if v != 'skip')
        if all_present_match and checked > 0:
            best_match = True; best_p = p_idx + 1; per_tag_final = per_tag; break
    if not best_match:
        per_tag_final = per_tag
    return best_match, best_p, per_tag_final


def summarize(name, probs, f_out):
    n = len(probs)
    if n == 0:
        f_out.write(f'[{name:38s}] n=   0\n')
        return {'subset': name, 'n': 0}
    fp = int((probs >= THRESHOLD).sum())
    tn = n - fp
    spec = tn / n
    lo, hi = wilson(spec, n)
    f_out.write(
        f'[{name:38s}] n={n:4d}  TN={tn:4d}  FP={fp:4d}  '
        f'spec={spec:.4f}  mean={probs.mean():.4f}  med={np.median(probs):.4f}  '
        f'Wilson95%=({lo:.4f}, {hi:.4f})\n'
    )
    return {'subset': name, 'n': n, 'tn': tn, 'fp': fp, 'spec': spec,
            'mean_prob': float(probs.mean()), 'median_prob': float(np.median(probs)),
            'ci_lo': lo, 'ci_hi': hi}


def main():
    protocols = sa.load_protocols()
    ess = sa.load_essential_tags()
    ess_forms = set()
    for e in ess:
        ess_forms.add(e.zfill(8)); ess_forms.add(e.lstrip('0').zfill(8))

    pat = pd.read_csv(IN_PATIENT)
    pat = pat[pat['status'] == 'ok'].copy()

    les = pd.read_csv(IN_LESION)
    les = les[les['status'] == 'ok'].copy()
    les['prob'] = les['prob'].astype(float)

    rows = []
    with open(OUT_SUMMARY, 'w') as f:
        f.write('=== CRLM lesion-level subset SPECIFICITY (ST3DCN) ===\n')
        f.write(f'Threshold={THRESHOLD}    (all CRLM patients HCC-negative; spec = TN/n)\n')
        f.write(f'Excluded: SliceThickness (protocol ≤1.25mm vs CRLM 2.5-7.5mm), '
                f'AgentSequence / ContrastPhase / AcquisitionType / InjectionDelay.\n\n')
        f.write(f'Total lesions: {len(les)}  from {les["patient_id"].nunique()} patients\n')

        variants = [
            ('STRICT (space-normalized model name only)', None, False),
            ('MEDIUM (+10% tube-current tolerance)',      0.10, False),
            ('RELAXED (+MSKCC agent brand-stem match)',   0.10, True),
        ]
        for variant_name, tc_tol, agent_stem in variants:
            f.write(f'\n\n═════════ {variant_name} ═════════\n')
            pat_flags = {}
            for _, row in pat.iterrows():
                m_pc, p_pc, _ = check_pc_patched(row, protocols, essential_only=False,
                                                  tube_current_tol_pct=tc_tol,
                                                  agent_stem_match=agent_stem)
                m_pce, p_pce, _ = check_pc_patched(row, protocols, essential_only=True,
                                                    essential_set=ess_forms,
                                                    tube_current_tol_pct=tc_tol,
                                                    agent_stem_match=agent_stem)
                pat_flags[row['patient_id']] = (m_pc, p_pc, m_pce, p_pce)

            les_v = les.copy()
            les_v['in_pc']  = les_v['patient_id'].map(lambda p: pat_flags.get(p, (False,None,False,None))[0])
            les_v['pc_proto']  = les_v['patient_id'].map(lambda p: pat_flags.get(p, (False,None,False,None))[1])
            les_v['in_pce'] = les_v['patient_id'].map(lambda p: pat_flags.get(p, (False,None,False,None))[2])
            les_v['pce_proto'] = les_v['patient_id'].map(lambda p: pat_flags.get(p, (False,None,False,None))[3])

            pu, pc, pce = les_v, les_v[les_v['in_pc']], les_v[les_v['in_pce']]
            n_pc_pat  = sum(1 for v in pat_flags.values() if v[0])
            n_pce_pat = sum(1 for v in pat_flags.values() if v[2])

            f.write('\n--- Subset specificities (lesion = unit) ---\n')
            r1 = summarize('pu',                     pu['prob'].values, f);  r1['variant'] = variant_name
            r2 = summarize('pc  (ST3DCN-14, -thick)', pc['prob'].values, f); r2['variant'] = variant_name
            r3 = summarize('pce (ST3DCN∩essential, -thick)', pce['prob'].values, f); r3['variant'] = variant_name
            rows.extend([r1, r2, r3])
            f.write(f'\n--- patient counts: pc={n_pc_pat}/{len(pat_flags)},  '
                    f'pce={n_pce_pat}/{len(pat_flags)} ---\n')

            if len(pc):
                proto_counts = pc.groupby('patient_id')['pc_proto'].first().value_counts().sort_index()
                f.write('--- pc protocol assignments ---\n')
                for p, c in proto_counts.items():
                    n_les = int((pc['pc_proto'] == p).sum())
                    f.write(f'  P{int(p):02d}: {c} patients, {n_les} lesions\n')
            if len(pce):
                proto_counts = pce.groupby('patient_id')['pce_proto'].first().value_counts().sort_index()
                f.write('--- pce protocol assignments ---\n')
                for p, c in proto_counts.items():
                    n_les = int((pce['pce_proto'] == p).sum())
                    f.write(f'  P{int(p):02d}: {c} patients, {n_les} lesions\n')

    pd.DataFrame(rows).to_csv(OUT_DETAIL, index=False)
    print(f'summary → {OUT_SUMMARY}')
    print(f'detail  → {OUT_DETAIL}\n')
    with open(OUT_SUMMARY) as f:
        sys.stdout.write(f.read())


if __name__ == '__main__':
    main()
