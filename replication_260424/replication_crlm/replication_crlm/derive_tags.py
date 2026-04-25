"""Derive essential-tag values that are fully-missing but computable from
available DICOM fields. DICOM-standards-based derivations only.

Derivations (applied to both HCC-TACE-Seg and CRLM cohorts):

  1. SpacingBetweenSlices (0018,0088)
     ← median(diff(sorted ImagePositionPatient[2] across slices in CT series))
     Reference: DICOM PS3.3 §C.7.6.2.1.1

  2. Exposure_uAs (0018,1153)
     ← ExposureTime (0018,1150, ms) × XRayTubeCurrent (0018,1151, mA)
     Reference: DICOM PS3.3 §C.8.7.2.1 — standard itself gives this as the
     example derivation. Empirically chosen over `Exposure × 1000` because
     the Exposure field (0018,1152) exhibits vendor-inconsistent semantics
     (per-rotation vs total-scan), while Time × Current is physically
     well-defined and shows lower intra-scanner CV (0.2-0.3 vs 0.8-1.3).

  3. [HCC only] Manufacturer, ManufacturerModelName
     ← NBIA digest Excel columns (direct DICOM field fill <2% in HCC-TACE-Seg)

Output: writes `<cohort>_lesion_results_with_derived.csv` alongside input,
adding columns `SpacingBetweenSlices_derived`, `Exposure_uAs_derived`, and
(HCC only) `Manufacturer_derived`, `ManufacturerModelName_derived`.

Note: ContrastBolusAgentPhase is implicitly captured via the pipeline's
SEG-ref + Portal Vein HU check (all series are PORTAL_VENOUS in both cohorts),
so the derived value would be a constant → invariant, not written as a column.
"""
import glob
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pydicom


# Cohort configuration
COHORTS = {
    'hcc': {
        'dicom_root': '/home/ubuntu/hcc_data/dicom',
        'lesion_csv': '/home/ubuntu/hcc_workspace/replication/lesion_results_with_tags.csv',
        'patient_csv': '/home/ubuntu/hcc_workspace/replication/batch_results_with_tags.csv',
        'out_lesion': '/home/ubuntu/hcc_workspace/replication/lesion_results_with_derived.csv',
        'out_patient': '/home/ubuntu/hcc_workspace/replication/batch_results_with_derived.csv',
        'nbia_digest': '/home/ubuntu/hcc_data/HCC-TACE-Seg_v1_202201-nbia-digest.xlsx',
    },
    'crlm': {
        'dicom_root': '/home/ubuntu/non-hcc_data/Colorectal-Liver-Metastases',
        'lesion_csv': '/home/ubuntu/hcc_workspace/replication_crlm/lesion_results_with_tags.csv',
        'patient_csv': '/home/ubuntu/hcc_workspace/replication_crlm/batch_results_with_tags.csv',
        'out_lesion': '/home/ubuntu/hcc_workspace/replication_crlm/lesion_results_with_derived.csv',
        'out_patient': '/home/ubuntu/hcc_workspace/replication_crlm/batch_results_with_derived.csv',
        'nbia_digest': None,  # CRLM doesn't need NBIA join (direct fill is 100%)
    },
}


def derive_spacing_between_slices(ct_uid: str, dicom_root: str):
    """Return median inter-slice z-spacing (mm) from ImagePositionPatient."""
    zs = []
    for f in glob.glob(f'{dicom_root}/{ct_uid}/*.dcm'):
        try:
            d = pydicom.dcmread(f, stop_before_pixels=True,
                                specific_tags=['ImagePositionPatient'])
        except Exception:
            continue
        ipp = getattr(d, 'ImagePositionPatient', None)
        if ipp is None or len(ipp) < 3:
            continue
        try:
            zs.append(float(ipp[2]))
        except (TypeError, ValueError):
            pass
    if len(zs) < 2:
        return None
    zs = sorted(zs)
    diffs = np.diff(zs)
    diffs = diffs[diffs > 0]  # drop any duplicate positions
    if len(diffs) == 0:
        return None
    return float(np.median(np.abs(diffs)))


def derive_cohort(cohort_key: str):
    cfg = COHORTS[cohort_key]
    print(f'\n=== {cohort_key.upper()} ===')
    lesion = pd.read_csv(cfg['lesion_csv'])
    patient = pd.read_csv(cfg['patient_csv'])

    # --- (1) SpacingBetweenSlices per ct_uid (compute once, broadcast) ---
    unique_cts = patient[patient['status'] == 'ok']['ct_uid'].dropna().unique()
    print(f'[1/3] SpacingBetweenSlices ← IPP  — scanning {len(unique_cts)} CT series...')
    sbs_by_ct = {}
    for i, ct_uid in enumerate(unique_cts, 1):
        sbs = derive_spacing_between_slices(ct_uid, cfg['dicom_root'])
        sbs_by_ct[ct_uid] = sbs
        if i % 25 == 0:
            print(f'     {i}/{len(unique_cts)}')
    fill_sbs = sum(1 for v in sbs_by_ct.values() if v is not None) / max(1, len(sbs_by_ct))
    print(f'  derived fill: {fill_sbs*100:.1f}%')

    # --- (2) Exposure_uAs = ExposureTime × XRayTubeCurrent (row-level) ---
    print(f'[2/3] Exposure_uAs ← ExposureTime × XRayTubeCurrent')
    for df, label in [(lesion, 'lesion'), (patient, 'patient')]:
        et = pd.to_numeric(df.get('ExposureTime'), errors='coerce')
        cu = pd.to_numeric(df.get('XRayTubeCurrent'), errors='coerce')
        df['Exposure_uAs_derived'] = et * cu  # ms × mA = μAs
        fill = df['Exposure_uAs_derived'].notna().sum()
        ok_n = (df['status'] == 'ok').sum() if 'status' in df else len(df)
        print(f'  {label:8s}: {fill}/{ok_n}  ({fill/max(1,ok_n)*100:.1f}%)')

    # broadcast SpacingBetweenSlices_derived
    for df in (lesion, patient):
        df['SpacingBetweenSlices_derived'] = df['ct_uid'].map(sbs_by_ct)

    # --- (3) Manufacturer / ModelName (HCC only, via NBIA digest) ---
    if cfg.get('nbia_digest'):
        print(f'[3/3] Manufacturer / ModelName ← NBIA digest')
        digest = pd.read_excel(cfg['nbia_digest'])
        # digest cols usually: 'Patient ID', 'Study Instance UID', 'Series Instance UID',
        # 'Modality', 'Manufacturer', 'Manufacturer Model Name', ...
        key_cols = ['Series Instance UID', 'Manufacturer', 'Manufacturer Model Name']
        avail = [c for c in key_cols if c in digest.columns]
        print(f'  digest cols available: {avail}')
        if set(key_cols).issubset(digest.columns):
            mfr = digest.set_index('Series Instance UID')[['Manufacturer',
                                                             'Manufacturer Model Name']].to_dict('index')
            def lookup(ct_uid, field):
                row = mfr.get(ct_uid)
                if row is None: return None
                v = row.get(field)
                if pd.isna(v): return None
                return str(v).strip() or None
            for df in (lesion, patient):
                df['Manufacturer_derived']          = df['ct_uid'].apply(
                    lambda u: lookup(u, 'Manufacturer'))
                df['ManufacturerModelName_derived'] = df['ct_uid'].apply(
                    lambda u: lookup(u, 'Manufacturer Model Name'))
            fill_m = patient['Manufacturer_derived'].notna().sum()
            fill_mn = patient['ManufacturerModelName_derived'].notna().sum()
            ok_n = (patient['status']=='ok').sum()
            print(f'  patient: Manufacturer {fill_m}/{ok_n} ({fill_m/max(1,ok_n)*100:.1f}%)  '
                  f'ModelName {fill_mn}/{ok_n} ({fill_mn/max(1,ok_n)*100:.1f}%)')
        else:
            print('  digest columns missing — skipping')
    else:
        print(f'[3/3] NBIA digest derivation skipped for {cohort_key} (direct fill already ≥99%)')

    lesion.to_csv(cfg['out_lesion'], index=False)
    patient.to_csv(cfg['out_patient'], index=False)
    print(f'→ {cfg["out_lesion"]}')
    print(f'→ {cfg["out_patient"]}')


def main():
    derive_cohort('hcc')
    derive_cohort('crlm')


if __name__ == '__main__':
    main()
