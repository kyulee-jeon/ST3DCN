"""DICOM-standard derivations that rescue essential tags from fully_missing
status. Two exact derivations (both blessed by DICOM PS3.3):

  1. SpacingBetweenSlices (0018,0088)
       = median(diff(sorted ImagePositionPatient[2] across series))
     — PS3.3 §C.7.6.2.1.1

  2. Exposure_uAs (0018,1153)
       = ExposureTime (0018,1150, ms) × XRayTubeCurrent (0018,1151, mA)
     — PS3.3 §C.8.7.2.1 (standard-listed example calc).
     Chosen over `Exposure × 1000` because (0018,1152) Exposure has
     vendor-inconsistent semantics (per-rotation vs total-scan); empirically
     intra-scanner CV is 0.13-0.38 for Time×Current vs 0.74-1.31 for
     Exposure×1000.

  3. [HCC-TACE-Seg only] Manufacturer / ManufacturerModelName ← NBIA digest
     join (DICOM tags populated <2%, but external TCIA digest Excel may
     cover more — though empirically only 1/105 in HCC-TACE-Seg).

Outputs augmented CSV alongside input with `*_derived` columns.
"""
import glob

import numpy as np
import pandas as pd
import pydicom


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
    diffs = diffs[diffs > 0]
    if len(diffs) == 0:
        return None
    return float(np.median(np.abs(diffs)))


def derive_tags(lesion_csv: str, patient_csv: str, dicom_root: str,
                 out_lesion: str, out_patient: str,
                 nbia_digest: str = None,
                 nbia_series_col: str = 'Series Instance UID',
                 nbia_manu_col: str = 'Manufacturer',
                 nbia_model_col: str = 'Manufacturer Model Name'):
    """Add derived columns to both lesion-level and patient-level CSVs.

    Columns added:
      - SpacingBetweenSlices_derived   (float mm)
      - Exposure_uAs_derived            (float uAs)
      - Manufacturer_derived            (str, if nbia_digest given)
      - ManufacturerModelName_derived   (str, if nbia_digest given)
    """
    lesion = pd.read_csv(lesion_csv)
    patient = pd.read_csv(patient_csv)

    # SpacingBetweenSlices: compute once per unique ct_uid
    unique_cts = patient[patient['status'] == 'ok']['ct_uid'].dropna().unique()
    print(f'SpacingBetweenSlices: scanning {len(unique_cts)} CT series...')
    sbs = {uid: derive_spacing_between_slices(uid, dicom_root) for uid in unique_cts}
    sbs_fill = sum(1 for v in sbs.values() if v is not None) / max(1, len(sbs))
    print(f'  derived fill: {sbs_fill*100:.1f}%')

    for df in (lesion, patient):
        df['SpacingBetweenSlices_derived'] = df['ct_uid'].map(sbs)
        # Exposure_uAs = ExposureTime × XRayTubeCurrent
        et = pd.to_numeric(df.get('ExposureTime'), errors='coerce')
        cu = pd.to_numeric(df.get('XRayTubeCurrent'), errors='coerce')
        df['Exposure_uAs_derived'] = et * cu

    # NBIA digest join (optional, HCC-TACE-Seg specific)
    if nbia_digest:
        print(f'NBIA digest join: {nbia_digest}')
        digest = pd.read_excel(nbia_digest)
        if set([nbia_series_col, nbia_manu_col, nbia_model_col]).issubset(digest.columns):
            mfr = digest.set_index(nbia_series_col)[[nbia_manu_col, nbia_model_col]].to_dict('index')
            def lookup(uid, field):
                row = mfr.get(uid)
                if row is None: return None
                v = row.get(field)
                return None if pd.isna(v) else str(v).strip() or None
            for df in (lesion, patient):
                df['Manufacturer_derived']          = df['ct_uid'].apply(lambda u: lookup(u, nbia_manu_col))
                df['ManufacturerModelName_derived'] = df['ct_uid'].apply(lambda u: lookup(u, nbia_model_col))

    lesion.to_csv(out_lesion, index=False)
    patient.to_csv(out_patient, index=False)
    print(f'→ {out_lesion}')
    print(f'→ {out_patient}')
