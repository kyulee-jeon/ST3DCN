"""Extract essential-26 + helper CT DICOM tags from each patient's CT series.
Merges with inference results (lesion_results.csv) to create an augmented
`_with_tags.csv` for downstream per-tag heterogeneity analysis.

Input:
  - lesion_results.csv with columns including: patient_id, ct_uid
  - DICOM root path

Output:
  - lesion_results_with_tags.csv — per-lesion row, with essential-26 tag
    values (patient-level majority vote across DCMs, optionally filtered
    by AcquisitionNumber).
"""
import csv
import glob
from pathlib import Path

import pandas as pd
import pydicom


TAGS = [
    # essential-26
    ('Manufacturer',                  (0x0008, 0x0070)),
    ('ManufacturerModelName',         (0x0008, 0x1090)),
    ('BodyPartExamined',              (0x0018, 0x0015)),
    ('SliceThickness',                (0x0018, 0x0050)),
    ('KVP',                           (0x0018, 0x0060)),
    ('SpacingBetweenSlices',          (0x0018, 0x0088)),
    ('ContrastBolusAgent',            (0x0018, 0x0010)),
    ('ContrastFlowRate',              (0x0018, 0x1046)),
    ('ReconstructionDiameter',        (0x0018, 0x1100)),
    ('XRayTubeCurrent',               (0x0018, 0x1151)),
    ('Exposure_uAs',                  (0x0018, 0x1153)),
    ('ConvolutionKernel',             (0x0018, 0x1210)),
    ('PatientPosition',               (0x0018, 0x5100)),
    ('RevolutionTime',                (0x0018, 0x9305)),
    ('TotalCollimationWidth',         (0x0018, 0x9307)),
    ('TableSpeed',                    (0x0018, 0x9309)),
    ('SpiralPitchFactor',             (0x0018, 0x9311)),
    ('ContrastBolusAgentPhase',       (0x0018, 0x9344)),
    ('CTDIvol',                       (0x0018, 0x9345)),
    ('ReconstructionAlgorithm',       (0x0018, 0x9315)),
    ('Rows',                          (0x0028, 0x0010)),
    ('Columns',                       (0x0028, 0x0011)),
    ('PixelSpacing',                  (0x0028, 0x0030)),
    ('WindowCenter',                  (0x0028, 0x1050)),
    ('WindowWidth',                   (0x0028, 0x1051)),
    ('ReformattingThickness',         (0x0072, 0x0512)),
    # extras used for derivation / context
    ('ExposureTime',                  (0x0018, 0x1150)),
    ('Exposure',                      (0x0018, 0x1152)),
    ('ContrastBolusVolume',           (0x0018, 0x1041)),
    ('ContrastBolusRoute',            (0x0018, 0x1040)),
    ('AcquisitionNumber',             (0x0020, 0x0012)),
    ('SeriesDescription',             (0x0008, 0x103E)),
    ('StudyDate',                     (0x0008, 0x0020)),
    ('ProtocolName',                  (0x0018, 0x1030)),
]


def _extract_series_tags(dicom_root: str, series_uid: str, acq_num=None) -> dict:
    """Majority vote per tag across DCMs in the series. If acq_num given,
    restrict to files matching that AcquisitionNumber."""
    files = glob.glob(f'{dicom_root}/{series_uid}/*.dcm')
    counts = {name: {} for name, _ in TAGS}
    for f in files:
        try:
            d = pydicom.dcmread(f, stop_before_pixels=True)
        except Exception:
            continue
        if acq_num is not None:
            a = int(getattr(d, 'AcquisitionNumber', 1) or 1)
            if a != acq_num:
                continue
        for name, tag in TAGS:
            if tag in d:
                v = d[tag].value
                if hasattr(v, 'original_string'):
                    v = str(v)
                if isinstance(v, (bytes, bytearray)):
                    v = v.decode(errors='replace')
                if not isinstance(v, (str, int, float)):
                    v = str(v)
                counts[name][v] = counts[name].get(v, 0) + 1
    return {name: (max(counts[name], key=counts[name].get) if counts[name] else '')
            for name in counts}


def tag_dataframe(in_csv: str, out_csv: str, dicom_root: str,
                    uid_col: str = 'ct_uid', acq_col: str = None):
    """Tag every row's CT series. Caches per-series results."""
    df = pd.read_csv(in_csv)
    all_tag_names = [name for name, _ in TAGS]
    cache = {}
    out_rows = []
    for idx, row in df.iterrows():
        out = row.to_dict()
        uid = str(row.get(uid_col, '')).strip()
        if row.get('status') == 'ok' and uid:
            acq = None
            if acq_col and acq_col in df.columns:
                try: acq = int(row[acq_col])
                except (ValueError, TypeError): acq = None
            key = (uid, acq)
            if key not in cache:
                cache[key] = _extract_series_tags(dicom_root, uid, acq_num=acq)
            out.update(cache[key])
        else:
            for name in all_tag_names:
                out[name] = ''
        out_rows.append(out)
    fields = list(df.columns) + all_tag_names
    with open(out_csv, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(out_rows)
    print(f'tag-augmented → {out_csv}')
