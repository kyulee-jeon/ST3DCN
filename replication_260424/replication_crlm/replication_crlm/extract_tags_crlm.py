"""Extract essential-26 CT tags from CRLM CT series and merge with batch_results.csv.

Same tag set as /home/ubuntu/hcc_workspace/replication/extract_tags.py.
Since CRLM has 1 CT/patient with 1 AcqNum, no acq filtering needed.
"""
import csv
import glob

import pandas as pd
import pydicom


DICOM_ROOT = '/home/ubuntu/non-hcc_data/Colorectal-Liver-Metastases'
BATCH_CSV = '/home/ubuntu/hcc_workspace/replication_crlm/batch_results.csv'
OUT_CSV = '/home/ubuntu/hcc_workspace/replication_crlm/batch_results_with_tags.csv'
LESION_CSV = '/home/ubuntu/hcc_workspace/replication_crlm/lesion_results.csv'
LESION_OUT_CSV = '/home/ubuntu/hcc_workspace/replication_crlm/lesion_results_with_tags.csv'

TAGS = [
    ('Manufacturer', (0x0008, 0x0070)),
    ('ManufacturerModelName', (0x0008, 0x1090)),
    ('BodyPartExamined', (0x0018, 0x0015)),
    ('SliceThickness', (0x0018, 0x0050)),
    ('KVP', (0x0018, 0x0060)),
    ('SpacingBetweenSlices', (0x0018, 0x0088)),
    ('ContrastBolusAgent', (0x0018, 0x0010)),
    ('ContrastFlowRate', (0x0018, 0x1046)),
    ('ReconstructionDiameter', (0x0018, 0x1100)),
    ('XRayTubeCurrent', (0x0018, 0x1151)),
    ('Exposure_uAs', (0x0018, 0x1153)),
    ('ConvolutionKernel', (0x0018, 0x1210)),
    ('PatientPosition', (0x0018, 0x5100)),
    ('RevolutionTime', (0x0018, 0x9305)),
    ('TotalCollimationWidth', (0x0018, 0x9307)),
    ('TableSpeed', (0x0018, 0x9309)),
    ('SpiralPitchFactor', (0x0018, 0x9311)),
    ('ContrastBolusAgentPhase', (0x0018, 0x9344)),
    ('CTDIvol', (0x0018, 0x9345)),
    ('ReconstructionAlgorithm', (0x0018, 0x9315)),
    ('Rows', (0x0028, 0x0010)),
    ('Columns', (0x0028, 0x0011)),
    ('PixelSpacing', (0x0028, 0x0030)),
    ('WindowCenter', (0x0028, 0x1050)),
    ('WindowWidth', (0x0028, 0x1051)),
    ('ReformattingThickness', (0x0072, 0x0512)),
    ('ExposureTime', (0x0018, 0x1150)),
    ('Exposure', (0x0018, 0x1152)),
    ('ContrastBolusVolume', (0x0018, 0x1041)),
    ('AcquisitionNumber', (0x0020, 0x0012)),
    ('SeriesDescription', (0x0008, 0x103E)),
    ('StudyDate', (0x0008, 0x0020)),
    ('ContrastBolusRoute', (0x0018, 0x1040)),
    ('ProtocolName', (0x0018, 0x1030)),
]


def extract_ct_tags(ct_uid: str):
    files = glob.glob(f'{DICOM_ROOT}/{ct_uid}/*.dcm')
    counts = {name: {} for name, _ in TAGS}
    for f in files:
        try:
            d = pydicom.dcmread(f, stop_before_pixels=True)
        except Exception:
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


def tag_df(in_csv, out_csv, uid_col='ct_uid'):
    df = pd.read_csv(in_csv)
    print(f'{in_csv}: {len(df)} rows')
    all_tag_names = [name for name, _ in TAGS]
    cache = {}
    out_rows = []
    for idx, row in df.iterrows():
        out = row.to_dict()
        uid = str(row.get(uid_col, '')).strip()
        if row.get('status') == 'ok' and uid:
            if uid not in cache:
                cache[uid] = extract_ct_tags(uid)
            out.update(cache[uid])
        else:
            for name in all_tag_names:
                out[name] = ''
        out_rows.append(out)
        if (idx + 1) % 20 == 0:
            print(f'  tagged {idx+1}/{len(df)}')

    fields = list(df.columns) + all_tag_names
    with open(out_csv, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(out_rows)
    print(f'written → {out_csv}')


def main():
    import os
    tag_df(BATCH_CSV, OUT_CSV)
    if os.path.exists(LESION_CSV):
        tag_df(LESION_CSV, LESION_OUT_CSV)


if __name__ == '__main__':
    main()
