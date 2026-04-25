"""ST3DCN replication pipeline for Colorectal-Liver-Metastases (CRLM).

Differences from HCC-TACE-Seg pipeline (replication/st3dcn_pipeline.py):
  - Each patient: exactly 1 CT + 1 SEG, 1 study (pre-op per TCIA) → no baseline/offset logic
  - Each CT series: 1 AcquisitionNumber → no arterial/portal AcqNum selection
  - SEG segment layout:
      1 = Liver, 2 = Liver Remnant, 3 = Hepatic Vein, 4 = Portal Vein,
      5..9 = Tumor_1..Tumor_5 (PropType = Mass)
    → Tumor mask = union of segments 5..9 (per-patient), or individual segment (per-lesion)
  - Phase verification: empirical — portal vein (seg 4) median HU should be > 120
    (expected in portal venous phase; sanity flag, not exclusion)
Everything else (SOP-UID → CT alignment, 70³ crop, HU window [-160, 240], resize) is identical.
"""
import glob
import os
from collections import defaultdict
from dataclasses import dataclass

import numpy as np
import pandas as pd
import pydicom
from scipy.ndimage import zoom


DICOM_ROOT = '/home/ubuntu/non-hcc_data/Colorectal-Liver-Metastases'
MANIFEST_JSON = f'{DICOM_ROOT}/series_manifest.json'

TUMOR_SEGMENTS = {5, 6, 7, 8, 9}  # Tumor_1..Tumor_5
PORTAL_VEIN_SEGMENT = 4

WINDOW_LEVEL = 40
WINDOW_WIDTH = 400
HU_MIN = WINDOW_LEVEL - WINDOW_WIDTH / 2  # -160
HU_MAX = WINDOW_LEVEL + WINDOW_WIDTH / 2  # 240
CROP_SIZE = 70


def load_manifest():
    """Return DataFrame with columns analogous to HCC manifest."""
    import json
    with open(MANIFEST_JSON) as f:
        rows = json.load(f)
    df = pd.DataFrame(rows)
    return df


def patient_series(manifest: pd.DataFrame, patient_id: str):
    """Return (ct_uid, seg_uid, study_uid) for the one-per-patient CRLM case."""
    p = manifest[manifest['PatientID'] == patient_id]
    assert len(p) == 2, f'{patient_id}: expected 1 CT + 1 SEG, got {len(p)}'
    ct = p[p['Modality'] == 'CT']
    seg = p[p['Modality'] == 'SEG']
    assert len(ct) == 1 and len(seg) == 1, f'{patient_id}: CT={len(ct)} SEG={len(seg)}'
    return (ct.iloc[0]['SeriesInstanceUID'], seg.iloc[0]['SeriesInstanceUID'],
            ct.iloc[0]['StudyInstanceUID'])


def read_seg(seg_uid: str):
    """Return (ds, pixel_array, frame_meta, segment_labels).

    frame_meta: [{frame, segment, ref_sop (or None), z (float), mask (2D bool)}]
    segment_labels: {seg_num -> label string}
    """
    files = glob.glob(f'{DICOM_ROOT}/{seg_uid}/*.dcm')
    assert len(files) == 1, f'Expected 1 SEG file, got {len(files)}'
    ds = pydicom.dcmread(files[0])
    px = ds.pixel_array
    frame_meta = []
    for i, fg in enumerate(ds.PerFrameFunctionalGroupsSequence):
        seg_num = int(fg.SegmentIdentificationSequence[0].ReferencedSegmentNumber)
        ref_sop = None
        try:
            ref_sop = (fg.DerivationImageSequence[0]
                         .SourceImageSequence[0]
                         .ReferencedSOPInstanceUID)
        except (AttributeError, IndexError, KeyError):
            pass
        z = float(fg.PlanePositionSequence[0].ImagePositionPatient[2])
        frame_meta.append({'frame': i, 'segment': seg_num, 'ref_sop': ref_sop,
                           'z': z, 'mask': px[i].astype(bool)})
    seg_labels = {int(s.SegmentNumber): str(s.SegmentLabel) for s in ds.SegmentSequence}
    return ds, px, frame_meta, seg_labels


def index_ct_series(ct_uid: str):
    """Return dict: sop_uid -> {z, hu, pix_spacing, slice_thick, for_uid}"""
    out = {}
    for f in glob.glob(f'{DICOM_ROOT}/{ct_uid}/*.dcm'):
        d = pydicom.dcmread(f)
        slope = float(getattr(d, 'RescaleSlope', 1.0))
        intercept = float(getattr(d, 'RescaleIntercept', 0.0))
        out[d.SOPInstanceUID] = {
            'z': float(d.ImagePositionPatient[2]),
            'hu': d.pixel_array.astype(np.float32) * slope + intercept,
            'pix_spacing': tuple(float(x) for x in d.PixelSpacing),
            'slice_thick': float(d.SliceThickness),
            'for_uid': getattr(d, 'FrameOfReferenceUID', None),
        }
    return out


def validate_sop_linkage(frame_meta, ct_index):
    """SEG must reference CT SOPInstanceUIDs. Return (n_matched, n_total_frames_with_ref)."""
    refs = [fm['ref_sop'] for fm in frame_meta if fm['ref_sop'] is not None]
    matched = sum(1 for r in refs if r in ct_index)
    return matched, len(refs)


def build_ct_volume(ct_index):
    """Return (ct_vol, z_coords, pix_spacing) z-sorted ascending."""
    rows = sorted(ct_index.items(), key=lambda r: r[1]['z'])
    ct_vol = np.stack([m['hu'] for _, m in rows], axis=0)
    z_coords = np.array([m['z'] for _, m in rows])
    pix = rows[0][1]['pix_spacing']
    return ct_vol, z_coords, pix


def build_segment_mask(frame_meta, z_coords, segments, z_tol_mm=1.25):
    """Union of frames whose segment ∈ `segments`, aligned to z_coords nearest-neighbor.

    `segments` can be int or iterable of ints.
    Returns (Z, H, W) uint8 mask, or None if no frames matched.
    """
    if isinstance(segments, int):
        segments = {segments}
    else:
        segments = set(segments)
    z_arr = np.asarray(z_coords)
    H = W = None
    slices = {}
    for fm in frame_meta:
        if fm['segment'] not in segments:
            continue
        if not fm['mask'].any():
            continue
        diffs = np.abs(z_arr - fm['z'])
        idx = int(np.argmin(diffs))
        if diffs[idx] > z_tol_mm:
            continue
        if idx in slices:
            slices[idx] |= fm['mask']
        else:
            slices[idx] = fm['mask'].copy()
        H, W = fm['mask'].shape
    if not slices:
        return None
    mask = np.zeros((len(z_coords), H, W), dtype=np.uint8)
    for i, m in slices.items():
        mask[i] = m.astype(np.uint8)
    return mask


def portal_vein_median_hu(frame_meta, ct_index):
    """Sanity check: median HU inside portal vein segment.
    Expected > ~120 HU for portal venous phase (arterial phase: lower)."""
    vals = []
    for fm in frame_meta:
        if fm['segment'] != PORTAL_VEIN_SEGMENT or not fm['mask'].any():
            continue
        if fm['ref_sop'] and fm['ref_sop'] in ct_index:
            hu = ct_index[fm['ref_sop']]['hu']
        else:
            # fallback: z-coord
            candidates = [(abs(m['z'] - fm['z']), m['hu']) for m in ct_index.values()]
            candidates.sort()
            if candidates and candidates[0][0] < 1.25:
                hu = candidates[0][1]
            else:
                continue
        vals.extend(hu[fm['mask']].tolist())
    if not vals:
        return None
    return float(np.median(vals))


def crop_and_window(ct_vol, mass_mask, out_size=CROP_SIZE, pad_lo=10, pad_hi=9):
    """Same as HCC pipeline: tight bbox + asymmetric padding, HU window, 70³ resize."""
    nz = np.argwhere(mass_mask)
    assert nz.size, 'Empty tumor mask'
    z0, y0, x0 = nz.min(0)
    z1, y1, x1 = nz.max(0) + 1
    Z, Y, X = ct_vol.shape
    z0 = max(0, z0 - pad_lo); y0 = max(0, y0 - pad_lo); x0 = max(0, x0 - pad_lo)
    z1 = min(Z, z1 + pad_hi); y1 = min(Y, y1 + pad_hi); x1 = min(X, x1 + pad_hi)
    ct_crop = ct_vol[z0:z1, y0:y1, x0:x1].astype(np.float32)
    ct_crop = np.clip(ct_crop, HU_MIN, HU_MAX)
    ct_crop = ((ct_crop - HU_MIN) / (HU_MAX - HU_MIN) * 255.0).astype(np.uint8)
    ct_crop = np.transpose(ct_crop, (1, 2, 0))  # (Y, X, Z) paper convention
    arr = ct_crop.astype(np.float32)
    if max(arr.shape) > out_size:
        zf = tuple(out_size / s for s in arr.shape)
        arr = zoom(arr, zf, order=1)
    else:
        r = np.zeros((out_size, out_size, out_size), dtype=np.float32)
        o0 = (out_size - arr.shape[0]) // 2
        o1 = (out_size - arr.shape[1]) // 2
        o2 = (out_size - arr.shape[2]) // 2
        r[o0:o0+arr.shape[0], o1:o1+arr.shape[1], o2:o2+arr.shape[2]] = arr
        arr = r
    arr = np.clip(arr, 0, 255) / 255.0
    bbox = (int(z0), int(y0), int(x0), int(z1), int(y1), int(x1))
    return arr.astype(np.float32), bbox


def process_patient(patient_id: str, manifest=None):
    """Per-patient: tumor = UNION of segments 5..9 → single 70³ crop → one prob."""
    if manifest is None:
        manifest = load_manifest()
    ct_uid, seg_uid, study_uid = patient_series(manifest, patient_id)
    _ds, _px, frame_meta, seg_labels = read_seg(seg_uid)
    ct_index = index_ct_series(ct_uid)
    n_match, n_ref = validate_sop_linkage(frame_meta, ct_index)
    ct_vol, z_coords, pix = build_ct_volume(ct_index)
    tumor_mask = build_segment_mask(frame_meta, z_coords, TUMOR_SEGMENTS)
    assert tumor_mask is not None and tumor_mask.sum() > 0, \
        f'{patient_id}: empty tumor mask'
    pv_hu = portal_vein_median_hu(frame_meta, ct_index)
    crop, bbox = crop_and_window(ct_vol, tumor_mask)
    present_tumors = sorted(set(fm['segment'] for fm in frame_meta
                                 if fm['segment'] in TUMOR_SEGMENTS and fm['mask'].any()))
    return {
        'patient_id': patient_id,
        'study_uid': study_uid,
        'seg_uid': seg_uid,
        'ct_uid': ct_uid,
        'ct_shape': ct_vol.shape,
        'sop_matched': (n_match, n_ref),
        'portal_vein_median_hu': pv_hu,
        'tumor_segments_present': present_tumors,
        'tumor_nz_slices': int((tumor_mask.sum((1, 2)) > 0).sum()),
        'bbox': bbox,
        'pix_spacing_mm': pix,
        'slice_thick': list(ct_index.values())[0]['slice_thick'],
        'crop': crop,
    }


def process_patient_per_lesion(patient_id: str, manifest=None, min_voxels=50):
    """Per-lesion: EACH tumor segment 5..9 → separate 70³ crop → separate prob.

    Uses explicit segmentation labels (Tumor_1..Tumor_5) rather than connected
    components, since CRLM annotations already identify individual tumors.
    """
    if manifest is None:
        manifest = load_manifest()
    ct_uid, seg_uid, study_uid = patient_series(manifest, patient_id)
    _ds, _px, frame_meta, seg_labels = read_seg(seg_uid)
    ct_index = index_ct_series(ct_uid)
    n_match, n_ref = validate_sop_linkage(frame_meta, ct_index)
    ct_vol, z_coords, pix = build_ct_volume(ct_index)
    pv_hu = portal_vein_median_hu(frame_meta, ct_index)

    lesions = []
    for seg_num in sorted(TUMOR_SEGMENTS):
        m = build_segment_mask(frame_meta, z_coords, seg_num)
        if m is None or int(m.sum()) < min_voxels:
            continue
        crop, bbox = crop_and_window(ct_vol, m)
        lesions.append({
            'lesion_idx': len(lesions),
            'segment_number': seg_num,
            'segment_label': seg_labels.get(seg_num, f'Seg{seg_num}'),
            'voxel_count': int(m.sum()),
            'bbox': bbox,
            'crop': crop,
        })
    assert lesions, f'{patient_id}: no tumor lesions ≥ {min_voxels} voxels'
    return {
        'patient_id': patient_id,
        'study_uid': study_uid,
        'seg_uid': seg_uid,
        'ct_uid': ct_uid,
        'ct_shape': ct_vol.shape,
        'sop_matched': (n_match, n_ref),
        'portal_vein_median_hu': pv_hu,
        'pix_spacing_mm': pix,
        'slice_thick': list(ct_index.values())[0]['slice_thick'],
        'lesions': lesions,
    }
