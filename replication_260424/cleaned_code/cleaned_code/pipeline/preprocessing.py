"""ST3DCN preprocessing — paper-faithful, pixel-exact on paper sample data.

Pipeline:
  1. SEG DICOM → per-frame mask + ReferencedSOPInstanceUID
  2. Align SEG frames to CT slices (SOP-UID first; FrameOfReference + z-coord fallback)
  3. Tumor mask (union of tumor segments) → tight bbox
  4. Asymmetric padding: +10 low, +9 high per axis (reverse-engineered from paper's .npy)
  5. HU window [-160, 240] → uint8 [0, 255]
  6. Axis transpose (Z, Y, X) → (Y, X, Z) — paper's .npy convention
  7. Resize/pad to 70³, normalize to [0, 1] float32

Validated pixel-exact against ST3DCN repo's demo pre-cropped .npy files
(ID_0848 HCC, QEH032 non-HCC). |Δ| = 0 on all 10/10 observations.
"""
import glob
from dataclasses import dataclass
from typing import Iterable, Optional

import numpy as np
import pydicom
from scipy.ndimage import zoom


# Constants — match ST3DCN paper + repo
WINDOW_LEVEL = 40
WINDOW_WIDTH = 400
HU_MIN = WINDOW_LEVEL - WINDOW_WIDTH / 2      # -160
HU_MAX = WINDOW_LEVEL + WINDOW_WIDTH / 2      # 240
CROP_SIZE = 70
PAD_LOW = 10
PAD_HIGH = 9


@dataclass
class SeriesRef:
    """Pointer to a DICOM series on disk. dicom_root / series_uid / *.dcm"""
    dicom_root: str
    series_uid: str


def read_seg(seg: SeriesRef):
    """Parse a SEG DICOM. Returns:
        ds: pydicom Dataset
        pixel_array: (N_frames, H, W) uint8
        frame_meta: list of {'frame', 'segment', 'ref_sop' (or None), 'z', 'mask'}
        segment_labels: {segment_number -> label_string}
    """
    files = glob.glob(f'{seg.dicom_root}/{seg.series_uid}/*.dcm')
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
    seg_labels = {int(s.SegmentNumber): str(s.SegmentLabel)
                  for s in ds.SegmentSequence}
    return ds, px, frame_meta, seg_labels


def index_ct_series(ct: SeriesRef):
    """Read a CT series and index by SOPInstanceUID.
    Returns dict: sop_uid -> {acq, z, hu, pix_spacing, slice_thick, for_uid}"""
    out = {}
    for f in glob.glob(f'{ct.dicom_root}/{ct.series_uid}/*.dcm'):
        d = pydicom.dcmread(f)
        slope = float(getattr(d, 'RescaleSlope', 1.0))
        intercept = float(getattr(d, 'RescaleIntercept', 0.0))
        out[d.SOPInstanceUID] = {
            'acq': int(getattr(d, 'AcquisitionNumber', 1) or 1),
            'z': float(d.ImagePositionPatient[2]),
            'hu': d.pixel_array.astype(np.float32) * slope + intercept,
            'pix_spacing': tuple(float(x) for x in d.PixelSpacing),
            'slice_thick': float(d.SliceThickness),
            'for_uid': getattr(d, 'FrameOfReferenceUID', None),
        }
    return out


def validate_sop_linkage(frame_meta, ct_index):
    """Count how many SEG-frame ReferencedSOPInstanceUIDs resolve to this CT series.
    Returns (n_matched, n_total_frames_with_ref)."""
    refs = [fm['ref_sop'] for fm in frame_meta if fm['ref_sop'] is not None]
    matched = sum(1 for r in refs if r in ct_index)
    return matched, len(refs)


def pick_mask_aligned_ct(frame_meta, ct_candidates, seg_for_uid=None):
    """When multiple CT series exist for a patient (e.g., triphasic HCC-TACE-Seg),
    pick the one whose SOP-UIDs overlap SEG references (primary) or whose
    FrameOfReferenceUID + z-coordinates match (fallback).

    ct_candidates: list of (series_uid, series_description) or list of series_uids.
    Returns: (best_ct_uid, ct_index_dict, sop_overlap_set)
    """
    if ct_candidates and not isinstance(ct_candidates[0], (list, tuple)):
        ct_candidates = [(u, '') for u in ct_candidates]
    ref_sops = set(fm['ref_sop'] for fm in frame_meta if fm['ref_sop'] is not None)
    seg_z_arr = np.fromiter((fm['z'] for fm in frame_meta), dtype=np.float64)

    best = None
    fallback_best = None
    for ct_uid, ct_desc in ct_candidates:
        is_pre = 'PRE' in str(ct_desc).upper()
        ct_index = index_ct_series(SeriesRef('', ct_uid))  # dicom_root already baked
        overlap = ref_sops & set(ct_index.keys()) if ref_sops else set()
        for_match = bool(seg_for_uid) and any(
            m['for_uid'] == seg_for_uid for m in ct_index.values())
        ct_z_arr = np.array([m['z'] for m in ct_index.values()], dtype=np.float64)
        if ct_z_arr.size:
            diffs = np.min(np.abs(seg_z_arr[:, None] - ct_z_arr[None, :]), axis=1)
            z_overlap = int((diffs <= 1.25).sum())
        else:
            z_overlap = 0
        seg_z_count = len(seg_z_arr)
        z_overlap_frac = z_overlap / seg_z_count if seg_z_count else 0.0
        if len(overlap) > 0 and z_overlap_frac >= 0.5:
            rank_p = (len(overlap), z_overlap, not is_pre)
            if best is None or rank_p > (len(best[2]), best[3], best[4]):
                best = (ct_uid, ct_index, overlap, z_overlap, not is_pre)
        rank_f = (for_match, not is_pre, z_overlap)
        if fallback_best is None or rank_f > (fallback_best[4], fallback_best[5], fallback_best[3]):
            fallback_best = (ct_uid, ct_index, overlap, z_overlap, for_match, not is_pre)
    if best:
        return best[0], best[1], best[2]
    assert fallback_best and fallback_best[3] > 0, \
        f'No CT series matches SEG by z-coord across {len(ct_candidates)} candidates'
    return fallback_best[0], fallback_best[1], set()


def build_ct_volume(ct_index, acq_filter: Optional[int] = None):
    """Stack CT slices into a (Z, H, W) volume, z-sorted ascending.
    If acq_filter is given, only slices with that AcquisitionNumber are kept
    (used for HCC-TACE-Seg triphasic studies to isolate portal phase)."""
    rows = [(sop, m) for sop, m in ct_index.items()
            if acq_filter is None or m['acq'] == acq_filter]
    rows.sort(key=lambda r: r[1]['z'])
    ct_vol = np.stack([m['hu'] for _, m in rows], axis=0)
    z_coords = np.array([m['z'] for _, m in rows])
    pix = rows[0][1]['pix_spacing']
    return ct_vol, z_coords, pix


def build_segment_mask(frame_meta, z_coords, segments: Iterable[int],
                         z_tol_mm: float = 1.25):
    """Union of SEG frames whose segment ∈ `segments`, aligned to CT z_coords.

    Returns (Z, H, W) uint8 mask, or None if no frames matched.
    """
    segments = set(segments) if not isinstance(segments, int) else {segments}
    z_arr = np.asarray(z_coords)
    H = W = None
    slices = {}
    for fm in frame_meta:
        if fm['segment'] not in segments or not fm['mask'].any():
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


def crop_and_window(ct_vol, mass_mask, out_size: int = CROP_SIZE,
                     pad_lo: int = PAD_LOW, pad_hi: int = PAD_HIGH):
    """ST3DCN paper-faithful crop:
      1. Tight tumor bbox + asymmetric padding (+10, +9) per axis
      2. HU window [-160, 240] → uint8 [0, 255]
      3. Axis transpose (Z, Y, X) → (Y, X, Z) — paper's .npy convention
      4. padzero if max(shape) ≤ 70 else linear resize to 70³
      5. Normalize to [0, 1] float32
    """
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
    ct_crop = np.transpose(ct_crop, (1, 2, 0))
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


def portal_vein_median_hu(frame_meta, ct_index, portal_vein_seg: int):
    """Phase sanity check: median HU inside portal vein segment.
    Portal venous phase expected > ~120 HU (opacified vein); arterial <80."""
    vals = []
    for fm in frame_meta:
        if fm['segment'] != portal_vein_seg or not fm['mask'].any():
            continue
        if fm['ref_sop'] and fm['ref_sop'] in ct_index:
            hu = ct_index[fm['ref_sop']]['hu']
        else:
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


def identify_portal_acq_by_aorta(frame_meta, ct_index, aorta_seg: int):
    """For HCC-TACE-Seg where SEG has both mass and aorta segments and CT has
    multiple AcquisitionNumbers (arterial + portal), the portal acq has lower
    max aorta HU than the arterial acq.
    Returns (portal_acq_number, per_acq_stats_dict)."""
    from collections import defaultdict
    acq_max = defaultdict(list)
    for fm in frame_meta:
        if fm['segment'] != aorta_seg or not fm['mask'].any():
            continue
        for sop, meta in ct_index.items():
            if abs(meta['z'] - fm['z']) < 1e-3:
                acq_max[meta['acq']].append(meta['hu'][fm['mask']].max())
    if not acq_max:
        acqs = sorted(set(m['acq'] for m in ct_index.values()))
        return acqs[-1], None
    stats = {a: (np.percentile(v, 95), np.max(v)) for a, v in acq_max.items()}
    portal = min(stats, key=lambda a: stats[a][1])
    return portal, stats
