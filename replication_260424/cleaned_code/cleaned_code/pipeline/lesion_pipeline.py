"""High-level per-lesion processing: given an adapter and a patient id,
produce one crop per tumor lesion ready for ST3DCN inference.

Lesion definition:
  - If adapter exposes multiple tumor segment numbers (e.g., CRLM:
    Tumor_1..Tumor_5), each segment = one lesion.
  - Otherwise fall back to connected components on the union tumor mask
    (e.g., HCC-TACE-Seg has a single Mass segment with possibly multiple CCs).
"""
from typing import Optional

from . import preprocessing as pp
from .dataset_adapter import DatasetAdapter, ProcessedCase
from .preprocessing import SeriesRef


def process_patient_lesions(adapter: DatasetAdapter, patient_id: str,
                              min_voxels: int = 50) -> ProcessedCase:
    """Process a single patient → one crop per tumor lesion.

    Flow:
      1. adapter.resolve_series  →  CT uid(s), SEG uid, study uid
      2. Read SEG, build frame_meta
      3. If multi-CT candidates, pick mask-aligned CT (preprocessing.pick_mask_aligned_ct)
      4. adapter.acq_filter(...)  →  None (all slices) or portal AcqNum
      5. Build CT volume (z-sorted), portal_vein_hu sanity (optional)
      6. If adapter.tumor_segments has >1 element → each segment = one lesion
         Else → connected-component decomposition on union mask
      7. crop_and_window(70³, HU [-160,240]) per lesion
    """
    ct_ref, seg_uid, study_uid = adapter.resolve_series(patient_id)

    seg_dataset, _, frame_meta, seg_labels = pp.read_seg(
        SeriesRef(adapter.dicom_root, seg_uid))

    # CT selection (single vs multi-candidate)
    if isinstance(ct_ref, list):
        # Multi-candidate list of (uid, desc) — pick mask-aligned
        # Each index_ct_series call needs dicom_root; patch ct_candidates
        ct_candidates = [(u, d) for u, d in ct_ref]
        # Use a local pick_mask_aligned_ct that knows dicom_root
        seg_for_uid = getattr(seg_dataset, 'FrameOfReferenceUID', None)
        ct_uid, ct_index, _ = _pick_mask_aligned_ct_rooted(
            adapter.dicom_root, frame_meta, ct_candidates, seg_for_uid)
    else:
        ct_uid = ct_ref
        ct_index = pp.index_ct_series(SeriesRef(adapter.dicom_root, ct_uid))
        n_match, n_ref = pp.validate_sop_linkage(frame_meta, ct_index)
        assert n_match == n_ref or n_ref == 0, \
            f'{patient_id}: SEG refs {n_match}/{n_ref} matched CT'

    # Acq filter (None for CRLM; portal-acq pick for HCC-TACE-Seg)
    acq = adapter.acq_filter(patient_id, frame_meta, ct_index)

    ct_vol, z_coords, pix = pp.build_ct_volume(ct_index, acq_filter=acq)
    slice_thick = list(ct_index.values())[0]['slice_thick']

    # Phase sanity (optional)
    phase = {}
    if adapter.portal_vein_segment:
        pv_hu = pp.portal_vein_median_hu(frame_meta, ct_index,
                                            adapter.portal_vein_segment)
        phase['portal_vein_median_hu'] = pv_hu

    # Per-lesion crops
    lesions = []
    tumor_segs = list(adapter.tumor_segments)
    if len(tumor_segs) > 1:
        # explicit per-segment lesions (e.g., CRLM Tumor_1..Tumor_5)
        for seg_num in sorted(tumor_segs):
            m = pp.build_segment_mask(frame_meta, z_coords, [seg_num])
            if m is None or int(m.sum()) < min_voxels:
                continue
            crop, bbox = pp.crop_and_window(ct_vol, m)
            lesions.append({
                'lesion_idx': len(lesions),
                'segment_number': seg_num,
                'segment_label': seg_labels.get(seg_num, f'Seg{seg_num}'),
                'voxel_count': int(m.sum()),
                'bbox': bbox,
                'crop': crop,
            })
    else:
        # single mass segment → connected components on union mask
        from skimage.measure import label as cc_label, regionprops
        mass_mask = pp.build_segment_mask(frame_meta, z_coords, tumor_segs)
        assert mass_mask is not None and mass_mask.sum() > 0, \
            f'{patient_id}: empty tumor mask'
        labeled = cc_label(mass_mask > 0, connectivity=1)
        props = sorted(regionprops(labeled), key=lambda p: p.area, reverse=True)
        for p in props:
            if p.area < min_voxels:
                continue
            comp_mask = (labeled == p.label).astype('uint8')
            crop, bbox = pp.crop_and_window(ct_vol, comp_mask)
            lesions.append({
                'lesion_idx': len(lesions),
                'segment_number': tumor_segs[0],
                'segment_label': seg_labels.get(tumor_segs[0], 'Mass'),
                'voxel_count': int(p.area),
                'bbox': bbox,
                'crop': crop,
            })

    assert lesions, f'{patient_id}: no tumor lesions >= {min_voxels} voxels'
    return ProcessedCase(
        patient_id=patient_id,
        ct_uid=ct_uid,
        seg_uid=seg_uid,
        study_uid=study_uid,
        ct_shape=ct_vol.shape,
        pix_spacing_mm=pix,
        slice_thick_mm=slice_thick,
        phase_sanity=phase,
        lesions=lesions,
    )


def _pick_mask_aligned_ct_rooted(dicom_root, frame_meta, ct_candidates,
                                    seg_for_uid):
    """Like preprocessing.pick_mask_aligned_ct but knows dicom_root
    (the function in preprocessing.py takes a SeriesRef)."""
    import numpy as np
    ref_sops = set(fm['ref_sop'] for fm in frame_meta if fm['ref_sop'] is not None)
    seg_z_arr = np.fromiter((fm['z'] for fm in frame_meta), dtype=np.float64)
    best = None; fallback_best = None
    for ct_uid, ct_desc in ct_candidates:
        is_pre = 'PRE' in str(ct_desc).upper()
        ct_index = pp.index_ct_series(SeriesRef(dicom_root, ct_uid))
        overlap = ref_sops & set(ct_index.keys()) if ref_sops else set()
        for_match = bool(seg_for_uid) and any(
            m['for_uid'] == seg_for_uid for m in ct_index.values())
        ct_z_arr = np.array([m['z'] for m in ct_index.values()], dtype=np.float64)
        if ct_z_arr.size:
            diffs = np.min(np.abs(seg_z_arr[:, None] - ct_z_arr[None, :]), axis=1)
            z_overlap = int((diffs <= 1.25).sum())
        else:
            z_overlap = 0
        frac = z_overlap / len(seg_z_arr) if len(seg_z_arr) else 0.0
        if len(overlap) > 0 and frac >= 0.5:
            rank_p = (len(overlap), z_overlap, not is_pre)
            if best is None or rank_p > (len(best[2]), best[3], best[4]):
                best = (ct_uid, ct_index, overlap, z_overlap, not is_pre)
        rank_f = (for_match, not is_pre, z_overlap)
        if fallback_best is None or rank_f > (fallback_best[4], fallback_best[5], fallback_best[3]):
            fallback_best = (ct_uid, ct_index, overlap, z_overlap, for_match, not is_pre)
    if best:
        return best[0], best[1], best[2]
    assert fallback_best and fallback_best[3] > 0, 'No CT series matches SEG'
    return fallback_best[0], fallback_best[1], set()
