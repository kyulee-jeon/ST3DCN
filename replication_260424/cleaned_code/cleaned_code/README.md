# ST3DCN replication + essential-tag heterogeneity analyzer

Cohort-agnostic framework that:
1. Replicates ST3DCN (Yu et al., *JHEP Reports* 2025) preprocessing + inference pixel-exact against the paper's demo sample.
2. Runs per-lesion inference on any TCIA-style DICOM collection with SEG masks.
3. Tests whether each of 26 *essential* CT DICOM tags expresses acquisition-parameter heterogeneity on that cohort (Kruskal-Wallis + Benjamini-Hochberg FDR).
4. Applies two DICOM-standard derivations (`SpacingBetweenSlices` ← `ImagePositionPatient`; `Exposure_uAs` ← `ExposureTime × XRayTubeCurrent`) to rescue otherwise-missing essential tags.
5. Produces a cross-cohort essential-vs-mentioned report.

Validated on two TCIA cohorts:
- **HCC-TACE-Seg** (105 HCC+ patients, 146 lesions) — sensitivity path.
- **Colorectal-Liver-Metastases (CRLM)** (197 HCC− patients, 464 lesions) — specificity path.

---

## Directory layout

```
cleaned_code/
├── pipeline/
│   ├── preprocessing.py        Paper-faithful SEG→CT align, 70³ crop, HU [-160,240]
│   ├── dataset_adapter.py       Per-cohort rules (segment layout, study selection)
│   ├── lesion_pipeline.py       Patient → per-lesion crops
│   └── inference.py             ST3DCN model loading + predict
├── analysis/
│   ├── essential_tags.py        Essential-26 list (+ mentioned flag)
│   ├── extract_tags.py          DICOM tag extraction from CT series
│   ├── derive_tags.py           DICOM-standard derivations
│   ├── binners.py               Tag binning strategies
│   ├── per_tag_analysis.py      Kruskal-Wallis + BH-FDR per tag
│   └── combined_report.py       Cross-cohort essential-vs-mentioned report
├── configs/
│   ├── hcc_tace_seg.yaml
│   ├── crlm.yaml
│   └── new_dataset_template.yaml
└── runners/
    ├── run_lesion_inference.py
    ├── run_analysis.py
    └── build_cross_cohort_report.py
```

---

## Input specification

### 1. DICOM data layout
Each series must live at `<dicom_root>/<SeriesInstanceUID>/*.dcm`. Required:
- **CT** series (Modality=CT): axial slices with populated `ImagePositionPatient`, `PixelSpacing`, `SliceThickness`.
- **SEG** series (Modality=SEG): DICOM-SEG format with `PerFrameFunctionalGroupsSequence` referencing source CT slices via `DerivationImageSequence.SourceImageSequence.ReferencedSOPInstanceUID`.

### 2. Manifest
A table with one row per series. Columns:
- `PatientID`, `StudyInstanceUID`, `SeriesInstanceUID`, `Modality`, `SeriesDescription` (may be empty).

Supported formats:
- CSV (`manifest_csv` config key)
- JSON list-of-dicts (`manifest_json` config key)

### 3. (Optional) Longitudinal offset / clinical metadata
For multi-study cohorts (e.g., HCC-TACE-Seg has baseline + follow-up CTs):
- NBIA digest Excel with `Longitudinal Temporal Offset From Event` column.
- Clinical Excel for patient roster.

### 4. (Optional) NBIA digest for Manufacturer / ModelName rescue
If DICOM-level Manufacturer fill is low, provide the NBIA digest Excel — `derive_tags.py` will join by `Series Instance UID`.

### 5. Segmentation layout (cohort-specific)
Each cohort YAML declares:
- `tumor_segments`: integer list of SEG `ReferencedSegmentNumber` values representing tumor(s).
  - CRLM: `[5,6,7,8,9]` (Tumor_1..Tumor_5)
  - HCC-TACE-Seg: `[2]` (Mass)
- `portal_vein_segment`: (optional) used for portal-venous-phase sanity check (median HU > 120 expected).
- `aorta_segment`: (optional) used when CT series bundles multiple acquisitions and portal phase must be picked by HU level.

---

## Preprocessing pipeline (paper-faithful)

**Validation**: the below sequence was empirically verified pixel-exact (|Δ|=0 on all 10/10 observations) against the ST3DCN repo's demo pre-cropped `.npy` files for `ID_0848` (HCC) and `QEH032` (non-HCC).

1. **SEG→CT alignment**
   - Primary: `ReferencedSOPInstanceUID` in SEG frame → match to CT `SOPInstanceUID`. Require ≥50% z-coordinate overlap.
   - Fallback: `FrameOfReferenceUID` + z-coordinate nearest-neighbor (tol 1.25mm).
2. **Phase selection** (cohort-specific)
   - Single-phase cohorts (CRLM): use all slices of the single CT series.
   - Multi-phase cohorts (HCC-TACE-Seg): discriminate arterial vs portal by aorta-segment max HU across `AcquisitionNumber` groups; portal = lower max HU.
3. **Volume assembly**
   - Stack CT slices in ascending z-order; preserve original resolution.
4. **Tumor mask**
   - Per-lesion: for multi-segment cohorts, each tumor segment = one lesion. For single-segment cohorts, connected-component labeling on the union mask.
5. **Crop**
   - Tight bbox of tumor mask; asymmetric padding `+10` (low) / `+9` (high) per axis.
6. **HU windowing**
   - Clip to [-160, 240] HU, scale to uint8 [0, 255].
7. **Axis transpose**
   - (Z, Y, X) → (Y, X, Z) to match paper's `.npy` convention.
8. **Resize / pad**
   - If any axis > 70: linear zoom to 70³.
   - Otherwise: zero-pad centered to 70³.
9. **Normalize**
   - Divide by 255 → float32 in [0, 1].

**Scope decision**: **portovenous phase only, pre-treatment only, lesion-level analysis**. Justifications:
- Paper Table S4F validates portal-venous-only path (AUC 0.912 ≈ triphasic 0.919).
- SEG masks are drawn on portal venous phase in both benchmark cohorts (Morshid 2019 for HCC-TACE-Seg; Simpson 2024 for CRLM).
- Avoids mask-to-arterial projection and DICOM phase-metadata gaps (`ContrastBolusStartTime`, `ContrastBolusAgentPhase` are 99-100% null in both cohorts).
- Lesion-level avoids the multi-lesion union-bbox dilution artifact (see CRLM results: per-patient union-bbox specificity 0.30 vs per-lesion 0.79).

---

## Output specification

### `run_lesion_inference.py` outputs
- `lesion_results.csv`: one row per lesion, columns: `patient_id, lesion_idx, segment_number, segment_label, voxel_count, prob, pred, phase_sanity, ct_shape, bbox, seg_uid, ct_uid, study_uid, pix_spacing_mm, slice_thick, n_lesions_patient, status, error`.
- `patient_results.csv`: one row per patient, aggregating min/max/mean prob across lesions.
- `lesion_summary.txt`: overall sensitivity (HCC+ cohort) or specificity (HCC− cohort) with Wilson 95% CI.

### `run_analysis.py` outputs
- `lesion_results_with_tags.csv`: lesion rows + 26 essential DICOM tags + extras.
- `lesion_results_with_derived.csv`: above + `SpacingBetweenSlices_derived`, `Exposure_uAs_derived`, optionally `Manufacturer_derived` / `ManufacturerModelName_derived`.
- `per_tag_classification.csv`: one row per essential tag, columns: `tag, hex, class, mentioned_in_paper, derived_applied, fill_rate, n_unique_bins, n_testable_bins, status ∈ {fully_missing, invariant, insufficient, testable}, kw_H, kw_p, kw_p_fdr, significant_fdr, mean_min, mean_max`.
- `per_tag_binned.csv`: per-tag per-bin breakdown with metric (sens or spec), Wilson CI, mean/median prob.
- `per_tag_report.md`: readable writeup.

### `build_cross_cohort_report.py` outputs
- A single markdown document with:
  - Observations 1/2/3 (mentioned coverage, cross-cohort intersection, mentioned unusability).
  - Essential-26 classification side-by-side across cohorts.
  - Per-tag bin tables with each cohort's n / metric / mean prob.

---

## Adding a new dataset

1. **Copy** `configs/new_dataset_template.yaml` → `configs/<your_name>.yaml`. Fill in paths and segment layout.
2. **Subclass** `DatasetAdapter` in `pipeline/dataset_adapter.py` if the manifest format or series-selection logic differs from CRLM / HCC-TACE-Seg. Register in `ADAPTERS` dict.
3. **Verify phase**: if SEG has a portal-vein segment, set `portal_vein_segment` and confirm median HU > 120 on a few patients. Otherwise, document the phase source.
4. **Run**
   ```bash
   python runners/run_lesion_inference.py --cohort <your_name> --out-dir out_<your_name>
   python runners/run_analysis.py          --cohort <your_name> --out-dir out_<your_name> \\
       --lesion-csv out_<your_name>/lesion_results.csv \\
       --patient-csv out_<your_name>/patient_results.csv
   ```
5. **Compare to existing cohorts**
   ```bash
   python runners/build_cross_cohort_report.py --hcc-dir out_hcc --crlm-dir out_crlm \\
       --out results_cross_cohort.md
   ```

---

## Environment

- Python 3.10 (`hcc_venv`)
- TensorFlow 2.11 + cuDNN 8.1 (CUDA 11.2 at `/usr/local/cuda-11.2/lib64`)
- `pydicom`, `scipy`, `scikit-image`, `pandas`, `numpy`
- ST3DCN repo cloned at `/home/ubuntu/hcc_workspace/ST3DCN` (for `ST3DCN_Utils`)

Run commands must set `LD_LIBRARY_PATH=/usr/local/cuda-11.2/lib64:$LD_LIBRARY_PATH` to pick up the correct cuDNN.

---

## Thresholds and statistical decisions

| Decision | Value | Reference |
|---|---|---|
| Binary HCC prediction threshold | p ≥ 0.8 | ST3DCN paper Methods |
| Minimum lesion voxels (skip noise) | 50 | Empirical; Morshid 2019 SEG artifacts |
| Minimum n per bin for Kruskal-Wallis | 5 | Standard nonparametric sample-size floor |
| FDR α | 0.05 | Benjamini-Hochberg |
| z-coordinate tolerance for SEG↔CT alignment | 1.25 mm | ≤ smallest plausible slice thickness |
| Portal vein HU sanity threshold | 120 HU | Typical portal-venous opacification |

---

## License & citation

- Original ST3DCN: Yu et al., *JHEP Reports* 2025; repo `HKUMedicineLiverAI/ST3DCN`.
- HCC-TACE-Seg: TCIA, Morshid 2019.
- CRLM: TCIA, Simpson 2024.
- This replication framework is internal research code. If useful to you, cite the corresponding protocol2dicom abstract when available.
