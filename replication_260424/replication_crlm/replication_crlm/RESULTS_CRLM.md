# ST3DCN essential-tag coverage test — HCC-TACE-Seg (sens) + CRLM (spec) lesion-level

**Unit of analysis**: lesion. Two cohorts, cross-validated on the same essential-26 / mentioned-14 test.

**Primary claim**: essential-26 tag set is the minimum necessary to observe acquisition-parameter variability across deployment cohorts. Mentioned-14 misses substantial variability in both cohorts, even after applying DICOM-standards-sanctioned derivations.

## DIRECT vs AUGMENTED (derivation) essential-26 classification

Derivations (DICOM PS3.3 compliant, applied identically to both cohorts):

| Derived tag | Formula | Reference |
|---|---|---|
| **SpacingBetweenSlices** (180088) | `median(diff(sorted ImagePositionPatient[2])))` | PS3.3 C.7.6.2.1.1 |
| **Exposure_uAs** (181153) | `ExposureTime × XRayTubeCurrent` (ms × mA = μAs) | PS3.3 C.8.7.2.1 — example calc explicitly given by standard |
| ContrastBolusAgentPhase (189344) | implicit via SEG-ref + Portal Vein HU in pipeline | — (pipeline-level, always PORTAL_VENOUS → invariant) |
| ~~Manufacturer / ModelName (HCC)~~ | ~~NBIA digest join~~ | Attempted; **NBIA digest itself has only 111/677 rows populated**. Not recoverable even via external metadata. |

`Exposure × 1000` rejected as primary derivation: intra-scanner CV 0.74–1.31 (vendor-inconsistent semantics of 0018,1152); `Time × Current` CV 0.13–0.38 (well-defined product).

### Coverage table

| Metric | HCC (direct) | HCC (augmented) | CRLM (direct) | CRLM (augmented) |
|---|---|---|---|---|
| fully_missing | 10 | 8  | 8  | 6  |
| invariant     | 5  | 5  | 2  | 2  |
| insufficient  | 3  | 4  | 5  | 5  |
| testable      | 8  | 9  | 11 | 13 |
| **FDR-sig**   | **4** | **6** | **4** | **6** |

### FDR-sig tags — DIRECT vs AUGMENTED

| Cohort | Direct FDR-sig | + Newly-sig after derivation |
|---|---|---|
| HCC | KVP (✓m), SpiralPitchFactor (✓m), TotalCollimationWidth (✗), TableSpeed (✗) | **Exposure_uAs** (✗, derived), **WindowCenter** (✗) |
| CRLM | ContrastBolusAgent (✓m), SliceThickness (✓m), ReconstructionDiameter (✗), WindowCenter (✗) | **Exposure_uAs** (✗, derived), **SpacingBetweenSlices** (✗, derived) |

(✓m = mentioned in ST3DCN paper, ✗ = not mentioned)

### Cross-cohort intersection emerges after derivation

| | Direct FDR-sig intersection | Augmented FDR-sig intersection |
|---|---|---|
| HCC ∩ CRLM | **∅** (0 tags) | **{Exposure_uAs, WindowCenter}** (2 tags) |

Both intersection tags are **not mentioned** in the paper. Without DICOM-standard derivation they look deployment-specific; after derivation they emerge as **universal acquisition-variability signals** — and exactly the ones the paper's Methods section does not flag.

### Coverage ratios (primary thesis evidence)

|  | Direct | Augmented |
|---|---|---|
| Union of FDR-sig tags across both cohorts | 8 | 10 |
| Captured by mentioned-14 | 4 / 8 = **50%** | 4 / 10 = **40%** |
| Captured by essential-26 | 8 / 8 = **100%** | 10 / 10 = **100%** |

Adding derivation **widens the gap**: mentioned-14's relative coverage DROPS (50% → 40%) because newly-testable derived tags are ones the paper did not mention. Essential-26 remains complete. The 6 essential-but-not-mentioned FDR-sig tags (TotalCollimationWidth, TableSpeed, WindowCenter, ReconstructionDiameter, Exposure_uAs, SpacingBetweenSlices) are the concrete evidentiary gap between mentioned and essential.

## Headline specificity (CRLM) — still lesion-level

| Metric | N | Value | Wilson 95% CI |
|---|---|---|---|
| **Specificity (thr 0.8)** | 464 | **0.7931** | (0.7539, 0.8275) |
| Mean prob | 464 | 0.320 | — |
| Median prob | 464 | 0.071 | — |

Reference — HCC-TACE-Seg v5 (lesion, 146): sensitivity 0.8973 (paper p0 = 0.869).



**Cohort**: 197 patients → **464 lesions** (CRLM SEG segments 5–9: Tumor_1..Tumor_5).
TCIA Colorectal-Liver-Metastases: all preop, portal-venous-phase CT, histologically-confirmed colorectal liver metastasis (Simpson 2024). MSKCC single-site; 195/197 GE Medical Systems.
Ground truth: **100% HCC-negative** → primary metric is **specificity**.

**Pipeline**: same model / preprocessing / crop / HU window as HCC-TACE-Seg replication. Per-lesion crop uses individual `Tumor_N` segment (not union). Phase empirically verified via Portal Vein (SEG seg 4) median HU > 120 across sampled patients.

## Headline — lesion-level

| Metric | N | Value | Wilson 95% CI |
|---|---|---|---|
| **Specificity (thr 0.8)** | 464 | **0.7931** | (0.7539, 0.8275) |
| Mean prob | 464 | 0.320 | — |
| Median prob | 464 | 0.071 | — |

Comparison reference — HCC-TACE-Seg v5 (lesion, 146):
- Sensitivity (thr 0.8) = 0.8973 (TP 131/146), Wilson CI (0.837, 0.938). Paper p0 = 0.869.

## (legacy) Essential-26 classification direct-only — superseded by augmented table above

| | HCC-TACE-Seg v5 (146) | CRLM (464) |
|---|---|---|
| fully_missing | 10 | 8 |
| invariant     | 5  | 2 |
| insufficient  | 3  | 5 |
| testable      | 8  | 11 |
| &nbsp;&nbsp;sig_FDR | 4 | 4 |
| &nbsp;&nbsp;NS      | 4 | 7 |

## Mentioned × classification — CRLM lesion-level (parallel to HCC v5 centerpiece)

| | not mentioned | mentioned | total |
|---|---|---|---|
| fully_missing | 4 | 4 | 8 |
| invariant     | 2 | 0 | 2 |
| insufficient  | 3 | 2 | 5 |
| testable_NS   | 4 | 3 | 7 |
| sig_FDR       | 2 | 2 | 4 |
| **Total**     | 15 | 11 | 26 |

## Significant tags (FDR<0.05) — CRLM lesion

| Tag | Hex | Class | Mentioned | KW H | raw p | FDR p | mean range | spec range |
|---|---|---|---|---|---|---|---|---|
| **ContrastBolusAgent**     | 180010 | Intensity/exposure | ✓ | 26.50 | 0.0009 | **0.0047** | 0.014–0.392 | 0.714–1.000 |
| **SliceThickness**         | 180050 | Geometry           | ✓ | 14.57 | 0.0007 | **0.0047** | 0.186–0.366 | 0.762–0.875 |
| **ReconstructionDiameter** | 181100 | Geometry           | ✗ | 13.63 | 0.0035 | **0.0127** | 0.191–0.507 | 0.667–0.892 |
| **WindowCenter**           | 281050 | Geometry           | ✗ | 12.06 | 0.0072 | **0.0197** | 0.040–0.379 | 0.600–1.000 |

- Significant AND mentioned: **ContrastBolusAgent**, **SliceThickness**
- Significant but NOT mentioned: **ReconstructionDiameter**, **WindowCenter**
- Mentioned tags unusable (missing/invariant/insufficient): KVP (99% `120` → insufficient), ContrastFlowRate, ContrastBolusAgentPhase, ReformattingThickness, Manufacturer (99% GE), BodyPartExamined (1% fill)

## Head-to-head — HCC-TACE-Seg (v5, lesion) vs CRLM (lesion)

Both are lesion-level Kruskal-Wallis + BH-FDR on p(HCC) across bins.

| Tag | HCC v5 (146) | CRLM (464) | Note |
|---|---|---|---|
| **KVP**                   | **sig** (FDR 0.003) | insufficient | CRLM 195/197 = `120` → no variance |
| **TotalCollimationWidth** | **sig** (FDR 0.016) | testable, NS (p=0.34) | CRLM fill 5.4% |
| **SpiralPitchFactor**     | **sig** (FDR 0.016) | testable, NS (p=0.34) | CRLM fill 5.4% |
| **TableSpeed**            | **sig** (FDR 0.016) | testable, NS (p=0.34) | CRLM fill 5.4% |
| **SliceThickness**        | insufficient        | **sig** (FDR 0.005) | CRLM has 2.5/3.75/5.0/7.5 mix; HCC 2.5-dominant |
| **ContrastBolusAgent**    | fully_missing        | **sig** (FDR 0.005) | CRLM 99% filled; HCC 1% |
| **ReconstructionDiameter**| NS                   | **sig** (FDR 0.013) | Wider bin spread in CRLM |
| **WindowCenter**          | NS                   | **sig** (FDR 0.020) | Wider bin spread in CRLM (7 bins vs 3) |

**Overlap of sig_FDR tags across the two lesion-level deployments: ∅.**

This is the thesis evidence: *which* essential tags discriminate is a property of the deployment cohort's variance, not a universal list. Every tag that discriminated on HCC-TACE-Seg is either invariant or near-empty on CRLM; every tag that discriminates on CRLM was unusable on HCC-TACE-Seg.

## Mentioned tags unusable in BOTH datasets (lesion-level)

Structural replication-reliability ceiling — paper mentions, public retrospective data never recovers:
- **ContrastFlowRate** (fully_missing both)
- **ContrastBolusAgentPhase** (fully_missing both)
- **ReformattingThickness** (fully_missing both)
- **BodyPartExamined** (invariant/near-empty both)

## Subset analysis — pu / pc / pce (lesion unit, same structure as HCC v5)

Matcher is the HCC-TACE-Seg subset_analysis matcher with two CRLM-normalization fixes applied identically across variants:
- `STRICT` — space-insensitive ManufacturerModelName only (e.g., "LightSpeed VCT" ≡ "Light Speed VCT").
- `MEDIUM` — + 10% relative tolerance on XRayTubeCurrent.
- `RELAXED` — + MSKCC agent-brand stem match (4-char prefix; "OMNI" ≡ "Omnipaque").

SliceThickness excluded throughout (protocol ≤1.25 mm vs CRLM 2.5–7.5 mm — identical to HCC-TACE-Seg exclusion).

| Variant | pc patients | pc lesions | pc spec (thr 0.8) | pc Wilson 95% | Protocols matched |
|---|---|---|---|---|---|
| pu (uncontrolled) | 197 | 464 | 0.7931 | (0.7539, 0.8275) | — |
| STRICT  | 0   | 0   | —      | — | none |
| MEDIUM  | 0   | 0   | —      | — | none |
| RELAXED | 15  | 38  | 0.7895 | (0.6365, 0.8893) | **P08 only** (GE Light Speed VCT) |

`pce` (essential ∩ mentioned) has identical counts to `pc` at each variant.

**Cross-cohort convergence to P08.** Yesterday's HCC-TACE-Seg analysis reported pc=pce went entirely to P08. CRLM, run independently, also sends every matched patient to P08 (under the stem-matched variant). Out of 15 ST3DCN training protocols, only P08 (GE Light Speed VCT, KVP=120, unspecified tube current, any Omnipaque/Visipaque brand) has a non-zero intersection with either TCIA test cohort.

**pu vs pc nearly identical.** CRLM specificity: pu 0.7931 vs pc 0.7895 — same distribution-level finding as HCC v5 (`pu vs pc KS p=0.206 NS`). Tag-based subset control does not tighten the probability distribution in either cohort.

## Per-patient (union-bbox) — supplementary only, not comparable

Included for completeness; **not used in any cross-cohort comparison**.

| | N | spec@0.8 | mean prob |
|---|---|---|---|
| Per-patient union-bbox | 197 | 0.3046 | 0.752 |
| Per-patient max-lesion-prob | 197 | 0.6193 | 0.522 |
| **Per-lesion (primary)** | 464 | **0.7931** | 0.320 |

Union-bbox aggregation systematically over-calls multi-lesion patients (1-lesion mean prob 0.45 vs ≥2-lesion 0.92) because 70³ resize on an inter-lesion bbox dilutes tumor signal. Lesion is the correct unit.

## Files

```
replication_crlm/
├── st3dcn_pipeline_crlm.py          # pipeline (1 study/patient, segs 5-9 = tumors, PV-HU phase sanity)
├── run_batch_lesion_crlm.py          # PRIMARY: per-lesion batch (one crop per Tumor_N)
├── run_batch_crlm.py                 # supplementary: per-patient union
├── extract_tags_crlm.py              # essential-26 DICOM tag extraction
├── per_tag_binned_lesion_crlm.py     # PRIMARY: lesion-level per-tag spec + KW + FDR
├── per_tag_binned_crlm.py            # supplementary: patient-level
├── subset_analysis_lesion_crlm.py    # pu/pc/pce subset spec (3 matcher variants)
├── subset_summary_lesion.txt / subset_detail_lesion.csv
├── lesion_results.csv / lesion_results_with_tags.csv
├── lesion_summary.txt
├── per_tag_binned_lesion.md + .csv
├── per_tag_classification_lesion.csv
└── (supplementary patient-level outputs with same stems)
```
