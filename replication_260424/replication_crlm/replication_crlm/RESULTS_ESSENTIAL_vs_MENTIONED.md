# Essential vs Mentioned DICOM Tags — ST3DCN replication on HCC-TACE-Seg + CRLM

**Cohorts**:
- **HCC-TACE-Seg** — 105 HCC+ patients / 146 lesions (TCIA). ST3DCN sensitivity path.
- **CRLM** — 197 HCC− patients / 464 lesions (TCIA, Colorectal-Liver-Metastases). ST3DCN specificity path.

Same pipeline, same preprocessing (70³ crop, HU window [-160,240]), same threshold 0.8 on p(HCC). Lesion-level analysis in both. Per-essential-tag heterogeneity tested by Kruskal-Wallis on p(HCC), with **BH-FDR α = 0.05** across testable tags.

**Essential set**: 26 CT tags (protocol2dicom). **Mentioned set**: 11 of those 26 (in ST3DCN paper Methods / Table S1). Table S1 also mentions 3 non-CT attributes (AcquisitionSequence / InjectionDelay / AgentSequence) which are not value-typed DICOM fields; excluded from per-tag testing.

**Two DICOM-standard derivations applied** to rescue fields that are direct-populated <5%:
- `SpacingBetweenSlices (0018,0088)` ← `median(diff(sorted ImagePositionPatient[2]))` — DICOM PS3.3 §C.7.6.2.1.1
- `Exposure_uAs (0018,1153)` ← `ExposureTime × XRayTubeCurrent` — DICOM PS3.3 §C.8.7.2.1 (standard-listed example calc; chosen over `Exposure×1000` because of 0018,1152 vendor-inconsistent semantics, intra-scanner CV 0.74-1.31 vs 0.13-0.38)

Non-recoverable despite effort: Manufacturer / ModelName on HCC-TACE-Seg (NBIA digest itself has only 111/677 rows populated).

## 0. Headline — thesis numbers

| | Direct-only | + DICOM-standard derivation |
|---|---|---|
| HCC sig tags | 4 (2 mentioned + 2 not) | **6** (2 mentioned + 4 not) |
| CRLM sig tags | 4 (2 mentioned + 2 not) | **5** (2 mentioned + 3 not) |
| Union sig tags (both cohorts) | 8 | **10** |
| HCC ∩ CRLM sig tags | **∅** | **{ WindowCenter }** |
| Mentioned-11 coverage of union | 4/8 = **50%** | 4/10 = **40%** |
| Essential-26 coverage of union | 8/8 = **100%** | 10/10 = **100%** |

Derivation widens the essential-vs-mentioned gap — newly-testable derived tags are all *not mentioned* in the paper.

## 1. Observation 1 — Mentioned-11 captures only 2/6 (HCC) and 2/5 (CRLM) of FDR-sig tags

| Cohort | FDR-sig total | mentioned in paper | NOT mentioned in paper |
|---|---|---|---|
| HCC-TACE-Seg (all HCC+) | 6 | KVP, SpiralPitchFactor (**2/6**) | Exposure_uAs, TableSpeed, TotalCollimationWidth, WindowCenter |
| CRLM (all HCC−)         | 5 | ContrastBolusAgent, SliceThickness (**2/5**) | ReconstructionDiameter, SpacingBetweenSlices, WindowCenter |

Across both cohorts, 4 of 10 FDR-sig tags are mentioned; the remaining 6 (Exposure_uAs, ReconstructionDiameter, SpacingBetweenSlices, TableSpeed, TotalCollimationWidth, WindowCenter) are essential-set-only discoveries that the paper did not cite.

## 2. Observation 2 — Cross-cohort intersection: ∅ (direct) → {WindowCenter} (derived)

- **Direct-only**: HCC ∩ CRLM = ∅
- **After derivation**: HCC ∩ CRLM = {WindowCenter}
  - `WindowCenter`: **not mentioned** in paper; essential-26 ∈

Which tags express variability is deployment-specific. A single-cohort study (like the ST3DCN paper) cannot enumerate them a priori — it can only find the tags visible in its own distribution. Consequences:
1. The safety net needs to be wider than any single cohort's sig set → essential-26 exists for this reason.
2. Using just the paper's mentioned-11 on HCC-TACE-Seg would miss 4 tags; on CRLM, 3 tags — all in essential-26 but not in Methods.
3. Derivation reveals tags that would be invisible in naive fill-only analysis. In particular, Exposure_uAs is direct-0% in both cohorts but derivable via the DICOM-sanctioned formula — turning a "fully missing" row into a testable (and in HCC, sig) one.

## 3. Observation 3 — Most mentioned tags are simply unusable

Not every mentioned tag is populated. The paper mentioning a tag in Methods does not guarantee the tag is observable on deployment DICOMs. *Even after derivation attempts, most mentioned tags remain non-testable.*

| Cohort | Mentioned-11 testable | Mentioned-11 FDR-sig | Mentioned-11 unusable |
|---|---|---|---|
| HCC-TACE-Seg | 3 | 2 | **8** |
| CRLM         | 5 | 2 | **6** |

Unusable mentioned tags per cohort:

- **HCC**: ContrastBolusAgent (fully_missing), ContrastFlowRate (fully_missing), ContrastBolusAgentPhase (fully_missing), SliceThickness (insufficient), ReformattingThickness (fully_missing), ManufacturerModelName (fully_missing), Manufacturer (fully_missing), BodyPartExamined (invariant)
- **CRLM**: KVP (insufficient), ContrastFlowRate (fully_missing), ContrastBolusAgentPhase (fully_missing), ReformattingThickness (fully_missing), Manufacturer (insufficient), BodyPartExamined (fully_missing)

Mentioned-in-BOTH-cohorts-unusable (structural paper-mentioned-vs-data gap):

- `BodyPartExamined` — HCC invariant, CRLM fully_missing
- `ContrastBolusAgentPhase` — HCC fully_missing, CRLM fully_missing
- `ContrastFlowRate` — HCC fully_missing, CRLM fully_missing
- `Manufacturer` — HCC fully_missing, CRLM insufficient
- `ReformattingThickness` — HCC fully_missing, CRLM fully_missing

Essential-26 as a whole nearly doubles observable variability:

| Cohort | Essential-26 testable | Essential-26 FDR-sig | Mentioned-11 testable | Mentioned-11 sig |
|---|---|---|---|---|
| HCC-TACE-Seg | 9 | 6 | 3 | 2 |
| CRLM         | 13 | 5 | 5 | 2 |

## 4. Full classification — all 26 essential tags (HCC ‖ CRLM)

| Tag | Hex | Class | Mentioned | HCC fill% | HCC status | HCC FDR | HCC | CRLM fill% | CRLM status | CRLM FDR | CRLM |
|---|---|---|---|---|---|---|---|---|---|---|---|
| KVP | 180060 | I | ✓ | 100.0 | testable | 0.0021 | **★** | 100.0 | insufficient | — |  |
| XRayTubeCurrent | 181151 | I | ✓ | 99.3 | testable | 0.6483 |  | 44.4 | testable | 0.8182 |  |
| ConvolutionKernel | 181210 | I | ✗ | 99.3 | insufficient | — |  | 44.4 | insufficient | — |  |
| ContrastBolusAgent | 180010 | I | ✓ | 1.4 | fully_missing | — |  | 99.1 | testable | 0.0056 | **★** |
| CTDIvol | 189345 | I | ✗ | 0.7 | fully_missing | — |  | 0.0 | fully_missing | — |  |
| ContrastFlowRate | 181046 | I | ✓ | 0.0 | fully_missing | — |  | 0.0 | fully_missing | — |  |
| Exposure_uAs | 181153 | I | ✗ | 98.6 | testable | 0.0021 | **★**◆ | 44.4 | testable | 0.4668 | ◆ |
| RevolutionTime | 189305 | I | ✗ | 52.7 | invariant | — |  | 5.4 | insufficient | — |  |
| ReconstructionAlgorithm | 189315 | I | ✗ | 0.0 | fully_missing | — |  | 0.0 | fully_missing | — |  |
| ContrastBolusAgentPhase | 189344 | I | ✓ | 0.0 | fully_missing | — |  | 0.0 | fully_missing | — |  |
| SliceThickness | 180050 | G | ✓ | 100.0 | insufficient | — |  | 100.0 | testable | 0.0056 | **★** |
| Rows | 280010 | G | ✗ | 100.0 | invariant | — |  | 100.0 | invariant | — |  |
| Columns | 280011 | G | ✗ | 100.0 | invariant | — |  | 100.0 | invariant | — |  |
| PixelSpacing | 280030 | G | ✗ | 100.0 | testable | 0.1974 |  | 100.0 | testable | 0.4668 |  |
| ReconstructionDiameter | 181100 | G | ✗ | 100.0 | testable | 0.2182 |  | 100.0 | testable | 0.0150 | **★** |
| WindowCenter | 281050 | G | ✗ | 100.0 | testable | 0.0488 | **★** | 100.0 | testable | 0.0233 | **★** |
| WindowWidth | 281051 | G | ✗ | 100.0 | insufficient | — |  | 100.0 | testable | 0.1114 |  |
| TotalCollimationWidth | 189307 | G | ✗ | 52.7 | testable | 0.0140 | **★** | 5.4 | testable | 0.9535 |  |
| SpiralPitchFactor | 189311 | G | ✓ | 52.7 | testable | 0.0140 | **★** | 5.4 | testable | 0.9535 |  |
| TableSpeed | 189309 | G | ✗ | 52.7 | testable | 0.0140 | **★** | 5.4 | testable | 0.9535 |  |
| SpacingBetweenSlices | 180088 | G | ✗ | 100.0 | insufficient | — | ◆ | 100.0 | testable | 0.0274 | **★**◆ |
| ReformattingThickness | 720512 | G | ✓ | 0.0 | fully_missing | — |  | 0.0 | fully_missing | — |  |
| ManufacturerModelName | 081090 | D | ✓ | 1.4 | fully_missing | — |  | 100.0 | testable | 0.9535 |  |
| Manufacturer | 080070 | D | ✓ | 1.4 | fully_missing | — |  | 100.0 | insufficient | — |  |
| BodyPartExamined | 180015 | A | ✓ | 100.0 | invariant | — |  | 0.6 | fully_missing | — |  |
| PatientPosition | 185100 | A | ✗ | 100.0 | invariant | — |  | 100.0 | insufficient | — |  |

★ = FDR-significant (<0.05) ; ◆ = derivation-augmented (SpacingBetweenSlices ← IPP; Exposure_uAs ← ExposureTime × XRayTubeCurrent)

## 5. Per-tag bin tables — HCC sens / mean-p  vs  CRLM spec / mean-p

Tables grouped by class, in order of essential spec.


### ▸ Intensity/exposure

### KVP `(180060)` — Intensity/exposure · ✓ mentioned

- **HCC**  · fill 100.0% · status **★ FDR-sig** · FDR 0.0021
- **CRLM** · fill 100.0% · status insufficient · FDR —

| bin | HCC n | HCC sens | HCC mean p | CRLM n | CRLM spec | CRLM mean p |
|---|---|---|---|---|---|---|
| `100` | 10 | 0.700 | 0.757 | — | — | — |
| `120` | 79 | 0.797 | 0.786 | — | — | — |
| `140` | 57 | 0.877 | 0.885 | — | — | — |
| `120.0` | — | — | — | 461 | 0.794 | 0.320 |
| `130.0` | — | — | — | 1 | 1.000 | 0.014 |
| `140.0` | — | — | — | 2 | 0.500 | 0.485 |

### XRayTubeCurrent `(181151)` — Intensity/exposure · ✓ mentioned

- **HCC**  · fill 99.3% · status testable · FDR 0.6483
- **CRLM** · fill 44.4% · status testable · FDR 0.8182

| bin | HCC n | HCC sens | HCC mean p | CRLM n | CRLM spec | CRLM mean p |
|---|---|---|---|---|---|---|
| `[0, 200)` | 6 | 1.000 | 0.924 | 13 | 0.846 | 0.238 |
| `[200, 300)` | 22 | 0.909 | 0.875 | 56 | 0.804 | 0.312 |
| `[300, 400)` | 91 | 0.791 | 0.798 | 62 | 0.726 | 0.329 |
| `[400, 500)` | 4 | 1.000 | 0.974 | 55 | 0.891 | 0.221 |
| `[500, 10000)` | 22 | 0.773 | 0.810 | 20 | 0.750 | 0.386 |
| `(missing)` | 1 | 1.000 | 0.968 | 258 | 0.787 | 0.340 |

### ConvolutionKernel `(181210)` — Intensity/exposure · ✗ not mentioned

- **HCC**  · fill 99.3% · status insufficient · FDR —
- **CRLM** · fill 44.4% · status insufficient · FDR —

| bin | HCC n | HCC sens | HCC mean p | CRLM n | CRLM spec | CRLM mean p |
|---|---|---|---|---|---|---|
| `(missing)` | 1 | 1.000 | 0.968 | 258 | 0.787 | 0.340 |
| `B` | 1 | 1.000 | 0.998 | — | — | — |
| `SOFT` | 4 | 0.500 | 0.424 | — | — | — |
| `STANDARD` | 140 | 0.829 | 0.832 | 203 | 0.803 | 0.295 |
| `AA-1` | — | — | — | 2 | 0.500 | 0.485 |
| `B40` | — | — | — | 1 | 1.000 | 0.014 |

### ContrastBolusAgent `(180010)` — Intensity/exposure · ✓ mentioned

- **HCC**  · fill 1.4% · status fully_missing · FDR —
- **CRLM** · fill 99.1% · status **★ FDR-sig** · FDR 0.0056

| bin | HCC n | HCC sens | HCC mean p | CRLM n | CRLM spec | CRLM mean p |
|---|---|---|---|---|---|---|
| `(missing)` | 144 | 0.826 | 0.827 | 4 | 1.000 | 0.339 |
| `ORAL & 125CC` | 2 | 0.500 | 0.494 | — | — | — |
| `&` | — | — | — | 2 | 1.000 | 0.378 |
| `APPLIED` | — | — | — | 1 | 1.000 | 0.014 |
| `BARIUM  & OMNI` | — | — | — | 1 | 1.000 | 0.003 |
| `BARIUM & NON IONIC` | — | — | — | 3 | 0.667 | 0.474 |
| `BARIUM & OMNI` | — | — | — | 53 | 0.736 | 0.320 |
| `CONT` | — | — | — | 2 | 0.500 | 0.485 |
| `GASTRO & NON IONIC` | — | — | — | 4 | 1.000 | 0.276 |
| `GASTRO & OMNI` | — | — | — | 98 | 0.867 | 0.265 |
| `GASTRO & OMNI 300` | — | — | — | 5 | 0.800 | 0.349 |
| `GASTRO & OMNIPAQUE` | — | — | — | 2 | 1.000 | 0.003 |
| `GASTRO & OMNIPAQUE 300` | — | — | — | 3 | 0.333 | 0.910 |
| `H20 & OMNI` | — | — | — | 3 | 0.667 | 0.333 |
| `H20 & OMNI300` | — | — | — | 1 | 1.000 | 0.051 |
| `H20 & OMNIPAQUE` | — | — | — | 3 | 1.000 | 0.046 |
| `H20 & OPTIRAY` | — | — | — | 2 | 0.000 | 0.860 |
| `H2O & OMNI` | — | — | — | 31 | 0.871 | 0.348 |
| `NON-IONIC` | — | — | — | 2 | 0.500 | 0.723 |
| `OMNI` | — | — | — | 14 | 0.714 | 0.392 |
| `OMNIPAQUE` | — | — | — | 2 | 1.000 | 0.008 |
| `ORAL & OMNI` | — | — | — | 17 | 1.000 | 0.026 |
| `ORAL & OMNI300` | — | — | — | 1 | 1.000 | 0.002 |
| `REDI-CAT & OMNIPAQUE` | — | — | — | 5 | 1.000 | 0.033 |
| `WATER  & OMNI` | — | — | — | 3 | 1.000 | 0.424 |
| `WATER & NON IONIC` | — | — | — | 2 | 0.000 | 0.956 |
| `WATER & OMNI` | — | — | — | 193 | 0.741 | 0.363 |
| `WATER & OMNIPAQUE` | — | — | — | 5 | 1.000 | 0.014 |
| `oral  & omni` | — | — | — | 2 | 1.000 | 0.515 |

### CTDIvol `(189345)` — Intensity/exposure · ✗ not mentioned

- **HCC**  · fill 0.7% · status fully_missing · FDR —
- **CRLM** · fill 0.0% · status fully_missing · FDR —

| bin | HCC n | HCC sens | HCC mean p | CRLM n | CRLM spec | CRLM mean p |
|---|---|---|---|---|---|---|
| `17.6` | 1 | 1.000 | 0.998 | — | — | — |
| `(missing)` | 145 | 0.821 | 0.821 | 464 | 0.793 | 0.320 |

### ContrastFlowRate `(181046)` — Intensity/exposure · ✓ mentioned

- **HCC**  · fill 0.0% · status fully_missing · FDR —
- **CRLM** · fill 0.0% · status fully_missing · FDR —

| bin | HCC n | HCC sens | HCC mean p | CRLM n | CRLM spec | CRLM mean p |
|---|---|---|---|---|---|---|
| `(missing)` | 146 | 0.822 | 0.822 | 464 | 0.793 | 0.320 |

### Exposure_uAs `(181153)` — Intensity/exposure · ✗ not mentioned

- **HCC**  · fill 98.6% · status **★ FDR-sig** · FDR 0.0021 — derivation-augmented (ExposureTime×XRayTubeCurrent / IPP z-diff)
- **CRLM** · fill 44.4% · status testable · FDR 0.4668

| bin | HCC n | HCC sens | HCC mean p | CRLM n | CRLM spec | CRLM mean p |
|---|---|---|---|---|---|---|
| `[0, 150000)` | 2 | 1.000 | 0.913 | 1 | 1.000 | 0.014 |
| `[150000, 250000)` | 14 | 0.786 | 0.751 | 29 | 0.828 | 0.287 |
| `[250000, 350000)` | 69 | 0.768 | 0.762 | 29 | 0.724 | 0.420 |
| `[350000, 450000)` | 42 | 0.905 | 0.913 | 83 | 0.843 | 0.236 |
| `[450000, 1000000)` | 17 | 0.824 | 0.872 | 64 | 0.766 | 0.324 |
| `(missing)` | 2 | 1.000 | 0.983 | 258 | 0.787 | 0.340 |

### RevolutionTime `(189305)` — Intensity/exposure · ✗ not mentioned

- **HCC**  · fill 52.7% · status invariant · FDR —
- **CRLM** · fill 5.4% · status insufficient · FDR —

| bin | HCC n | HCC sens | HCC mean p | CRLM n | CRLM spec | CRLM mean p |
|---|---|---|---|---|---|---|
| `0.8` | 77 | 0.844 | 0.851 | 2 | 1.000 | 0.005 |
| `(missing)` | 69 | 0.797 | 0.791 | 439 | 0.795 | 0.318 |
| `0.7` | — | — | — | 23 | 0.739 | 0.392 |

### ReconstructionAlgorithm `(189315)` — Intensity/exposure · ✗ not mentioned

- **HCC**  · fill 0.0% · status fully_missing · FDR —
- **CRLM** · fill 0.0% · status fully_missing · FDR —

| bin | HCC n | HCC sens | HCC mean p | CRLM n | CRLM spec | CRLM mean p |
|---|---|---|---|---|---|---|
| `(missing)` | 146 | 0.822 | 0.822 | 464 | 0.793 | 0.320 |

### ContrastBolusAgentPhase `(189344)` — Intensity/exposure · ✓ mentioned

- **HCC**  · fill 0.0% · status fully_missing · FDR —
- **CRLM** · fill 0.0% · status fully_missing · FDR —

| bin | HCC n | HCC sens | HCC mean p | CRLM n | CRLM spec | CRLM mean p |
|---|---|---|---|---|---|---|
| `(missing)` | 146 | 0.822 | 0.822 | 464 | 0.793 | 0.320 |

### ▸ Geometry

### SliceThickness `(180050)` — Geometry · ✓ mentioned

- **HCC**  · fill 100.0% · status insufficient · FDR —
- **CRLM** · fill 100.0% · status **★ FDR-sig** · FDR 0.0056

| bin | HCC n | HCC sens | HCC mean p | CRLM n | CRLM spec | CRLM mean p |
|---|---|---|---|---|---|---|
| `2.5` | 143 | 0.818 | 0.820 | 210 | 0.762 | 0.366 |
| `3.0` | 1 | 1.000 | 0.998 | — | — | — |
| `5.0` | 2 | 1.000 | 0.929 | 197 | 0.802 | 0.308 |
| `1.25` | — | — | — | 1 | 1.000 | 0.602 |
| `7.5` | — | — | — | 56 | 0.875 | 0.186 |

### Rows `(280010)` — Geometry · ✗ not mentioned

- **HCC**  · fill 100.0% · status invariant · FDR —
- **CRLM** · fill 100.0% · status invariant · FDR —

| bin | HCC n | HCC sens | HCC mean p | CRLM n | CRLM spec | CRLM mean p |
|---|---|---|---|---|---|---|
| `512` | 146 | 0.822 | 0.822 | 464 | 0.793 | 0.320 |

### Columns `(280011)` — Geometry · ✗ not mentioned

- **HCC**  · fill 100.0% · status invariant · FDR —
- **CRLM** · fill 100.0% · status invariant · FDR —

| bin | HCC n | HCC sens | HCC mean p | CRLM n | CRLM spec | CRLM mean p |
|---|---|---|---|---|---|---|
| `512` | 146 | 0.822 | 0.822 | 464 | 0.793 | 0.320 |

### PixelSpacing `(280030)` — Geometry · ✗ not mentioned

- **HCC**  · fill 100.0% · status testable · FDR 0.1974
- **CRLM** · fill 100.0% · status testable · FDR 0.4668

| bin | HCC n | HCC sens | HCC mean p | CRLM n | CRLM spec | CRLM mean p |
|---|---|---|---|---|---|---|
| `[0.50, 0.65)` | 4 | 1.000 | 0.973 | 22 | 0.727 | 0.382 |
| `[0.65, 0.75)` | 62 | 0.742 | 0.752 | 135 | 0.785 | 0.293 |
| `[0.75, 0.85)` | 54 | 0.870 | 0.864 | 197 | 0.807 | 0.312 |
| `[0.85, 1.00)` | 26 | 0.885 | 0.881 | 110 | 0.791 | 0.357 |

### ReconstructionDiameter `(181100)` — Geometry · ✗ not mentioned

- **HCC**  · fill 100.0% · status testable · FDR 0.2182
- **CRLM** · fill 100.0% · status **★ FDR-sig** · FDR 0.0150

| bin | HCC n | HCC sens | HCC mean p | CRLM n | CRLM spec | CRLM mean p |
|---|---|---|---|---|---|---|
| `[0, 320)` | 1 | 1.000 | 0.984 | 2 | 0.000 | 0.873 |
| `[320, 360)` | 10 | 1.000 | 0.966 | 65 | 0.892 | 0.191 |
| `[360, 400)` | 55 | 0.709 | 0.725 | 154 | 0.734 | 0.363 |
| `[400, 500)` | 78 | 0.872 | 0.867 | 228 | 0.820 | 0.311 |
| `[500+)` | 2 | 1.000 | 0.988 | 15 | 0.667 | 0.507 |

### WindowCenter `(281050)` — Geometry · ✗ not mentioned

- **HCC**  · fill 100.0% · status **★ FDR-sig** · FDR 0.0488
- **CRLM** · fill 100.0% · status **★ FDR-sig** · FDR 0.0233

| bin | HCC n | HCC sens | HCC mean p | CRLM n | CRLM spec | CRLM mean p |
|---|---|---|---|---|---|---|
| `40` | 105 | 0.867 | 0.857 | 438 | 0.790 | 0.329 |
| `55` | 2 | 1.000 | 0.909 | — | — | — |
| `70` | 39 | 0.692 | 0.726 | — | — | — |
| `[40, -500]` | — | — | — | 1 | 1.000 | 0.014 |
| `[40, 40.000000, 100.000000]` | — | — | — | 2 | 0.500 | 0.485 |
| `25` | — | — | — | 1 | 0.000 | 0.971 |
| `30` | — | — | — | 5 | 1.000 | 0.040 |
| `35` | — | — | — | 5 | 0.600 | 0.379 |
| `50` | — | — | — | 12 | 1.000 | 0.050 |

### WindowWidth `(281051)` — Geometry · ✗ not mentioned

- **HCC**  · fill 100.0% · status insufficient · FDR —
- **CRLM** · fill 100.0% · status testable · FDR 0.1114

| bin | HCC n | HCC sens | HCC mean p | CRLM n | CRLM spec | CRLM mean p |
|---|---|---|---|---|---|---|
| `400` | 144 | 0.819 | 0.821 | 431 | 0.794 | 0.326 |
| `500` | 2 | 1.000 | 0.909 | — | — | — |
| `[400, 1500]` | — | — | — | 1 | 1.000 | 0.014 |
| `[440, 440.000000, 150.000000]` | — | — | — | 2 | 0.500 | 0.485 |
| `350` | — | — | — | 3 | 0.667 | 0.346 |
| `370` | — | — | — | 1 | 0.000 | 0.971 |
| `450` | — | — | — | 26 | 0.846 | 0.199 |

### TotalCollimationWidth `(189307)` — Geometry · ✗ not mentioned

- **HCC**  · fill 52.7% · status **★ FDR-sig** · FDR 0.0140
- **CRLM** · fill 5.4% · status testable · FDR 0.9535

| bin | HCC n | HCC sens | HCC mean p | CRLM n | CRLM spec | CRLM mean p |
|---|---|---|---|---|---|---|
| `[20, 40)` | 40 | 0.900 | 0.901 | 17 | 0.824 | 0.340 |
| `[40, 80)` | 37 | 0.784 | 0.797 | 8 | 0.625 | 0.405 |
| `(missing)` | 69 | 0.797 | 0.791 | 439 | 0.795 | 0.318 |

### SpiralPitchFactor `(189311)` — Geometry · ✓ mentioned

- **HCC**  · fill 52.7% · status **★ FDR-sig** · FDR 0.0140
- **CRLM** · fill 5.4% · status testable · FDR 0.9535

| bin | HCC n | HCC sens | HCC mean p | CRLM n | CRLM spec | CRLM mean p |
|---|---|---|---|---|---|---|
| `0.938` | 40 | 0.900 | 0.901 | 17 | 0.824 | 0.340 |
| `0.984` | 37 | 0.784 | 0.797 | 8 | 0.625 | 0.405 |
| `(missing)` | 69 | 0.797 | 0.791 | 439 | 0.795 | 0.318 |

### TableSpeed `(189309)` — Geometry · ✗ not mentioned

- **HCC**  · fill 52.7% · status **★ FDR-sig** · FDR 0.0140
- **CRLM** · fill 5.4% · status testable · FDR 0.9535

| bin | HCC n | HCC sens | HCC mean p | CRLM n | CRLM spec | CRLM mean p |
|---|---|---|---|---|---|---|
| `[20, 40)` | 40 | 0.900 | 0.901 | 17 | 0.824 | 0.340 |
| `[40, 60)` | 37 | 0.784 | 0.797 | 8 | 0.625 | 0.405 |
| `(missing)` | 69 | 0.797 | 0.791 | 439 | 0.795 | 0.318 |

### SpacingBetweenSlices `(180088)` — Geometry · ✗ not mentioned

- **HCC**  · fill 100.0% · status insufficient · FDR — — derivation-augmented (ExposureTime×XRayTubeCurrent / IPP z-diff)
- **CRLM** · fill 100.0% · status **★ FDR-sig** · FDR 0.0274

| bin | HCC n | HCC sens | HCC mean p | CRLM n | CRLM spec | CRLM mean p |
|---|---|---|---|---|---|---|
| `[0.0, 1.5)` | 4 | 1.000 | 0.984 | 5 | 1.000 | 0.405 |
| `[1.5, 3.0)` | 139 | 0.813 | 0.815 | 233 | 0.747 | 0.371 |
| `[3.0, 5.0)` | 1 | 1.000 | 0.998 | — | — | — |
| `[5.0, 8.0)` | 2 | 1.000 | 0.929 | 226 | 0.836 | 0.266 |

### ReformattingThickness `(720512)` — Geometry · ✓ mentioned

- **HCC**  · fill 0.0% · status fully_missing · FDR —
- **CRLM** · fill 0.0% · status fully_missing · FDR —

| bin | HCC n | HCC sens | HCC mean p | CRLM n | CRLM spec | CRLM mean p |
|---|---|---|---|---|---|---|
| `(missing)` | 146 | 0.822 | 0.822 | 464 | 0.793 | 0.320 |

### ▸ Device

### ManufacturerModelName `(081090)` — Device · ✓ mentioned

- **HCC**  · fill 1.4% · status fully_missing · FDR —
- **CRLM** · fill 100.0% · status testable · FDR 0.9535

| bin | HCC n | HCC sens | HCC mean p | CRLM n | CRLM spec | CRLM mean p |
|---|---|---|---|---|---|---|
| `(missing)` | 144 | 0.826 | 0.827 | — | — | — |
| `LightSpeed Plus` | 2 | 0.500 | 0.494 | 11 | 0.909 | 0.212 |
| `Emotion Duo` | — | — | — | 1 | 1.000 | 0.014 |
| `LightSpeed QX/i` | — | — | — | 50 | 0.780 | 0.324 |
| `LightSpeed Ultra` | — | — | — | 21 | 0.810 | 0.236 |
| `LightSpeed VCT` | — | — | — | 38 | 0.789 | 0.330 |
| `LightSpeed16` | — | — | — | 341 | 0.792 | 0.327 |
| `Philips CT Secura` | — | — | — | 2 | 0.500 | 0.485 |

### Manufacturer `(080070)` — Device · ✓ mentioned

- **HCC**  · fill 1.4% · status fully_missing · FDR —
- **CRLM** · fill 100.0% · status insufficient · FDR —

| bin | HCC n | HCC sens | HCC mean p | CRLM n | CRLM spec | CRLM mean p |
|---|---|---|---|---|---|---|
| `(missing)` | 144 | 0.826 | 0.827 | — | — | — |
| `GE MEDICAL SYSTEMS` | 2 | 0.500 | 0.494 | 461 | 0.794 | 0.320 |
| `Philips Medical Systems` | — | — | — | 2 | 0.500 | 0.485 |
| `SIEMENS` | — | — | — | 1 | 1.000 | 0.014 |

### ▸ Anatomy/position

### BodyPartExamined `(180015)` — Anatomy/position · ✓ mentioned

- **HCC**  · fill 100.0% · status invariant · FDR —
- **CRLM** · fill 0.6% · status fully_missing · FDR —

| bin | HCC n | HCC sens | HCC mean p | CRLM n | CRLM spec | CRLM mean p |
|---|---|---|---|---|---|---|
| `LIVER` | 146 | 0.822 | 0.822 | — | — | — |
| `(missing)` | — | — | — | 461 | 0.794 | 0.320 |
| `ABDOMEN` | — | — | — | 3 | 0.667 | 0.328 |

### PatientPosition `(185100)` — Anatomy/position · ✗ not mentioned

- **HCC**  · fill 100.0% · status invariant · FDR —
- **CRLM** · fill 100.0% · status insufficient · FDR —

| bin | HCC n | HCC sens | HCC mean p | CRLM n | CRLM spec | CRLM mean p |
|---|---|---|---|---|---|---|
| `FFS` | 146 | 0.822 | 0.822 | 463 | 0.793 | 0.321 |
| `FFP` | — | — | — | 1 | 1.000 | 0.002 |
