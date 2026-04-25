# Per-essential-tag SPECIFICITY + probability heterogeneity (CRLM, ST3DCN)

**Cohort**: 197 CRLM patients, 100% HCC-negative (specificity-only, no HCC ground truth).

**Threshold**: 0.8  ·  **Binary pred**: `prob ≥ 0.8` → pred HCC (FP here).

**Baseline pu**: n=197, TN=60, FP=137, specificity = **0.3046**, Wilson 95% CI = (0.2445, 0.3721).

Primary analysis: Kruskal-Wallis on probabilities across bins, BH-FDR across testable tags (same framework as HCC-TACE-Seg RESULTS_v5).

## Essential-26 classification

- fully_missing: 8
- invariant:     2
- insufficient:  5
- testable:      11
  - significant (FDR<0.05): 1
  - not significant:                10

## Significant tags (FDR<0.05)

| Tag | Hex | Class | Mentioned | KW H | raw p | FDR p | mean range | spec range |
|---|---|---|---|---|---|---|---|---|
| **SliceThickness** | 180050 | Geometry | ✓ | 12.72 | 0.0017 | **0.0190** | 0.569–0.792 | 0.253–0.481 |

## Tag availability overview

| Tag | Hex | Class | Fill% | # bins | Spec range | Mean prob range | Status |
|---|---|---|---|---|---|---|---|
| KVP | 180060 | I | 100.0 | 3 | 0.303–0.303 | 0.755–0.755 | insufficient |
| XRayTubeCurrent | 181151 | I | 48.2 | 5 | 0.000–0.455 | 0.644–0.940 | testable |
| ConvolutionKernel | 181210 | I | 48.2 | 3 | 0.355–0.355 | 0.688–0.688 | insufficient |
| ContrastBolusAgent | 180010 | I | 99.0 | 27 | 0.200–0.400 | 0.603–0.935 | testable |
| CTDIvol | 189345 | I | 0.0 | 0 | — | — | fully_missing |
| ContrastFlowRate | 181046 | I | 0.0 | 0 | — | — | fully_missing |
| Exposure_uAs | 181153 | I | 0.0 | 0 | — | — | fully_missing |
| RevolutionTime | 189305 | I | 6.1 | 2 | 0.182–0.182 | 0.778–0.778 | insufficient |
| ReconstructionAlgorithm | 189315 | I | 0.0 | 0 | — | — | fully_missing |
| ContrastBolusAgentPhase | 189344 | I | 0.0 | 0 | — | — | fully_missing |
| SliceThickness | 180050 | G | 100.0 | 4 | 0.253–0.481 | 0.569–0.792 | testable |
| Rows | 280010 | G | 100.0 | 1 | 0.305–0.305 | 0.752–0.752 | invariant |
| Columns | 280011 | G | 100.0 | 1 | 0.305–0.305 | 0.752–0.752 | invariant |
| PixelSpacing | 280030 | G | 100.0 | 4 | 0.205–0.385 | 0.699–0.804 | testable |
| ReconstructionDiameter | 181100 | G | 100.0 | 5 | 0.275–0.423 | 0.642–0.776 | testable |
| WindowCenter | 281050 | G | 100.0 | 7 | 0.297–0.667 | 0.663–0.754 | testable |
| WindowWidth | 281051 | G | 100.0 | 6 | 0.297–0.500 | 0.690–0.754 | testable |
| TotalCollimationWidth | 189307 | G | 6.1 | 2 | 0.143–0.200 | 0.751–0.825 | testable |
| SpiralPitchFactor | 189311 | G | 6.1 | 2 | 0.143–0.200 | 0.751–0.825 | testable |
| TableSpeed | 189309 | G | 6.1 | 2 | 0.143–0.200 | 0.751–0.825 | testable |
| SpacingBetweenSlices | 180088 | G | 0.0 | 0 | — | — | fully_missing |
| ReformattingThickness | 720512 | G | 0.0 | 0 | — | — | fully_missing |
| ManufacturerModelName | 081090 | D | 100.0 | 7 | 0.222–0.429 | 0.656–0.793 | testable |
| Manufacturer | 080070 | D | 100.0 | 3 | 0.303–0.303 | 0.755–0.755 | insufficient |
| BodyPartExamined | 180015 | A | 1.0 | 1 | — | — | fully_missing |
| PatientPosition | 185100 | A | 100.0 | 2 | 0.301–0.301 | 0.756–0.756 | insufficient |

*Class*: **i** = intensity/exposure, **g** = geometry, **d** = device, **a** = anatomy/position.

## Intensity/exposure

### KVP (`180060`)

| Value / bin | n | TN | FP | Spec | Wilson 95% CI | mean p |
|---|---|---|---|---|---|---|
| `120.0` | 195 | 59 | 136 | 0.3026 | (0.2424, 0.3703) | 0.755 |
| `130.0` | 1 | 1 | 0 | 1.0000 | (0.2065, 1.0000) | 0.014 |
| `140.0` | 1 | 0 | 1 | 0.0000 | (0.0000, 0.7935) | 0.978 |

### XRayTubeCurrent (`181151`)

| Value / bin | n | TN | FP | Spec | Wilson 95% CI | mean p |
|---|---|---|---|---|---|---|
| `[0, 200)` | 7 | 2 | 5 | 0.2857 | (0.0822, 0.6411) | 0.672 |
| `[200, 300)` | 26 | 11 | 15 | 0.4231 | (0.2554, 0.6105) | 0.651 |
| `[300, 400)` | 33 | 11 | 22 | 0.3333 | (0.1975, 0.5039) | 0.684 |
| `[400, 500)` | 22 | 10 | 12 | 0.4545 | (0.2692, 0.6534) | 0.644 |
| `[500, 10000)` | 7 | 0 | 7 | 0.0000 | (0.0000, 0.3543) | 0.940 |
| `(missing)` | 102 | 26 | 76 | 0.2549 | (0.1803, 0.3473) | 0.816 |

### ConvolutionKernel (`181210`)

| Value / bin | n | TN | FP | Spec | Wilson 95% CI | mean p |
|---|---|---|---|---|---|---|
| `(missing)` | 102 | 26 | 76 | 0.2549 | (0.1803, 0.3473) | 0.816 |
| `AA-1` | 1 | 0 | 1 | 0.0000 | (0.0000, 0.7935) | 0.978 |
| `B40` | 1 | 1 | 0 | 1.0000 | (0.2065, 1.0000) | 0.014 |
| `STANDARD` | 93 | 33 | 60 | 0.3548 | (0.2651, 0.4561) | 0.688 |

### ContrastBolusAgent (`180010`)

| Value / bin | n | TN | FP | Spec | Wilson 95% CI | mean p |
|---|---|---|---|---|---|---|
| `&` | 1 | 0 | 1 | 0.0000 | (0.0000, 0.7935) | 0.913 |
| `(missing)` | 2 | 0 | 2 | 0.0000 | (0.0000, 0.6576) | 0.990 |
| `APPLIED` | 1 | 1 | 0 | 1.0000 | (0.2065, 1.0000) | 0.014 |
| `BARIUM  & OMNI` | 1 | 1 | 0 | 1.0000 | (0.2065, 1.0000) | 0.003 |
| `BARIUM & NON IONIC` | 2 | 0 | 2 | 0.0000 | (0.0000, 0.6576) | 0.912 |
| `BARIUM & OMNI` | 26 | 10 | 16 | 0.3846 | (0.2243, 0.5747) | 0.603 |
| `CONT` | 1 | 0 | 1 | 0.0000 | (0.0000, 0.7935) | 0.978 |
| `GASTRO & NON IONIC` | 3 | 2 | 1 | 0.6667 | (0.2077, 0.9385) | 0.581 |
| `GASTRO & OMNI` | 43 | 16 | 27 | 0.3721 | (0.2438, 0.5214) | 0.724 |
| `GASTRO & OMNI 300` | 1 | 0 | 1 | 0.0000 | (0.0000, 0.7935) | 0.986 |
| `GASTRO & OMNIPAQUE` | 2 | 1 | 1 | 0.5000 | (0.0945, 0.9055) | 0.431 |
| `GASTRO & OMNIPAQUE 300` | 1 | 0 | 1 | 0.0000 | (0.0000, 0.7935) | 0.996 |
| `H20 & OMNI` | 1 | 0 | 1 | 0.0000 | (0.0000, 0.7935) | 0.921 |
| `H20 & OMNI300` | 1 | 1 | 0 | 1.0000 | (0.2065, 1.0000) | 0.051 |
| `H20 & OMNIPAQUE` | 1 | 0 | 1 | 0.0000 | (0.0000, 0.7935) | 0.985 |
| `H20 & OPTIRAY` | 1 | 0 | 1 | 0.0000 | (0.0000, 0.7935) | 0.994 |
| `H2O & OMNI` | 13 | 5 | 8 | 0.3846 | (0.1771, 0.6448) | 0.743 |
| `NON-IONIC` | 2 | 1 | 1 | 0.5000 | (0.0945, 0.9055) | 0.723 |
| `OMNI` | 5 | 1 | 4 | 0.2000 | (0.0362, 0.6245) | 0.935 |
| `OMNIPAQUE` | 1 | 0 | 1 | 0.0000 | (0.0000, 0.7935) | 0.992 |
| `ORAL & OMNI` | 5 | 2 | 3 | 0.4000 | (0.1176, 0.7693) | 0.858 |
| `ORAL & OMNI300` | 1 | 1 | 0 | 1.0000 | (0.2065, 1.0000) | 0.002 |
| `REDI-CAT & OMNIPAQUE` | 1 | 1 | 0 | 1.0000 | (0.2065, 1.0000) | 0.595 |
| `WATER  & OMNI` | 1 | 0 | 1 | 0.0000 | (0.0000, 0.7935) | 0.989 |
| `WATER & NON IONIC` | 1 | 0 | 1 | 0.0000 | (0.0000, 0.7935) | 0.937 |
| `WATER & OMNI` | 77 | 17 | 60 | 0.2208 | (0.1427, 0.3254) | 0.812 |
| `WATER & OMNIPAQUE` | 1 | 0 | 1 | 0.0000 | (0.0000, 0.7935) | 0.997 |
| `oral  & omni` | 1 | 0 | 1 | 0.0000 | (0.0000, 0.7935) | 0.988 |

- **CTDIvol** (`189345`): fully missing in CRLM (fill=0%).

- **ContrastFlowRate** (`181046`): fully missing in CRLM (fill=0%).

- **Exposure_uAs** (`181153`): fully missing in CRLM (fill=0%).

### RevolutionTime (`189305`)

| Value / bin | n | TN | FP | Spec | Wilson 95% CI | mean p |
|---|---|---|---|---|---|---|
| `0.7` | 11 | 2 | 9 | 0.1818 | (0.0514, 0.4770) | 0.778 |
| `0.8` | 1 | 0 | 1 | 0.0000 | (0.0000, 0.7935) | 0.968 |
| `(missing)` | 185 | 58 | 127 | 0.3135 | (0.2510, 0.3836) | 0.749 |

- **ReconstructionAlgorithm** (`189315`): fully missing in CRLM (fill=0%).

- **ContrastBolusAgentPhase** (`189344`): fully missing in CRLM (fill=0%).

## Geometry

### SliceThickness (`180050`)

| Value / bin | n | TN | FP | Spec | Wilson 95% CI | mean p |
|---|---|---|---|---|---|---|
| `1.25` | 1 | 1 | 0 | 1.0000 | (0.2065, 1.0000) | 0.602 |
| `2.5` | 79 | 20 | 59 | 0.2532 | (0.1703, 0.3589) | 0.792 |
| `5.0` | 90 | 26 | 64 | 0.2889 | (0.2054, 0.3896) | 0.773 |
| `7.5` | 27 | 13 | 14 | 0.4815 | (0.3074, 0.6601) | 0.569 |

- **Rows** (`280010`): single value `512` (n=197, spec=0.305, mean p=0.752) — invariant.

- **Columns** (`280011`): single value `512` (n=197, spec=0.305, mean p=0.752) — invariant.

### PixelSpacing (`280030`)

| Value / bin | n | TN | FP | Spec | Wilson 95% CI | mean p |
|---|---|---|---|---|---|---|
| `[0.50, 0.65)` | 13 | 5 | 8 | 0.3846 | (0.1771, 0.6448) | 0.699 |
| `[0.65, 0.75)` | 57 | 18 | 39 | 0.3158 | (0.2100, 0.4448) | 0.725 |
| `[0.75, 0.85)` | 83 | 28 | 55 | 0.3373 | (0.2448, 0.4442) | 0.752 |
| `[0.85, 1.00)` | 44 | 9 | 35 | 0.2045 | (0.1115, 0.3450) | 0.804 |

### ReconstructionDiameter (`181100`)

| Value / bin | n | TN | FP | Spec | Wilson 95% CI | mean p |
|---|---|---|---|---|---|---|
| `[0, 320)` | 2 | 0 | 2 | 0.0000 | (0.0000, 0.6576) | 0.873 |
| `[320, 360)` | 26 | 11 | 15 | 0.4231 | (0.2554, 0.6105) | 0.642 |
| `[360, 400)` | 69 | 21 | 48 | 0.3043 | (0.2085, 0.4208) | 0.776 |
| `[400, 500)` | 91 | 25 | 66 | 0.2747 | (0.1936, 0.3741) | 0.773 |
| `[500+)` | 9 | 3 | 6 | 0.3333 | (0.1206, 0.6458) | 0.651 |

### WindowCenter (`281050`)

| Value / bin | n | TN | FP | Spec | Wilson 95% CI | mean p |
|---|---|---|---|---|---|---|
| `[40, -500]` | 1 | 1 | 0 | 1.0000 | (0.2065, 1.0000) | 0.014 |
| `[40, 40.000000, 100.000000]` | 1 | 0 | 1 | 0.0000 | (0.0000, 0.7935) | 0.978 |
| `25` | 1 | 0 | 1 | 0.0000 | (0.0000, 0.7935) | 0.971 |
| `30` | 2 | 0 | 2 | 0.0000 | (0.0000, 0.6576) | 0.908 |
| `35` | 1 | 0 | 1 | 0.0000 | (0.0000, 0.7935) | 0.979 |
| `40` | 185 | 55 | 130 | 0.2973 | (0.2361, 0.3667) | 0.754 |
| `50` | 6 | 4 | 2 | 0.6667 | (0.3000, 0.9032) | 0.663 |

### WindowWidth (`281051`)

| Value / bin | n | TN | FP | Spec | Wilson 95% CI | mean p |
|---|---|---|---|---|---|---|
| `[400, 1500]` | 1 | 1 | 0 | 1.0000 | (0.2065, 1.0000) | 0.014 |
| `[440, 440.000000, 150.000000]` | 1 | 0 | 1 | 0.0000 | (0.0000, 0.7935) | 0.978 |
| `350` | 2 | 0 | 2 | 0.0000 | (0.0000, 0.6576) | 0.997 |
| `370` | 1 | 0 | 1 | 0.0000 | (0.0000, 0.7935) | 0.971 |
| `400` | 182 | 54 | 128 | 0.2967 | (0.2351, 0.3667) | 0.754 |
| `450` | 10 | 5 | 5 | 0.5000 | (0.2366, 0.7634) | 0.690 |

### TotalCollimationWidth (`189307`)

| Value / bin | n | TN | FP | Spec | Wilson 95% CI | mean p |
|---|---|---|---|---|---|---|
| `[20, 40)` | 7 | 1 | 6 | 0.1429 | (0.0257, 0.5131) | 0.825 |
| `[40, 80)` | 5 | 1 | 4 | 0.2000 | (0.0362, 0.6245) | 0.751 |
| `(missing)` | 185 | 58 | 127 | 0.3135 | (0.2510, 0.3836) | 0.749 |

### SpiralPitchFactor (`189311`)

| Value / bin | n | TN | FP | Spec | Wilson 95% CI | mean p |
|---|---|---|---|---|---|---|
| `0.938` | 7 | 1 | 6 | 0.1429 | (0.0257, 0.5131) | 0.825 |
| `0.984` | 5 | 1 | 4 | 0.2000 | (0.0362, 0.6245) | 0.751 |
| `(missing)` | 185 | 58 | 127 | 0.3135 | (0.2510, 0.3836) | 0.749 |

### TableSpeed (`189309`)

| Value / bin | n | TN | FP | Spec | Wilson 95% CI | mean p |
|---|---|---|---|---|---|---|
| `[20, 40)` | 7 | 1 | 6 | 0.1429 | (0.0257, 0.5131) | 0.825 |
| `[40, 60)` | 5 | 1 | 4 | 0.2000 | (0.0362, 0.6245) | 0.751 |
| `(missing)` | 185 | 58 | 127 | 0.3135 | (0.2510, 0.3836) | 0.749 |

- **SpacingBetweenSlices** (`180088`): fully missing in CRLM (fill=0%).

- **ReformattingThickness** (`720512`): fully missing in CRLM (fill=0%).

## Device

### ManufacturerModelName (`081090`)

| Value / bin | n | TN | FP | Spec | Wilson 95% CI | mean p |
|---|---|---|---|---|---|---|
| `Emotion Duo` | 1 | 1 | 0 | 1.0000 | (0.2065, 1.0000) | 0.014 |
| `LightSpeed Plus` | 7 | 3 | 4 | 0.4286 | (0.1582, 0.7495) | 0.656 |
| `LightSpeed QX/i` | 22 | 7 | 15 | 0.3182 | (0.1636, 0.5268) | 0.791 |
| `LightSpeed Ultra` | 9 | 2 | 7 | 0.2222 | (0.0632, 0.5474) | 0.743 |
| `LightSpeed VCT` | 15 | 4 | 11 | 0.2667 | (0.1090, 0.5195) | 0.793 |
| `LightSpeed16` | 142 | 43 | 99 | 0.3028 | (0.2333, 0.3828) | 0.751 |
| `Philips CT Secura` | 1 | 0 | 1 | 0.0000 | (0.0000, 0.7935) | 0.978 |

### Manufacturer (`080070`)

| Value / bin | n | TN | FP | Spec | Wilson 95% CI | mean p |
|---|---|---|---|---|---|---|
| `GE MEDICAL SYSTEMS` | 195 | 59 | 136 | 0.3026 | (0.2424, 0.3703) | 0.755 |
| `Philips Medical Systems` | 1 | 0 | 1 | 0.0000 | (0.0000, 0.7935) | 0.978 |
| `SIEMENS` | 1 | 1 | 0 | 1.0000 | (0.2065, 1.0000) | 0.014 |

## Anatomy/position

- **BodyPartExamined** (`180015`): single value `ABDOMEN` (n=2, spec=0.500, mean p=0.496); missing in 195/197.

### PatientPosition (`185100`)

| Value / bin | n | TN | FP | Spec | Wilson 95% CI | mean p |
|---|---|---|---|---|---|---|
| `FFP` | 1 | 1 | 0 | 1.0000 | (0.2065, 1.0000) | 0.002 |
| `FFS` | 196 | 59 | 137 | 0.3010 | (0.2411, 0.3686) | 0.756 |

