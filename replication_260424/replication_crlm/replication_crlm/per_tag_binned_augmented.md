# Augmented per-tag specificity + probability heterogeneity (cohort=CRLM, ST3DCN)

Derivation substitutions applied:
- **SpacingBetweenSlices**: substituted from SpacingBetweenSlices_derived (direct fill 0/464)
- **Exposure_uAs**: substituted from Exposure_uAs_derived (direct fill 0/464)

## Essential-26 classification (augmented)

- fully_missing: 6
- invariant: 2
- insufficient: 5
- testable: 13
  - significant (FDR<0.05): 5

## Significant tags (FDR<0.05)

| Tag | Hex | Class | Mentioned | Derived? | KW H | raw p | FDR p | mean range |
|---|---|---|---|---|---|---|---|---|
| **SliceThickness** | 180050 | Geometry | ✓ |  | 14.57 | 0.0007 | **0.0056** | 0.186–0.366 |
| **ContrastBolusAgent** | 180010 | Intensity/exposure | ✓ |  | 26.50 | 0.0009 | **0.0056** | 0.014–0.392 |
| **ReconstructionDiameter** | 181100 | Geometry | ✗ |  | 13.63 | 0.0035 | **0.0150** | 0.191–0.507 |
| **WindowCenter** | 281050 | Geometry | ✗ |  | 12.06 | 0.0072 | **0.0233** | 0.040–0.379 |
| **SpacingBetweenSlices** | 180088 | Geometry | ✗ | ✓ | 9.10 | 0.0105 | **0.0274** | 0.266–0.405 |

## Tag classification after derivation

| Tag | Hex | Mentioned | Derived? | Fill% | status | KW FDR |
|---|---|---|---|---|---|---|
| KVP | 180060 | ✓ |  | 100.0 | insufficient | — |
| XRayTubeCurrent | 181151 | ✓ |  | 44.4 | testable | 0.8182 |
| ConvolutionKernel | 181210 | ✗ |  | 44.4 | insufficient | — |
| ContrastBolusAgent | 180010 | ✓ |  | 99.1 | testable | 0.0056 |
| CTDIvol | 189345 | ✗ |  | 0.0 | fully_missing | — |
| ContrastFlowRate | 181046 | ✓ |  | 0.0 | fully_missing | — |
| Exposure_uAs | 181153 | ✗ | ✓ | 44.4 | testable | 0.4668 |
| RevolutionTime | 189305 | ✗ |  | 5.4 | insufficient | — |
| ReconstructionAlgorithm | 189315 | ✗ |  | 0.0 | fully_missing | — |
| ContrastBolusAgentPhase | 189344 | ✓ |  | 0.0 | fully_missing | — |
| SliceThickness | 180050 | ✓ |  | 100.0 | testable | 0.0056 |
| Rows | 280010 | ✗ |  | 100.0 | invariant | — |
| Columns | 280011 | ✗ |  | 100.0 | invariant | — |
| PixelSpacing | 280030 | ✗ |  | 100.0 | testable | 0.4668 |
| ReconstructionDiameter | 181100 | ✗ |  | 100.0 | testable | 0.0150 |
| WindowCenter | 281050 | ✗ |  | 100.0 | testable | 0.0233 |
| WindowWidth | 281051 | ✗ |  | 100.0 | testable | 0.1114 |
| TotalCollimationWidth | 189307 | ✗ |  | 5.4 | testable | 0.9535 |
| SpiralPitchFactor | 189311 | ✓ |  | 5.4 | testable | 0.9535 |
| TableSpeed | 189309 | ✗ |  | 5.4 | testable | 0.9535 |
| SpacingBetweenSlices | 180088 | ✗ | ✓ | 100.0 | testable | 0.0274 |
| ReformattingThickness | 720512 | ✓ |  | 0.0 | fully_missing | — |
| ManufacturerModelName | 081090 | ✓ |  | 100.0 | testable | 0.9535 |
| Manufacturer | 080070 | ✓ |  | 100.0 | insufficient | — |
| BodyPartExamined | 180015 | ✓ |  | 0.6 | fully_missing | — |
| PatientPosition | 185100 | ✗ |  | 100.0 | insufficient | — |
