# Augmented per-tag sensitivity + probability heterogeneity (cohort=HCC, ST3DCN)

Derivation substitutions applied:
- **SpacingBetweenSlices**: substituted from SpacingBetweenSlices_derived (direct fill 5/146)
- **Exposure_uAs**: substituted from Exposure_uAs_derived (direct fill 0/146)

## Essential-26 classification (augmented)

- fully_missing: 8
- invariant: 5
- insufficient: 4
- testable: 9
  - significant (FDR<0.05): 6

## Significant tags (FDR<0.05)

| Tag | Hex | Class | Mentioned | Derived? | KW H | raw p | FDR p | mean range |
|---|---|---|---|---|---|---|---|---|
| **KVP** | 180060 | Intensity/exposure | ✓ |  | 15.83 | 0.0004 | **0.0021** | 0.757–0.885 |
| **Exposure_uAs** | 181153 | Intensity/exposure | ✗ | ✓ | 17.89 | 0.0005 | **0.0021** | 0.751–0.913 |
| **TotalCollimationWidth** | 189307 | Geometry | ✗ |  | 7.08 | 0.0078 | **0.0140** | 0.797–0.901 |
| **SpiralPitchFactor** | 189311 | Geometry | ✓ |  | 7.08 | 0.0078 | **0.0140** | 0.797–0.901 |
| **TableSpeed** | 189309 | Geometry | ✗ |  | 7.08 | 0.0078 | **0.0140** | 0.797–0.901 |
| **WindowCenter** | 281050 | Geometry | ✗ |  | 4.57 | 0.0325 | **0.0488** | 0.726–0.857 |

## Tag classification after derivation

| Tag | Hex | Mentioned | Derived? | Fill% | status | KW FDR |
|---|---|---|---|---|---|---|
| KVP | 180060 | ✓ |  | 100.0 | testable | 0.0021 |
| XRayTubeCurrent | 181151 | ✓ |  | 99.3 | testable | 0.6483 |
| ConvolutionKernel | 181210 | ✗ |  | 99.3 | insufficient | — |
| ContrastBolusAgent | 180010 | ✓ |  | 1.4 | fully_missing | — |
| CTDIvol | 189345 | ✗ |  | 0.7 | fully_missing | — |
| ContrastFlowRate | 181046 | ✓ |  | 0.0 | fully_missing | — |
| Exposure_uAs | 181153 | ✗ | ✓ | 98.6 | testable | 0.0021 |
| RevolutionTime | 189305 | ✗ |  | 52.7 | invariant | — |
| ReconstructionAlgorithm | 189315 | ✗ |  | 0.0 | fully_missing | — |
| ContrastBolusAgentPhase | 189344 | ✓ |  | 0.0 | fully_missing | — |
| SliceThickness | 180050 | ✓ |  | 100.0 | insufficient | — |
| Rows | 280010 | ✗ |  | 100.0 | invariant | — |
| Columns | 280011 | ✗ |  | 100.0 | invariant | — |
| PixelSpacing | 280030 | ✗ |  | 100.0 | testable | 0.1974 |
| ReconstructionDiameter | 181100 | ✗ |  | 100.0 | testable | 0.2182 |
| WindowCenter | 281050 | ✗ |  | 100.0 | testable | 0.0488 |
| WindowWidth | 281051 | ✗ |  | 100.0 | insufficient | — |
| TotalCollimationWidth | 189307 | ✗ |  | 52.7 | testable | 0.0140 |
| SpiralPitchFactor | 189311 | ✓ |  | 52.7 | testable | 0.0140 |
| TableSpeed | 189309 | ✗ |  | 52.7 | testable | 0.0140 |
| SpacingBetweenSlices | 180088 | ✗ | ✓ | 100.0 | insufficient | — |
| ReformattingThickness | 720512 | ✓ |  | 0.0 | fully_missing | — |
| ManufacturerModelName | 081090 | ✓ |  | 1.4 | fully_missing | — |
| Manufacturer | 080070 | ✓ |  | 1.4 | fully_missing | — |
| BodyPartExamined | 180015 | ✓ |  | 100.0 | invariant | — |
| PatientPosition | 185100 | ✗ |  | 100.0 | invariant | — |
