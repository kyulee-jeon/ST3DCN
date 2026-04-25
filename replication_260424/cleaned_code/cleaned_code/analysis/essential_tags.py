"""The protocol2dicom essential-26 CT-tag set, with mentioned-flag and
binning hints used by the per-tag heterogeneity analyzer.

- **Essential**: 26 CT tags deemed necessary to capture acquisition
  parameter variability for CT-based AI replication.
- **Mentioned**: subset that ST3DCN paper (Yu et al., JHEP Rep 2025)
  lists in Methods / Table S1 as acquisition parameters.
"""
# (column, dicom_hex, class_code, class_name, binner_key, mentioned_in_st3dcn)
ESSENTIAL_26 = [
    ('KVP',                     '180060', 'i', 'Intensity/exposure', 'cat',              True),
    ('XRayTubeCurrent',         '181151', 'i', 'Intensity/exposure', 'num_tube_current', True),
    ('ConvolutionKernel',       '181210', 'i', 'Intensity/exposure', 'kernel_family',    False),
    ('ContrastBolusAgent',      '180010', 'i', 'Intensity/exposure', 'cat',              True),
    ('CTDIvol',                 '189345', 'i', 'Intensity/exposure', 'cat',              False),
    ('ContrastFlowRate',        '181046', 'i', 'Intensity/exposure', 'cat',              True),
    ('Exposure_uAs',            '181153', 'i', 'Intensity/exposure', 'num_exposure_uas', False),
    ('RevolutionTime',          '189305', 'i', 'Intensity/exposure', 'num_revtime',      False),
    ('ReconstructionAlgorithm', '189315', 'i', 'Intensity/exposure', 'cat',              False),
    ('ContrastBolusAgentPhase', '189344', 'i', 'Intensity/exposure', 'cat',              True),
    ('SliceThickness',          '180050', 'g', 'Geometry',           'cat',              True),
    ('Rows',                    '280010', 'g', 'Geometry',           'cat',              False),
    ('Columns',                 '280011', 'g', 'Geometry',           'cat',              False),
    ('PixelSpacing',            '280030', 'g', 'Geometry',           'num_pixspacing',   False),
    ('ReconstructionDiameter',  '181100', 'g', 'Geometry',           'num_recondiam',    False),
    ('WindowCenter',            '281050', 'g', 'Geometry',           'cat_multival',     False),
    ('WindowWidth',             '281051', 'g', 'Geometry',           'cat_multival',     False),
    ('TotalCollimationWidth',   '189307', 'g', 'Geometry',           'num_tcw',          False),
    ('SpiralPitchFactor',       '189311', 'g', 'Geometry',           'num_pitch',        True),
    ('TableSpeed',              '189309', 'g', 'Geometry',           'num_tablespeed',   False),
    ('SpacingBetweenSlices',    '180088', 'g', 'Geometry',           'num_spacing',      False),
    ('ReformattingThickness',   '720512', 'g', 'Geometry',           'cat',              True),
    ('ManufacturerModelName',   '081090', 'd', 'Device',             'cat',              True),
    ('Manufacturer',            '080070', 'd', 'Device',             'cat',              True),
    ('BodyPartExamined',        '180015', 'a', 'Anatomy/position',   'cat',              True),
    ('PatientPosition',         '185100', 'a', 'Anatomy/position',   'cat',              False),
]


# Three sequence-type ST3DCN-mentioned "tags" that aren't value-typed DICOM fields
# → excluded from per-tag testing (documented, non-testable):
MENTIONED_SEQUENCES_EXCLUDED = [
    ('ContrastBolusAgentSequence',     '180012'),
    ('ContrastBolusInjectionDelay',    '1811B7'),
    ('CTAcquisitionTypeSequence',      '189301'),
]


def as_dataframe():
    import pandas as pd
    return pd.DataFrame(
        ESSENTIAL_26,
        columns=['tag','hex','class_code','class_name','binner','mentioned']
    )
