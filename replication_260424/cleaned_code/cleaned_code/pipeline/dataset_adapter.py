"""Dataset-specific adapter — tells the generic pipeline how each cohort
encodes its SEG segments and how to pick the right CT.

Add a new dataset by subclassing `DatasetAdapter` and registering it in
`ADAPTERS`. See `CRLMAdapter` and `HCCTACESegAdapter` as templates.
"""
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Iterable, Optional, List, Tuple

import pandas as pd


@dataclass
class ProcessedCase:
    """Output of adapter.process_patient (per-lesion) or process_patient_aggregate."""
    patient_id: str
    ct_uid: str
    seg_uid: str
    study_uid: str
    ct_shape: Tuple[int, int, int]
    pix_spacing_mm: Tuple[float, float]
    slice_thick_mm: float
    phase_sanity: dict = field(default_factory=dict)  # e.g., {'portal_vein_median_hu': 148}
    lesions: List[dict] = field(default_factory=list)   # one dict per lesion with 'crop', 'bbox', ...


class DatasetAdapter(ABC):
    """One subclass per cohort. Encapsulates cohort-specific rules:
    - Which SEG segment numbers are tumors
    - Which SEG segment is used for phase sanity (portal vein or aorta)
    - How to identify the right CT series when multiple exist
    - How to identify pre-treatment study (for longitudinal cohorts)
    """

    # Required subclass attributes
    dicom_root: str
    tumor_segments: Iterable[int]
    portal_vein_segment: Optional[int] = None  # for phase sanity check
    aorta_segment: Optional[int] = None          # for HCC-TACE-Seg-style portal-acq picking
    min_portal_vein_hu: float = 120.0

    @abstractmethod
    def list_patients(self) -> List[str]: ...

    @abstractmethod
    def resolve_series(self, patient_id: str) -> Tuple[str, str, str]:
        """Return (ct_series_uid, seg_series_uid, study_uid). Pre-treatment-only.
        For multi-CT cohorts, pick the mask-aligned CT."""

    def acq_filter(self, patient_id: str, frame_meta, ct_index) -> Optional[int]:
        """Return AcquisitionNumber to keep (e.g., portal) or None for all."""
        return None


# ======================================================================
# CRLM (Colorectal-Liver-Metastases) — MSKCC, preop portal venous only
# ======================================================================
class CRLMAdapter(DatasetAdapter):
    tumor_segments = (5, 6, 7, 8, 9)     # Tumor_1..Tumor_5 (PropType=Mass)
    portal_vein_segment = 4
    aorta_segment = None
    min_portal_vein_hu = 120.0

    def __init__(self, dicom_root='/home/ubuntu/non-hcc_data/Colorectal-Liver-Metastases',
                 manifest_json=None):
        self.dicom_root = dicom_root
        self.manifest_json = manifest_json or f'{dicom_root}/series_manifest.json'
        import json
        with open(self.manifest_json) as f:
            rows = json.load(f)
        self.manifest = pd.DataFrame(rows)

    def list_patients(self):
        return sorted(self.manifest['PatientID'].unique().tolist())

    def resolve_series(self, patient_id: str):
        p = self.manifest[self.manifest['PatientID'] == patient_id]
        ct = p[p['Modality'] == 'CT']
        seg = p[p['Modality'] == 'SEG']
        assert len(ct) == 1 and len(seg) == 1, \
            f'{patient_id}: expected 1 CT + 1 SEG, got {len(ct)} + {len(seg)}'
        return (ct.iloc[0]['SeriesInstanceUID'], seg.iloc[0]['SeriesInstanceUID'],
                ct.iloc[0]['StudyInstanceUID'])

    def acq_filter(self, patient_id, frame_meta, ct_index):
        return None  # CRLM: each CT series has a single AcqNum


# ======================================================================
# HCC-TACE-Seg — multi-phase, multi-study, SEG with aorta segment
# ======================================================================
class HCCTACESegAdapter(DatasetAdapter):
    tumor_segments = (2,)                # Mass
    portal_vein_segment = None           # no explicit portal vein segment
    aorta_segment = 4
    min_portal_vein_hu = None

    def __init__(self,
                 dicom_root='/home/ubuntu/hcc_data/dicom',
                 manifest_csv='/home/ubuntu/hcc_data/series_manifest.csv',
                 nbia_digest='/home/ubuntu/hcc_data/HCC-TACE-Seg_v1_202201-nbia-digest.xlsx',
                 clinical='/home/ubuntu/hcc_data/HCC-TACE-Seg_clinical_data-V2.xlsx'):
        self.dicom_root = dicom_root
        self.manifest = pd.read_csv(manifest_csv)
        self.digest = pd.read_excel(nbia_digest)
        self.clinical = pd.read_excel(clinical, sheet_name='data table')

    def list_patients(self):
        return sorted(self.clinical['TCIA_ID'].unique().tolist())

    def _baseline_study(self, patient_id):
        """Earliest study that contains a SEG (pre-TACE baseline)."""
        d = self.digest[self.digest['Patient ID'] == patient_id].copy()
        assert not d.empty, f'Patient {patient_id} not in NBIA digest'
        off = pd.to_numeric(
            d['Longitudinal Temporal Offset From Event']
                .astype(str).str.lstrip("'"),
            errors='coerce')
        d['_offset'] = off
        studies = (d.groupby('Study Instance UID')['_offset'].min()
                    .reset_index().sort_values('_offset'))
        seg_studies = set(self.manifest[(self.manifest['PatientID'] == patient_id) &
                                         (self.manifest['Modality'] == 'SEG')]['StudyInstanceUID'])
        candidates = studies[studies['Study Instance UID'].isin(seg_studies)]
        if not candidates.empty:
            return candidates.iloc[0]['Study Instance UID']
        return studies.iloc[0]['Study Instance UID']

    def resolve_series(self, patient_id: str):
        study = self._baseline_study(patient_id)
        s = self.manifest[(self.manifest['PatientID'] == patient_id) &
                          (self.manifest['StudyInstanceUID'] == study)]
        segs = s[s['Modality'] == 'SEG']['SeriesInstanceUID'].tolist()
        cts = s[s['Modality'] == 'CT'][['SeriesInstanceUID', 'SeriesDescription']].values.tolist()
        assert segs, f'{patient_id}: no SEG in baseline study {study}'
        # Caller must pick mask-aligned CT via preprocessing.pick_mask_aligned_ct
        # For this adapter, we return (None, seg_uid, study) and let caller iterate.
        return (cts, segs[0], study)  # `cts` is list of (uid, desc) for picker

    def acq_filter(self, patient_id, frame_meta, ct_index):
        """HCC-TACE-Seg: one CT series can contain arterial+portal+delayed as
        separate AcquisitionNumbers. Use aorta (seg 4) max HU — portal is lower."""
        from .preprocessing import identify_portal_acq_by_aorta
        portal_acq, _ = identify_portal_acq_by_aorta(frame_meta, ct_index,
                                                       aorta_seg=self.aorta_segment)
        return portal_acq


ADAPTERS = {
    'crlm':          CRLMAdapter,
    'hcc_tace_seg':  HCCTACESegAdapter,
}


def get_adapter(name: str, **kwargs) -> DatasetAdapter:
    cls = ADAPTERS.get(name)
    assert cls is not None, f'Unknown dataset: {name}. Known: {list(ADAPTERS)}'
    return cls(**kwargs)
