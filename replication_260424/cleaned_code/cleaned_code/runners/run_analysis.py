"""Full per-tag analysis for one cohort: extract tags → derive → per-tag KW+FDR.

Usage:
    python run_analysis.py \\
        --cohort crlm \\
        --lesion-csv  out_crlm/lesion_results.csv \\
        --patient-csv out_crlm/patient_results.csv \\
        --out-dir     out_crlm
"""
import argparse
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cohort', required=True, choices=['hcc_tace_seg','crlm'])
    parser.add_argument('--lesion-csv',  required=True)
    parser.add_argument('--patient-csv', required=True)
    parser.add_argument('--out-dir',     required=True)
    args = parser.parse_args()

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from pipeline.dataset_adapter import get_adapter
    from analysis.extract_tags import tag_dataframe
    from analysis.derive_tags import derive_tags
    from analysis.per_tag_analysis import analyze

    adapter = get_adapter(args.cohort)
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    les_tagged     = out/'lesion_results_with_tags.csv'
    pat_tagged     = out/'patient_results_with_tags.csv'
    les_derived    = out/'lesion_results_with_derived.csv'
    pat_derived    = out/'patient_results_with_derived.csv'
    cls_csv        = out/'per_tag_classification.csv'
    bin_csv        = out/'per_tag_binned.csv'
    md_path        = out/'per_tag_report.md'

    print('=== extract DICOM tags ===')
    tag_dataframe(args.lesion_csv,  les_tagged, adapter.dicom_root)
    tag_dataframe(args.patient_csv, pat_tagged, adapter.dicom_root)

    print('\n=== derive standards-based tags ===')
    nbia = getattr(adapter, 'nbia_digest_xlsx', None) or getattr(adapter, 'digest', None)
    # only pass NBIA if it's a path
    nbia_path = None
    if args.cohort == 'hcc_tace_seg':
        nbia_path = '/home/ubuntu/hcc_data/HCC-TACE-Seg_v1_202201-nbia-digest.xlsx'
    derive_tags(str(les_tagged), str(pat_tagged), adapter.dicom_root,
                  str(les_derived), str(pat_derived),
                  nbia_digest=nbia_path)

    print('\n=== per-tag heterogeneity (KW + BH-FDR) ===')
    all_positive = args.cohort == 'hcc_tace_seg'
    analyze(str(les_derived), str(cls_csv), str(bin_csv), str(md_path),
             all_positive_label=all_positive, apply_derivations=True)


if __name__ == '__main__':
    main()
