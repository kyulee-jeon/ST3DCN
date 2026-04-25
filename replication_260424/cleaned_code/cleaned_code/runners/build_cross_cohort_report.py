"""Cross-cohort essential-vs-mentioned report.

Usage:
    python build_cross_cohort_report.py \\
        --hcc-dir out_hcc \\
        --crlm-dir out_crlm \\
        --out     RESULTS_ESSENTIAL_vs_MENTIONED.md
"""
import argparse
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--hcc-dir')
    parser.add_argument('--crlm-dir')
    parser.add_argument('--out', required=True)
    args = parser.parse_args()

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from analysis.combined_report import build_report

    configs = {}
    if args.hcc_dir:
        configs['hcc_tace_seg'] = {
            'cls_csv': f'{args.hcc_dir}/per_tag_classification.csv',
            'bin_csv': f'{args.hcc_dir}/per_tag_binned.csv',
            'label':   'HCC-TACE-Seg (HCC+)',
            'metric':  'sens',
        }
    if args.crlm_dir:
        configs['crlm'] = {
            'cls_csv': f'{args.crlm_dir}/per_tag_classification.csv',
            'bin_csv': f'{args.crlm_dir}/per_tag_binned.csv',
            'label':   'CRLM (HCC−)',
            'metric':  'spec',
        }
    build_report(configs, args.out,
                   thesis_title='Essential vs Mentioned DICOM Tags — ST3DCN replication (cross-cohort)')


if __name__ == '__main__':
    main()
