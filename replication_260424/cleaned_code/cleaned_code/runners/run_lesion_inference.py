"""Per-lesion batch inference runner. Cohort-agnostic via DatasetAdapter.

Usage:
    LD_LIBRARY_PATH=/usr/local/cuda-11.2/lib64:$LD_LIBRARY_PATH \\
    CUDA_VISIBLE_DEVICES=2 \\
    python run_lesion_inference.py --cohort crlm --out-dir ./out_crlm

Outputs:
    <out_dir>/lesion_results.csv
    <out_dir>/patient_results.csv
    <out_dir>/lesion_summary.txt     (sens or spec + Wilson CI)
"""
import argparse
import csv
import os
import sys
import time
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cohort', required=True,
                        choices=['hcc_tace_seg', 'crlm'],
                        help='Registered adapter name')
    parser.add_argument('--st3dcn-repo', default='/home/ubuntu/hcc_workspace/ST3DCN',
                        help='Path to ST3DCN repo clone (provides ST3DCN_Utils)')
    parser.add_argument('--weights',
                        default='/home/ubuntu/hcc_workspace/ST3DCN/weights/ST3DCN_Model.h5')
    parser.add_argument('--out-dir', required=True)
    parser.add_argument('--threshold', type=float, default=0.8)
    parser.add_argument('--min-voxels', type=int, default=50)
    args = parser.parse_args()

    # Ensure cleaned_code package import
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    sys.path.insert(0, args.st3dcn_repo)

    from pipeline.dataset_adapter import get_adapter
    from pipeline.lesion_pipeline import process_patient_lesions
    from pipeline.inference import load_model, predict, binary_pred

    os.makedirs(args.out_dir, exist_ok=True)
    adapter = get_adapter(args.cohort)
    patients = adapter.list_patients()
    print(f'==> {args.cohort}: {len(patients)} patients')

    model = load_model(args.weights)

    lesion_fields = ['patient_id','lesion_idx','segment_number','segment_label',
                      'voxel_count','prob','pred','phase_sanity',
                      'ct_shape','bbox','seg_uid','ct_uid','study_uid',
                      'pix_spacing_mm','slice_thick','n_lesions_patient',
                      'status','error']
    patient_fields = ['patient_id','status','n_lesions','max_prob','min_prob',
                       'mean_prob','any_pred_positive','study_uid','seg_uid','ct_uid',
                       'phase_sanity','elapsed_s','error']
    lesion_rows = []; patient_rows = []
    t_start = time.time()

    for i, pid in enumerate(patients, 1):
        t0 = time.time()
        try:
            rec = process_patient_lesions(adapter, pid, min_voxels=args.min_voxels)
            phase = rec.phase_sanity
            phase_str = ','.join(f'{k}={v:.1f}' if isinstance(v,(int,float)) else f'{k}={v}'
                                  for k,v in phase.items())
            probs = []
            for les in rec.lesions:
                p = predict(model, les['crop'])
                probs.append(p)
                lesion_rows.append({
                    'patient_id': pid, 'lesion_idx': les['lesion_idx'],
                    'segment_number': les['segment_number'],
                    'segment_label': les['segment_label'],
                    'voxel_count': les['voxel_count'],
                    'prob': f'{p:.6f}', 'pred': binary_pred(p, args.threshold),
                    'phase_sanity': phase_str,
                    'ct_shape': 'x'.join(str(s) for s in rec.ct_shape),
                    'bbox': ','.join(str(v) for v in les['bbox']),
                    'seg_uid': rec.seg_uid, 'ct_uid': rec.ct_uid,
                    'study_uid': rec.study_uid,
                    'pix_spacing_mm': ','.join(f'{v:.4f}' for v in rec.pix_spacing_mm),
                    'slice_thick': f'{rec.slice_thick_mm:.2f}',
                    'n_lesions_patient': len(rec.lesions),
                    'status': 'ok', 'error': '',
                })
            patient_rows.append({
                'patient_id': pid, 'status': 'ok',
                'n_lesions': len(rec.lesions),
                'max_prob': f'{max(probs):.6f}',
                'min_prob': f'{min(probs):.6f}',
                'mean_prob': f'{sum(probs)/len(probs):.6f}',
                'any_pred_positive': int(max(probs) >= args.threshold),
                'study_uid': rec.study_uid, 'seg_uid': rec.seg_uid,
                'ct_uid': rec.ct_uid, 'phase_sanity': phase_str,
                'elapsed_s': f'{time.time()-t0:.2f}', 'error': '',
            })
            print(f'[{i:3}/{len(patients)}] {pid}  n_les={len(rec.lesions)}  '
                  f'probs={[f"{p:.3f}" for p in probs]}  '
                  f't={time.time()-t0:.1f}s')
        except Exception as e:
            err = f'{type(e).__name__}: {e}'
            patient_rows.append({k: '' for k in patient_fields})
            patient_rows[-1].update({'patient_id': pid, 'status': 'fail',
                                     'error': err,
                                     'elapsed_s': f'{time.time()-t0:.2f}'})
            print(f'[{i:3}/{len(patients)}] {pid}  FAIL  {err}')

    # Write
    ld = Path(args.out_dir) / 'lesion_results.csv'
    pd_path = Path(args.out_dir) / 'patient_results.csv'
    with open(ld, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=lesion_fields); w.writeheader(); w.writerows(lesion_rows)
    with open(pd_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=patient_fields); w.writeheader(); w.writerows(patient_rows)

    # Summary
    ok = [r for r in lesion_rows if r['status'] == 'ok']
    probs = [float(r['prob']) for r in ok]
    preds = [int(r['pred']) for r in ok]
    n = len(probs)
    hits = sum(preds)
    metric = hits/n if adapter.__class__.__name__.startswith('HCC') else (n-hits)/n
    metric_name = 'sensitivity' if adapter.__class__.__name__.startswith('HCC') else 'specificity'

    def wilson(p,n,z=1.96):
        if n==0: return (0,0)
        d=1+z*z/n; c=(p+z*z/(2*n))/d; h=z*((p*(1-p)/n+z*z/(4*n*n))**0.5)/d
        return max(0,c-h),min(1,c+h)
    lo, hi = wilson(metric, n)
    sum_path = Path(args.out_dir) / 'lesion_summary.txt'
    with open(sum_path, 'w') as f:
        f.write(f'Cohort: {args.cohort}\n')
        f.write(f'N patients: {len(patients)}  N lesions: {n}\n')
        f.write(f'{metric_name} (thr {args.threshold}): {metric:.4f}  Wilson 95% CI ({lo:.4f}, {hi:.4f})\n')
        f.write(f'Mean p(HCC): {np.mean(probs):.4f}  Median: {np.median(probs):.4f}\n')
        f.write(f'Elapsed: {time.time()-t_start:.1f}s\n')
    print(f'\n→ {ld}\n→ {pd_path}\n→ {sum_path}')


if __name__ == '__main__':
    main()
