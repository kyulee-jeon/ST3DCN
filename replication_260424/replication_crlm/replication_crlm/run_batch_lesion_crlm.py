"""Per-lesion batch inference over all 197 CRLM patients.

Each tumor segment (Tumor_1..Tumor_5) = one row. All lesions are ground-truth
NEGATIVE for HCC → lesion-level specificity.
"""
import os, sys, time, csv
os.environ['CUDA_VISIBLE_DEVICES'] = '2'
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

import numpy as np
import tensorflow as tf
for g in tf.config.list_physical_devices('GPU'):
    tf.config.experimental.set_memory_growth(g, True)

sys.path.insert(0, '/home/ubuntu/hcc_workspace/replication_crlm')
sys.path.insert(0, '/home/ubuntu/hcc_workspace/ST3DCN')

from st3dcn_pipeline_crlm import process_patient_per_lesion, load_manifest
from ST3DCN_Utils import multi_scale_get_model_DCN

OUT_DIR = '/home/ubuntu/hcc_workspace/replication_crlm'
RESULTS_CSV = f'{OUT_DIR}/lesion_results.csv'
SUMMARY_TXT = f'{OUT_DIR}/lesion_summary.txt'
WEIGHTS = '/home/ubuntu/hcc_workspace/ST3DCN/weights/ST3DCN_Model.h5'
THRESHOLD = 0.8
MIN_VOX = 50

manifest = load_manifest()
patients = sorted(manifest['PatientID'].unique())
model = multi_scale_get_model_DCN(width=70, height=70, depth=70,
                                  batch_size=16, factor=8, num_class=2)
model.load_weights(WEIGHTS)

fields = ['patient_id', 'lesion_idx', 'segment_number', 'segment_label',
          'voxel_count', 'prob', 'pred', 'portal_vein_median_hu',
          'sop_matched', 'ct_shape', 'bbox', 'seg_uid', 'ct_uid',
          'pix_spacing_mm', 'slice_thick', 'n_lesions_patient',
          'status', 'error']

rows = []
t_start = time.time()
for i, pid in enumerate(patients, 1):
    t0 = time.time()
    try:
        rec = process_patient_per_lesion(pid, manifest, min_voxels=MIN_VOX)
        probs_here = []
        for les in rec['lesions']:
            x = les['crop'][None, :, :, :, None]
            prob = float(model.predict(x, verbose=0).ravel()[0])
            probs_here.append(prob)
            rows.append({
                'patient_id': pid,
                'lesion_idx': les['lesion_idx'],
                'segment_number': les['segment_number'],
                'segment_label': les['segment_label'],
                'voxel_count': les['voxel_count'],
                'prob': f'{prob:.6f}',
                'pred': int(prob >= THRESHOLD),
                'portal_vein_median_hu': (f"{rec['portal_vein_median_hu']:.1f}"
                                          if rec['portal_vein_median_hu'] is not None else ''),
                'sop_matched': f"{rec['sop_matched'][0]}/{rec['sop_matched'][1]}",
                'ct_shape': 'x'.join(str(s) for s in rec['ct_shape']),
                'bbox': ','.join(str(v) for v in les['bbox']),
                'seg_uid': rec['seg_uid'],
                'ct_uid': rec['ct_uid'],
                'pix_spacing_mm': ','.join(f'{v:.4f}' for v in rec['pix_spacing_mm']),
                'slice_thick': f"{rec['slice_thick']:.2f}",
                'n_lesions_patient': len(rec['lesions']),
                'status': 'ok', 'error': '',
            })
        elapsed = time.time() - t0
        print(f'[{i:3}/{len(patients)}] {pid}  n_les={len(rec["lesions"])}  '
              f'probs={[f"{p:.3f}" for p in probs_here]}  t={elapsed:.1f}s')
    except Exception as e:
        rows.append({
            'patient_id': pid, 'lesion_idx': -1, 'segment_number': '', 'segment_label': '',
            'voxel_count': 0, 'prob': '', 'pred': '', 'portal_vein_median_hu': '',
            'sop_matched': '', 'ct_shape': '', 'bbox': '', 'seg_uid': '', 'ct_uid': '',
            'pix_spacing_mm': '', 'slice_thick': '', 'n_lesions_patient': 0,
            'status': 'fail', 'error': f'{type(e).__name__}: {e}',
        })
        print(f'[{i:3}/{len(patients)}] {pid}  FAIL  {type(e).__name__}: {e}')

with open(RESULTS_CSV, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(rows)

ok = [r for r in rows if r['status'] == 'ok']
n_patients_ok = len(set(r['patient_id'] for r in ok))
n_les = len(ok)
fp_les = sum(1 for r in ok if int(r['pred']) == 1)
tn_les = n_les - fp_les
spec_les = tn_les / n_les if n_les else 0.0

def wilson(p, n, z=1.96):
    if n == 0: return (0.0, 0.0)
    denom = 1 + z**2 / n
    centre = (p + z**2 / (2*n)) / denom
    half = z * ((p*(1-p)/n + z**2/(4*n**2))**0.5) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))

lo, hi = wilson(spec_les, n_les)
probs = [float(r['prob']) for r in ok]

summary = f"""=== CRLM lesion-level specificity (portovenous-only, pre-op) ===
Patients processed:         {n_patients_ok}/{len(patients)}
Total lesions:              {n_les}
lesions / patient (mean):   {n_les/max(1,n_patients_ok):.2f}

Threshold for binary HCC pred: {THRESHOLD}
False positives (pred==1):  {fp_les} / {n_les}
True negatives  (pred==0):  {tn_les} / {n_les}
SPECIFICITY (lesion, raw):  {spec_les:.4f}
Wilson 95% CI:              ({lo:.4f}, {hi:.4f})

Probability distribution (lesion):
  mean   {np.mean(probs):.4f}
  median {np.median(probs):.4f}
  p05    {np.percentile(probs, 5):.4f}
  p95    {np.percentile(probs, 95):.4f}
"""
print(summary)
with open(SUMMARY_TXT, 'w') as f:
    f.write(summary)
print(f'summary → {SUMMARY_TXT}')
print(f'elapsed: {time.time()-t_start:.1f}s')
