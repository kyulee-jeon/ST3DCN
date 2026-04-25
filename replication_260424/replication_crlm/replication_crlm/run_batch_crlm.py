"""Per-patient batch inference over all 197 CRLM patients.

Each CRLM patient is ground-truth NEGATIVE for HCC (metastasis, not HCC) →
primary outcome: specificity = TN / (TN + FP) where FP = prob ≥ 0.8.

Outputs:
  - batch_results.csv: per-patient probability + metadata + failure reason
  - batch_summary.txt: specificity + Wilson 95% CI
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

from st3dcn_pipeline_crlm import process_patient, load_manifest
from ST3DCN_Utils import multi_scale_get_model_DCN

OUT_DIR = '/home/ubuntu/hcc_workspace/replication_crlm'
RESULTS_CSV = f'{OUT_DIR}/batch_results.csv'
SUMMARY_TXT = f'{OUT_DIR}/batch_summary.txt'
WEIGHTS = '/home/ubuntu/hcc_workspace/ST3DCN/weights/ST3DCN_Model.h5'
THRESHOLD = 0.8

print('==> loading manifest + model')
manifest = load_manifest()
patients = sorted(manifest['PatientID'].unique())
print(f'   {len(patients)} patients')

model = multi_scale_get_model_DCN(width=70, height=70, depth=70,
                                  batch_size=16, factor=8, num_class=2)
model.load_weights(WEIGHTS)

fields = ['patient_id', 'status', 'prob', 'pred',
          'portal_vein_median_hu', 'sop_matched', 'tumor_segments_present',
          'ct_shape', 'bbox', 'tumor_nz_slices', 'seg_uid', 'ct_uid',
          'pix_spacing_mm', 'slice_thick', 'elapsed_s', 'error']

results = []
t_start = time.time()
for i, pid in enumerate(patients, 1):
    t0 = time.time()
    row = {k: '' for k in fields}
    row['patient_id'] = pid
    try:
        rec = process_patient(pid, manifest)
        x = rec['crop'][None, :, :, :, None]
        prob = float(model.predict(x, verbose=0).ravel()[0])
        row.update({
            'status': 'ok',
            'prob': f'{prob:.6f}',
            'pred': int(prob >= THRESHOLD),
            'portal_vein_median_hu': (f"{rec['portal_vein_median_hu']:.1f}"
                                      if rec['portal_vein_median_hu'] is not None else ''),
            'sop_matched': f"{rec['sop_matched'][0]}/{rec['sop_matched'][1]}",
            'tumor_segments_present': ','.join(str(s) for s in rec['tumor_segments_present']),
            'ct_shape': 'x'.join(str(s) for s in rec['ct_shape']),
            'bbox': ','.join(str(v) for v in rec['bbox']),
            'tumor_nz_slices': rec['tumor_nz_slices'],
            'seg_uid': rec['seg_uid'],
            'ct_uid': rec['ct_uid'],
            'pix_spacing_mm': ','.join(f'{v:.4f}' for v in rec['pix_spacing_mm']),
            'slice_thick': f"{rec['slice_thick']:.2f}",
        })
        elapsed = time.time() - t0
        row['elapsed_s'] = f'{elapsed:.2f}'
        print(f'[{i:3}/{len(patients)}] {pid}  p={prob:.4f}  '
              f'pv_hu={rec["portal_vein_median_hu"]}  tumors={rec["tumor_segments_present"]}  '
              f't={elapsed:.1f}s')
    except Exception as e:
        row['status'] = 'fail'
        row['error'] = f'{type(e).__name__}: {e}'
        row['elapsed_s'] = f'{time.time() - t0:.2f}'
        print(f'[{i:3}/{len(patients)}] {pid}  FAIL  {row["error"]}')
    results.append(row)

with open(RESULTS_CSV, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(results)
print(f'\nresults → {RESULTS_CSV}')

ok = [r for r in results if r['status'] == 'ok']
n_ok = len(ok)
fp = sum(1 for r in ok if int(r['pred']) == 1)   # false positives (pred HCC, truth negative)
tn = n_ok - fp
spec = tn / n_ok if n_ok else 0.0

def wilson(p, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    denom = 1 + z**2 / n
    centre = (p + z**2 / (2*n)) / denom
    half = z * ((p*(1-p)/n + z**2/(4*n**2))**0.5) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))

lo, hi = wilson(spec, n_ok)
total = len(results)
failures = total - n_ok
probs = [float(r['prob']) for r in ok]

summary = f"""=== ST3DCN specificity on CRLM (portovenous-only, pre-op) ===
Total patients in manifest:   {total}
Successfully processed:       {n_ok}
Failures:                     {failures}

Threshold for binary HCC pred: {THRESHOLD}
False positives (pred==1):    {fp} / {n_ok}
True negatives  (pred==0):    {tn} / {n_ok}
SPECIFICITY (pu, raw):        {spec:.4f}
Wilson 95% CI:                ({lo:.4f}, {hi:.4f})

Probability distribution:
  mean   {np.mean(probs):.4f}
  median {np.median(probs):.4f}
  p05    {np.percentile(probs, 5):.4f}
  p25    {np.percentile(probs, 25):.4f}
  p75    {np.percentile(probs, 75):.4f}
  p95    {np.percentile(probs, 95):.4f}
"""
print(summary)
with open(SUMMARY_TXT, 'w') as f:
    f.write(summary)
    if failures:
        f.write('\n--- failures ---\n')
        from collections import Counter
        err_counter = Counter(r['error'].split(':')[0] for r in results if r['status']=='fail')
        for e, c in err_counter.most_common():
            f.write(f'  {e}: {c}\n')
print(f'summary → {SUMMARY_TXT}')
print(f'total wall-time: {time.time() - t_start:.1f}s')
