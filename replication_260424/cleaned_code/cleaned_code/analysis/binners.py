"""Binning functions for per-tag heterogeneity analysis.

Numeric tags are discretized into ranges; categorical tags are normalized
and left as-is. Every bin key must be a stable string for GroupBy.
"""
import re

import numpy as np
import pandas as pd


def try_float(v):
    if v is None: return None
    try:
        if pd.isna(v): return None
    except (TypeError, ValueError): pass
    try:
        x = float(str(v).strip("'\""))
        return None if np.isnan(x) else x
    except (ValueError, TypeError, AttributeError):
        return None


def parse_pixspacing_mean(v):
    if pd.isna(v): return None
    nums = [float(x) for x in re.findall(r'\d+\.?\d*', str(v))]
    return None if not nums else sum(nums) / len(nums)


def _num_canon(s):
    try:
        x = float(s)
        return str(int(x)) if x == int(x) else f'{x:g}'
    except (ValueError, TypeError):
        return s


def normalize_cat(v):
    if pd.isna(v): return '(missing)'
    s = str(v).strip().strip("'\"")
    return s if s and s.lower() != 'nan' else '(missing)'


def normalize_multival_cat(v):
    """Handle DICOM MultiValue like '[40, 40]' → '40' if uniform."""
    if pd.isna(v): return '(missing)'
    s = str(v).strip().strip("'\"")
    if not s or s.lower() == 'nan': return '(missing)'
    m = re.match(r'\[\s*([^,\]]+)\s*(?:,\s*([^,\]]+))?\s*\]', s)
    if m:
        parts = [p.strip() for p in m.groups() if p]
        if parts and all(p == parts[0] for p in parts):
            return _num_canon(parts[0])
    return _num_canon(s)


def kernel_family(v):
    if pd.isna(v): return '(missing)'
    s = str(v).upper().strip().strip("'\"")
    if not s: return '(missing)'
    for fam in ('STANDARD', 'SOFT', 'B30', 'B40', 'B50', 'FC', 'BONE', 'LUNG'):
        if fam in s: return fam
    return s[:8]


def _bin_num(v, edges):
    x = try_float(v)
    if x is None: return '(missing)'
    for lo, hi in edges:
        if lo <= x < hi:
            return f'[{lo}, {hi})'
    return f'[{edges[-1][1]}+)'


def bin_tube_current(v):    return _bin_num(v, [(0,200),(200,300),(300,400),(400,500),(500,10_000)])
def bin_recondiam(v):       return _bin_num(v, [(0,320),(320,360),(360,400),(400,500)])
def bin_tcw(v):             return _bin_num(v, [(0,20),(20,40),(40,80),(80,160)])
def bin_tablespeed(v):      return _bin_num(v, [(0,20),(20,40),(40,60),(60,100)])
def bin_exposure_uas(v):    return _bin_num(v, [(0,150_000),(150_000,250_000),
                                                 (250_000,350_000),(350_000,450_000),
                                                 (450_000,1_000_000),(1_000_000,10_000_000)])

def bin_pixspacing(v):
    x = parse_pixspacing_mean(v)
    if x is None: return '(missing)'
    edges = [(0.50, 0.65), (0.65, 0.75), (0.75, 0.85), (0.85, 1.00)]
    for lo, hi in edges:
        if lo <= x < hi:
            return f'[{lo:.2f}, {hi:.2f})'
    return '(out-of-range)'


def bin_pitch(v):
    x = try_float(v); return '(missing)' if x is None else f'{round(x,3)}'


def bin_revtime(v):
    x = try_float(v); return '(missing)' if x is None else f'{round(x,2)}'


def bin_spacing(v):
    x = try_float(v)
    if x is None: return '(missing)'
    edges = [(0, 1.5), (1.5, 3.0), (3.0, 5.0), (5.0, 8.0), (8.0, 20.0)]
    for lo, hi in edges:
        if lo <= x < hi:
            return f'[{lo:.1f}, {hi:.1f})'
    return '[20+)'


BINNERS = {
    'cat':              normalize_cat,
    'cat_multival':     normalize_multival_cat,
    'kernel_family':    kernel_family,
    'num_tube_current': bin_tube_current,
    'num_recondiam':    bin_recondiam,
    'num_pixspacing':   bin_pixspacing,
    'num_pitch':        bin_pitch,
    'num_tcw':          bin_tcw,
    'num_tablespeed':   bin_tablespeed,
    'num_revtime':      bin_revtime,
    'num_exposure_uas': bin_exposure_uas,
    'num_spacing':      bin_spacing,
}
