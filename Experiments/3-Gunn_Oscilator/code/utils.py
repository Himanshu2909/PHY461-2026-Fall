"""
utils.py - Shared utilities for Gunn Oscillator Experiment Data Analysis
Provides:
  - Dynamic least-count detection from data strings / numbers
  - Linear and polynomial regression with full parameter uncertainty propagation
  - Weighted average calculation with errors
  - Consistent publication-quality plotting styles
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Output directories
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, "data")
PLOTS_DIR = os.path.join(BASE_DIR, "plots")
os.makedirs(PLOTS_DIR, exist_ok=True)

# Set global matplotlib style
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 14,
    'lines.linewidth': 1.8,
    'lines.markersize': 6,
    'errorbar.capsize': 3,
    'grid.alpha': 0.4,
    'grid.linestyle': '--',
    'figure.autolayout': True
})


def detect_least_count(values):
    """
    Infers the least count based on the smallest figure reported in a collection of values.
    Rule:
      If a reading is 34.12, least count is 0.01 (2 decimal places -> 10^-2).
      If readings are integers (e.g. 27, 501), least count is 1.0 (0 decimal places -> 10^0).
    """
    max_decimals = 0
    has_float = False
    
    for v in values:
        if pd.isna(v):
            continue
        s = str(v).strip()
        if '.' in s:
            has_float = True
            mantissa = s.split('e')[0].split('E')[0]
            dec_part = mantissa.split('.')[1]
            max_decimals = max(max_decimals, len(dec_part))
        else:
            try:
                val = float(s)
                if val != int(val):
                    has_float = True
                    dec_str = f"{val:.10f}".rstrip('0').split('.')[1]
                    max_decimals = max(max_decimals, len(dec_str))
            except ValueError:
                pass

    if not has_float or max_decimals == 0:
        return 1.0
    return 10.0 ** (-max_decimals)


def linear_fit_with_errors(x, y, x_err=None, y_err=None):
    """
    Performs linear regression y = slope * x + intercept with error estimation.
    Returns:
      slope, slope_err, intercept, intercept_err, r_squared, fit_fn
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    
    if y_err is not None and np.all(y_err > 0):
        weights = 1.0 / (y_err ** 2)
        poly, cov = np.polyfit(x, y, 1, w=np.sqrt(weights), cov=True)
    else:
        poly, cov = np.polyfit(x, y, 1, cov=True)
        
    slope, intercept = poly[0], poly[1]
    slope_err = np.sqrt(cov[0, 0])
    intercept_err = np.sqrt(cov[1, 1])
    
    # Calculate R-squared
    y_pred = slope * x + intercept
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r_squared = 1.0 - (ss_res / ss_tot) if ss_tot != 0 else 1.0
    
    return slope, slope_err, intercept, intercept_err, r_squared, lambda x_val: slope * x_val + intercept


def weighted_average(values, errors):
    """
    Computes weighted average and uncertainty of independent measurements.
    w_i = 1 / (sigma_i^2)
    mean = sum(w_i * x_i) / sum(w_i)
    err = 1 / sqrt(sum(w_i))
    """
    v = np.asarray(values, dtype=float)
    e = np.asarray(errors, dtype=float)
    weights = 1.0 / (e ** 2)
    sum_w = np.sum(weights)
    mean_val = np.sum(weights * v) / sum_w
    err_val = 1.0 / np.sqrt(sum_w)
    return mean_val, err_val


def format_val_err(val, err, unit="", precision=None):
    """
    Formats value +/- error with suitable precision.
    """
    if err is None or np.isnan(err) or err == 0:
        return f"{val:.3f} {unit}".strip()
    if precision is not None:
        return f"({val:.{precision}f} ± {err:.{precision}f}) {unit}".strip()
    if err < 1e-4:
        return f"({val:.6f} ± {err:.6f}) {unit}".strip()
    elif err < 1e-2:
        return f"({val:.4f} ± {err:.4f}) {unit}".strip()
    elif err < 1.0:
        return f"({val:.3f} ± {err:.3f}) {unit}".strip()
    else:
        return f"({val:.2f} ± {err:.2f}) {unit}".strip()
