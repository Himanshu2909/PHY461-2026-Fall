"""
part_b.py - Analysis of Gunn as a Source of Microwaves (Part b)
Data:
  - b-i-microwaves-6.7V.csv
  - b-i-microwaves-8V.csv
Calculations:
  - Least count detection for micrometer reading and frequency
  - Mechanical cavity tuning sensitivity df/dx with uncertainty
  - Electronic tuning shift Delta f = f(8V) - f(6.7V) and pushing factor df/dV
  - Theoretical free-space wavelength lambda_0 and guide wavelength lambda_g with propagated errors
  - Publication-quality tuning and dispersion plots
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from utils import DATA_DIR, PLOTS_DIR, detect_least_count, linear_fit_with_errors, format_val_err

# Physical constants
C_CMS = 2.99792458e10   # Speed of light in cm/s
A_WAVEGUIDE = 2.286     # WR-90 waveguide broad dimension in cm
LAMBDA_C = 2.0 * A_WAVEGUIDE  # Cutoff wavelength = 4.572 cm


def calculate_wavelengths(f_ghz, f_err_ghz):
    """
    Computes free-space wavelength lambda_0 (cm) and theoretical guide wavelength lambda_g (cm)
    with analytical error propagation.
    """
    f_hz = f_ghz * 1.0e9
    f_err_hz = f_err_ghz * 1.0e9
    
    # lambda_0 = c / f
    lam_0 = C_CMS / f_hz
    lam_0_err = lam_0 * (f_err_hz / f_hz)
    
    # lambda_g = lam_0 / sqrt(1 - (lam_0 / lam_c)^2)
    cutoff_ratio = lam_0 / LAMBDA_C
    denom = np.sqrt(1.0 - cutoff_ratio ** 2)
    lam_g = lam_0 / denom
    
    # d(lam_g)/d(lam_0) = (lam_g / lam_0)^3
    lam_g_err = ((lam_g / lam_0) ** 3) * lam_0_err
    
    return lam_0, lam_0_err, lam_g, lam_g_err


def analyze_part_b():
    file_67 = os.path.join(DATA_DIR, "b-i-microwaves-6.7V.csv")
    file_80 = os.path.join(DATA_DIR, "b-i-microwaves-8V.csv")
    
    df_67 = pd.read_csv(file_67)
    df_80 = pd.read_csv(file_80)
    
    # Clean column names
    col_x = df_67.columns[0]
    col_f = df_67.columns[1]
    
    lc_x = detect_least_count(df_67[col_x])
    lc_f = detect_least_count(df_67[col_f])
    
    x_67 = df_67[col_x].values
    f_67 = df_67[col_f].values
    
    x_80 = df_80[col_x].values
    f_80 = df_80[col_f].values
    
    # 1. Tuning range
    range_67 = np.max(f_67) - np.min(f_67)
    range_80 = np.max(f_80) - np.min(f_80)
    
    # 2. Linear fits for mechanical tuning sensitivity df/dx
    slope_67, err_s67, int_67, _, r2_67, fit_67 = linear_fit_with_errors(x_67, f_67, y_err=np.full_like(f_67, lc_f))
    slope_80, err_s80, int_80, _, r2_80, fit_80 = linear_fit_with_errors(x_80, f_80, y_err=np.full_like(f_80, lc_f))
    
    # 3. Electronic tuning shift Delta f_bias at each micrometer position
    delta_f = (f_80 - f_67) * 1000.0  # in MHz
    delta_f_err = np.sqrt(2) * (lc_f * 1000.0)  # in MHz
    mean_delta_f = np.mean(delta_f)
    std_delta_f = np.std(delta_f, ddof=1) if len(delta_f) > 1 else 0.0
    
    delta_v = 8.0 - 6.7  # 1.3 V
    pushing_factor = mean_delta_f / delta_v  # MHz / V
    pushing_factor_err = delta_f_err / delta_v
    
    # 4. Wavelengths calculation
    lam0_67, lam0_err_67, lamg_67, lamg_err_67 = calculate_wavelengths(f_67, lc_f)
    lam0_80, lam0_err_80, lamg_80, lamg_err_80 = calculate_wavelengths(f_80, lc_f)
    
    results = {
        "lc_x": lc_x,
        "lc_f": lc_f,
        "x": x_67,
        "f_67": f_67,
        "f_80": f_80,
        "range_67": range_67,
        "range_80": range_80,
        "slope_67": slope_67,
        "err_s67": err_s67,
        "r2_67": r2_67,
        "slope_80": slope_80,
        "err_s80": err_s80,
        "r2_80": r2_80,
        "delta_f": delta_f,
        "delta_f_err": delta_f_err,
        "mean_delta_f": mean_delta_f,
        "pushing_factor": pushing_factor,
        "pushing_factor_err": pushing_factor_err,
        "lam0_67": lam0_67,
        "lam0_err_67": lam0_err_67,
        "lamg_67": lamg_67,
        "lamg_err_67": lamg_err_67,
        "lam0_80": lam0_80,
        "lam0_err_80": lam0_err_80,
        "lamg_80": lamg_80,
        "lamg_err_80": lamg_err_80
    }
    return results


def plot_part_b(res):
    x = res["x"]
    f_67 = res["f_67"]
    f_80 = res["f_80"]
    lc_f = res["lc_f"]
    
    # Plot 1: Tuning curves
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.errorbar(x, f_67, yerr=lc_f, fmt='o-', color='#1f77b4', label=f'Gunn Bias = 6.7 V ($df/dx = {res["slope_67"]:.3f}$ GHz/mm)', capsize=4)
    ax.errorbar(x, f_80, yerr=lc_f, fmt='s--', color='#d62728', label=f'Gunn Bias = 8.0 V ($df/dx = {res["slope_80"]:.3f}$ GHz/mm)', capsize=4)
    
    ax.set_xlabel('Cavity Tuning Micrometer Reading (mm)')
    ax.set_ylabel('Oscillation Frequency $f$ (GHz)')
    ax.set_title('Gunn Oscillator Cavity & Electronic Frequency Tuning')
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend(loc='upper right', frameon=True)
    
    info_text = (
        f"Tuning Range: {res['range_67']:.3f} GHz\n"
        f"Bias Shift: $\\Delta f = +{res['mean_delta_f']:.1f}$ MHz\n"
        f"Pushing Factor: ${res['pushing_factor']:.2f} \\pm {res['pushing_factor_err']:.2f}$ MHz/V"
    )
    ax.text(0.05, 0.15, info_text, transform=ax.transAxes,
            bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.85, edgecolor='gray'))
    
    save_path1 = os.path.join(PLOTS_DIR, "part_b_tuning_curves.png")
    fig.savefig(save_path1, dpi=300)
    plt.close(fig)
    print(f"Saved: {save_path1}")
    
    # Plot 2: Theoretical Wavelengths
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.plot(x, res["lam0_67"], 'o-', color='#2ca02c', label=r'Free Space $\lambda_0$ (6.7 V)')
    ax.plot(x, res["lamg_67"], 's-', color='#9467bd', label=r'Guide Wavelength $\lambda_g$ (6.7 V)')
    ax.plot(x, res["lamg_80"], '^--', color='#8c564b', label=r'Guide Wavelength $\lambda_g$ (8.0 V)')
    
    ax.axhline(LAMBDA_C, color='gray', linestyle=':', label=r'Cutoff $\lambda_c = 4.572$ cm')
    
    ax.set_xlabel('Cavity Tuning Micrometer Reading (mm)')
    ax.set_ylabel('Wavelength (cm)')
    ax.set_title(r'Calculated Wavelengths ($\lambda_0$ and $\lambda_g$) vs Cavity Micrometer')
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend(loc='lower right', frameon=True)
    
    save_path2 = os.path.join(PLOTS_DIR, "part_b_wavelengths.png")
    fig.savefig(save_path2, dpi=300)
    plt.close(fig)
    print(f"Saved: {save_path2}")


def run():
    print("\n" + "="*50)
    print("      PART (b): GUNN MICROWAVE SOURCE TUNING")
    print("="*50)
    
    res = analyze_part_b()
    
    print(f"  Least Counts: Micrometer = {res['lc_x']} mm, Frequency = {res['lc_f']} GHz")
    print(f"  Mechanical Tuning Sensitivity (6.7 V): {format_val_err(res['slope_67'], res['err_s67'], 'GHz/mm')} (R^2 = {res['r2_67']:.5f})")
    print(f"  Mechanical Tuning Sensitivity (8.0 V): {format_val_err(res['slope_80'], res['err_s80'], 'GHz/mm')} (R^2 = {res['r2_80']:.5f})")
    print(f"  Total Tuning Range: {res['range_67']:.3f} GHz (6.7 V), {res['range_80']:.3f} GHz (8.0 V)")
    print(f"  Electronic Frequency Pulling (8V - 6.7V): {res['mean_delta_f']:.1f} ± {res['delta_f_err']:.1f} MHz")
    print(f"  Frequency Pushing Figure: {format_val_err(res['pushing_factor'], res['pushing_factor_err'], 'MHz/V')}")
    print(f"  Representative Guide Wavelength (at 5 mm, 10.625 GHz): {format_val_err(res['lamg_67'][0], res['lamg_err_67'][0], 'cm')}")
    print(f"  Representative Guide Wavelength (at 10 mm, 9.425 GHz): {format_val_err(res['lamg_67'][-1], res['lamg_err_67'][-1], 'cm')}")
    
    plot_part_b(res)
    return res


if __name__ == "__main__":
    run()
