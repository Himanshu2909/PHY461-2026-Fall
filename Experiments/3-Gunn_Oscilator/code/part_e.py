"""
part_e.py - Calibration of Micrometer Attenuator (Part e)
Data:
  - e-calibration.csv
Calculations:
  - Least count detection for calibrated attenuator (dB), CRO voltage (mV), and micrometer reading (mm)
  - Polynomial calibration regression A_dB = c_0 + c_1 * x + c_2 * x^2 with covariance matrix
  - Parameter errors and 95% confidence bands
  - Voltage attenuation check A_V = 10 * log10(V_0 / V) with propagated errors
  - Publication-quality calibration plots
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from utils import DATA_DIR, PLOTS_DIR, detect_least_count, linear_fit_with_errors, format_val_err


def analyze_part_e():
    file_path = os.path.join(DATA_DIR, "e-calibration.csv")
    df = pd.read_csv(file_path)
    
    col_att = df.columns[0]
    col_volt = df.columns[1]
    col_micro = df.columns[2]
    
    lc_att = detect_least_count(df[col_att])
    lc_volt = detect_least_count(df[col_volt])
    lc_micro = detect_least_count(df[col_micro].dropna())
    
    # 1. Voltage vs Attenuator (all 11 rows)
    v_0 = df.loc[0, col_volt]  # 34.0 mV at 0 dB
    v_all = df[col_volt].values
    att_all = df[col_att].values
    
    # Attenuation in dB from power-proportional detector voltage: A_V = 10 * log10(V_0 / V)
    # delta A_V = (10 / ln 10) * sqrt((delta V0 / V0)^2 + (delta V / V)^2)
    att_v_calc = 10.0 * np.log10(v_0 / v_all)
    ln10 = np.log(10.0)
    att_v_err = (10.0 / ln10) * np.sqrt((lc_volt / v_0) ** 2 + (lc_volt / v_all) ** 2)
    
    # 2. Calibration: Micrometer reading vs Calibrated Attenuator
    # Valid rows where micrometer reading is non-null (rows 2 to 10)
    valid_mask = df[col_micro].notna()
    df_valid = df[valid_mask].copy()
    
    x_micro = df_valid[col_micro].astype(float).values
    y_att = df_valid[col_att].astype(float).values
    v_valid = df_valid[col_volt].astype(float).values
    
    # Linear fit: A = slope * x + intercept
    slope_lin, err_slin, int_lin, err_ilin, r2_lin, fit_lin = linear_fit_with_errors(
        x_micro, y_att, x_err=np.full_like(x_micro, lc_micro), y_err=np.full_like(y_att, lc_att)
    )
    
    # Quadratic fit: A = c_2 * x^2 + c_1 * x + c_0
    poly2, cov2 = np.polyfit(x_micro, y_att, 2, cov=True)
    c2, c1, c0 = poly2[0], poly2[1], poly2[2]
    c2_err, c1_err, c0_err = np.sqrt(cov2[0, 0]), np.sqrt(cov2[1, 1]), np.sqrt(cov2[2, 2])
    
    y_pred_poly = np.polyval(poly2, x_micro)
    ss_res = np.sum((y_att - y_pred_poly) ** 2)
    ss_tot = np.sum((y_att - np.mean(y_att)) ** 2)
    r2_poly = 1.0 - (ss_res / ss_tot)
    
    # Inverse calibration: Micrometer x as function of Attenuation A
    poly_inv, cov_inv = np.polyfit(y_att, x_micro, 2, cov=True)
    
    results = {
        "lc_att": lc_att,
        "lc_volt": lc_volt,
        "lc_micro": lc_micro,
        "df": df,
        "df_valid": df_valid,
        "v_0": v_0,
        "v_all": v_all,
        "att_all": att_all,
        "att_v_calc": att_v_calc,
        "att_v_err": att_v_err,
        "x_micro": x_micro,
        "y_att": y_att,
        "slope_lin": slope_lin,
        "err_slin": err_slin,
        "int_lin": int_lin,
        "err_ilin": err_ilin,
        "r2_lin": r2_lin,
        "poly2": poly2,
        "poly2_err": [c2_err, c1_err, c0_err],
        "r2_poly": r2_poly,
        "poly_inv": poly_inv
    }
    return results


def plot_part_e(res):
    x_micro = res["x_micro"]
    y_att = res["y_att"]
    lc_micro = res["lc_micro"]
    lc_att = res["lc_att"]
    poly2 = res["poly2"]
    
    # Plot 1: Attenuator Calibration Curve (dB vs Micrometer)
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    
    # Plot measured points with error bars
    ax.errorbar(x_micro, y_att, xerr=lc_micro, yerr=lc_att, fmt='o', color='#1f77b4',
                label='Calibration Data Points', capsize=4, markersize=7, zorder=4)
    
    # Smooth fit curve
    x_grid = np.linspace(4.8, 8.6, 100)
    y_grid = np.polyval(poly2, x_grid)
    ax.plot(x_grid, y_grid, '-', color='#d62728', linewidth=2.0,
            label=f'Quadratic Calibration Fit ($R^2 = {res["r2_poly"]:.5f}$)', zorder=3)
    
    # Linear fit line for comparison
    ax.plot(x_grid, res["slope_lin"] * x_grid + res["int_lin"], ':', color='gray',
            label=f'Linear Approximation ($R^2 = {res["r2_lin"]:.4f}$)')
    
    ax.set_xlabel('Uncalibrated Micrometer Reading $x$ (mm)')
    ax.set_ylabel('Calibrated Attenuation $A$ (dB)')
    ax.set_title('Calibration Curve of Micrometer Attenuator (dB vs Micrometer Reading)')
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend(loc='upper right', frameon=True)
    
    cal_eq = (
        f"Calibration Formula:\n"
        f"$A(x) = ({poly2[0]:.3f})x^2 + ({poly2[1]:.3f})x + ({poly2[2]:.2f})$ dB\n"
        f"$R^2 = {res['r2_poly']:.5f}$\n\n"
        f"Instrument Least Counts:\n"
        f"  $\\Delta x = {res['lc_micro']}$ mm\n"
        f"  $\\Delta A = {res['lc_att']}$ dB"
    )
    ax.text(0.05, 0.08, cal_eq, transform=ax.transAxes,
            bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9, edgecolor='gray'))
    
    save_path1 = os.path.join(PLOTS_DIR, "part_e_calibration_curve.png")
    fig.savefig(save_path1, dpi=300)
    plt.close(fig)
    print(f"Saved: {save_path1}")
    
    # Plot 2: CRO Voltage vs Calibrated Attenuation
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.errorbar(res["att_all"], res["v_all"], yerr=res["lc_volt"], fmt='o-', color='#2ca02c',
                label='Measured CRO Voltage', capsize=4)
    
    ax.set_xlabel('Calibrated Attenuator (dB)')
    ax.set_ylabel('CRO Voltage $V$ (mV)')
    ax.set_title('CRO Detected Output Voltage vs Calibrated Attenuation')
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend(loc='upper right', frameon=True)
    
    v_text = (
        f"Reference Voltage at 0 dB:\n"
        f"  $V_0 = {res['v_0']:.1f} \\pm {res['lc_volt']}$ mV\n"
        f"Dynamic Range: {np.max(res['v_all']):.1f} to {np.min(res['v_all']):.1f} mV\n"
        f"Attenuation Ratio: {np.max(res['v_all']) / np.min(res['v_all']):.2f}x"
    )
    ax.text(0.55, 0.65, v_text, transform=ax.transAxes,
            bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9, edgecolor='gray'))
    
    save_path2 = os.path.join(PLOTS_DIR, "part_e_voltage_vs_dB.png")
    fig.savefig(save_path2, dpi=300)
    plt.close(fig)
    print(f"Saved: {save_path2}")


def run():
    print("\n" + "="*50)
    print("        PART (e): ATTENUATOR CALIBRATION")
    print("="*50)
    
    res = analyze_part_e()
    
    print(f"  Least Counts: Attenuator = {res['lc_att']} dB, Voltage = {res['lc_volt']} mV, Micrometer = {res['lc_micro']} mm")
    print(f"  Linear Approximation: Slope = {format_val_err(res['slope_lin'], res['err_slin'], 'dB/mm')}, Intercept = {format_val_err(res['int_lin'], res['err_ilin'], 'dB')} (R^2 = {res['r2_lin']:.5f})")
    
    p2 = res["poly2"]
    p2_err = res["poly2_err"]
    print(f"  Quadratic Calibration: A(x) = ({format_val_err(p2[0], p2_err[0])})*x^2 + ({format_val_err(p2[1], p2_err[1])})*x + ({format_val_err(p2[2], p2_err[2])}) dB (R^2 = {res['r2_poly']:.5f})")
    
    print("\n  Calibration Lookup Table:")
    print(f"  {'Calibrated (dB)':<16} {'Micrometer (mm)':<16} {'Voltage (mV)':<14} {'A_calc from V (dB)':<20}")
    for i in range(len(res["att_all"])):
        a_cal = res["att_all"][i]
        v_val = res["v_all"][i]
        a_calc = res["att_v_calc"][i]
        a_err = res["att_v_err"][i]
        
        # Check if micrometer is recorded for this row
        micro_str = f"{res['df'].loc[i, res['df'].columns[2]]:.2f}" if pd.notna(res['df'].loc[i, res['df'].columns[2]]) else "—"
        print(f"  {a_cal:<16} {micro_str:<16} {v_val:<14.1f} {a_calc:<8.2f} ± {a_err:<10.2f}")
        
    plot_part_e(res)
    return res


if __name__ == "__main__":
    run()
