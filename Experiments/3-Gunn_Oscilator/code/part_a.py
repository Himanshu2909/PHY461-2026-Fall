"""
part_a.py - Analysis of Gunn Diode I-V Characteristics (Part a)
Data:
  - a-i-IV-5mm.csv (Cavity micrometer = 5 mm)
  - a-ii-IV-8mm.csv (Cavity micrometer = 8 mm)
Calculations:
  - Least count detection for V and I
  - Ohmic region dynamic resistance R_ohmic with uncertainty
  - Threshold voltage V_0 and peak current I_peak with uncertainty
  - GaAs active region thickness d = V_0 / E_th and transit frequency nu = v / d
  - Oscillatory/domain transit region characterization (fluctuation spread)
  - Sustained oscillation region negative differential resistance R_diff
  - Publication-quality I-V curve plots
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from utils import DATA_DIR, PLOTS_DIR, detect_least_count, linear_fit_with_errors, format_val_err

# Physical constants
ETH_MANUAL = 36000.0  # V/cm (from Lab Manual Question 2)
ETH_LIT = 3200.0      # V/cm (Standard bulk GaAs threshold field)
V_DRIFT = 1.0e7       # cm/s (Electron peak drift velocity in GaAs)


def analyze_dataset(file_path, micrometer_label):
    df = pd.read_csv(file_path)
    
    # Extract least counts
    lc_v = detect_least_count(df["V"])
    lc_i = detect_least_count(df["I"])
    
    v = df["V"].values
    current = df["I"].values
    
    # 1. Ohmic region (initial linear rise up to V ~ 2.5 V)
    ohmic_mask = (v <= 2.5) & (df.index < 30)
    v_ohmic = v[ohmic_mask]
    i_ohmic = current[ohmic_mask]
    
    slope_g, g_err, int_g, int_err, r2_ohmic, fit_ohmic = linear_fit_with_errors(
        v_ohmic, i_ohmic, y_err=np.full_like(i_ohmic, lc_i)
    )
    # R_ohmic = 1000 / G (since I is in mA, V in V -> G in mA/V = mS -> R in Ohms)
    r_ohmic = 1000.0 / slope_g
    r_ohmic_err = (1000.0 / (slope_g ** 2)) * g_err
    
    # 2. Threshold region (plateau and peak current)
    # Peak current is 501 mA
    peak_i = np.max(current)
    peak_mask = (current == peak_i) & (df.index <= 32)
    v_plateau = v[peak_mask]
    v_thresh_mid = np.mean(v_plateau)
    v_thresh_err = (np.max(v_plateau) - np.min(v_plateau)) / 2.0
    if v_thresh_err == 0:
        v_thresh_err = lc_v
    v_thresh_drop = np.max(v_plateau)  # point just before current drops/fluctuates
    
    # 3. Slab thickness d = V_0 / E_th
    # In microns: d (cm) * 1e4 = d (um)
    d_cm_manual = v_thresh_mid / ETH_MANUAL
    d_um_manual = d_cm_manual * 1.0e4
    d_um_err_manual = (v_thresh_err / ETH_MANUAL) * 1.0e4
    
    d_cm_lit = v_thresh_mid / ETH_LIT
    d_um_lit = d_cm_lit * 1.0e4
    d_um_err_lit = (v_thresh_err / ETH_LIT) * 1.0e4
    
    # 4. Transit-time frequency nu = v / d
    nu_ghz_manual = (V_DRIFT / d_cm_manual) / 1.0e9
    nu_ghz_err_manual = nu_ghz_manual * (v_thresh_err / v_thresh_mid)
    
    nu_ghz_lit = (V_DRIFT / d_cm_lit) / 1.0e9
    nu_ghz_err_lit = nu_ghz_lit * (v_thresh_err / v_thresh_mid)
    
    # 5. Oscillatory / domain transit fluctuating region
    # Starts where current drops/fluctuates up to smooth high-bias region (V < 5.5 V)
    osc_mask = (df.index >= 33) & (v < 5.5)
    v_osc = v[osc_mask]
    i_osc = current[osc_mask]
    mean_i_osc = np.mean(i_osc)
    std_i_osc = np.std(i_osc, ddof=1)
    min_i_osc, max_i_osc = np.min(i_osc), np.max(i_osc)
    
    # 6. Sustained oscillation high-bias region (V >= 5.5 V)
    high_mask = (v >= 5.5)
    v_high = v[high_mask]
    i_high = current[high_mask]
    
    slope_ndr, ndr_err, int_ndr, _, r2_ndr, fit_ndr = linear_fit_with_errors(
        v_high, i_high, y_err=np.full_like(i_high, lc_i)
    )
    # R_diff = 1000 / slope_ndr in Ohms (negative differential resistance)
    r_diff = 1000.0 / slope_ndr
    r_diff_err = (1000.0 / (slope_ndr ** 2)) * ndr_err
    
    results = {
        "label": micrometer_label,
        "lc_v": lc_v,
        "lc_i": lc_i,
        "v_thresh_mid": v_thresh_mid,
        "v_thresh_err": v_thresh_err,
        "v_thresh_drop": v_thresh_drop,
        "peak_i": peak_i,
        "r_ohmic": r_ohmic,
        "r_ohmic_err": r_ohmic_err,
        "r2_ohmic": r2_ohmic,
        "d_um_manual": d_um_manual,
        "d_um_err_manual": d_um_err_manual,
        "nu_ghz_manual": nu_ghz_manual,
        "nu_ghz_err_manual": nu_ghz_err_manual,
        "d_um_lit": d_um_lit,
        "d_um_err_lit": d_um_err_lit,
        "nu_ghz_lit": nu_ghz_lit,
        "nu_ghz_err_lit": nu_ghz_err_lit,
        "mean_i_osc": mean_i_osc,
        "std_i_osc": std_i_osc,
        "min_i_osc": min_i_osc,
        "max_i_osc": max_i_osc,
        "r_diff": r_diff,
        "r_diff_err": r_diff_err,
        "r2_ndr": r2_ndr,
        "df": df,
        "fit_ohmic": fit_ohmic,
        "fit_ndr": fit_ndr,
        "v_ohmic": v_ohmic,
        "v_high": v_high
    }
    return results


def plot_part_a(res_5mm, res_8mm):
    # Plot individual 5mm
    for res in [res_5mm, res_8mm]:
        df = res["df"]
        v = df["V"].values
        i_curr = df["I"].values
        lbl = res["label"]
        
        fig, ax = plt.subplots(figsize=(8, 5.5))
        
        # 1. Ohmic region line
        ohmic_idx = df.index <= 32
        ax.plot(v[ohmic_idx], i_curr[ohmic_idx], 'o-', color='#1f77b4', label='Ohmic & Pre-threshold', zorder=4)
        
        # 2. Fluctuating / domain oscillation scatter points
        osc_idx = (df.index >= 33) & (v < 5.5)
        ax.scatter(v[osc_idx], i_curr[osc_idx], color='#d62728', marker='^', s=45, alpha=0.85,
                   label='Domain Oscillation / Instability', zorder=5)
        
        # 3. High-bias smooth region
        high_idx = (v >= 5.5)
        ax.plot(v[high_idx], i_curr[high_idx], 's-', color='#2ca02c', label='Sustained Oscillation (NDR)', zorder=4)
        
        # Shaded region for domain oscillation
        ax.axvspan(np.min(v[osc_idx]), np.max(v[osc_idx]), color='#ff7f0e', alpha=0.12,
                   label='Oscillatory Regime Window')
        
        # Mark threshold V_0
        ax.axvline(res["v_thresh_mid"], color='black', linestyle=':', linewidth=1.5,
                   label=f'$V_0 = {res["v_thresh_mid"]:.2f} \\pm {res["v_thresh_err"]:.2f}$ V')
        
        ax.set_xlabel('Bias Voltage $V$ (V)')
        ax.set_ylabel('Diode Current $I$ (mA)')
        ax.set_title(f'Gunn Diode I-V Characteristics ({lbl})')
        ax.set_xlim(-0.2, 8.6)
        ax.set_ylim(-10, 620)
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.legend(loc='upper left', frameon=True, framealpha=0.9, fontsize=9.5)
        
        # Text annotation of key values in clear upper-right region
        info_text = (
            f"$R_{{Ohmic}} = {res['r_ohmic']:.2f} \\pm {res['r_ohmic_err']:.2f}\\ \\Omega$\n"
            f"$V_0 = {res['v_thresh_mid']:.2f} \\pm {res['v_thresh_err']:.2f}$ V\n"
            f"$d_{{active}} = {res['d_um_manual']:.3f} \\pm {res['d_um_err_manual']:.3f}\\ \\mu$m\n"
            f"$R_{{diff}} = {res['r_diff']:.2f} \\pm {res['r_diff_err']:.2f}\\ \\Omega$"
        )
        ax.text(0.98, 0.97, info_text, transform=ax.transAxes, verticalalignment='top',
                horizontalalignment='right', fontsize=9.5,
                bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9, edgecolor='gray'))
        
        fname = f"part_a_IV_{lbl.replace(' ', '_').lower()}.png"
        save_path = os.path.join(PLOTS_DIR, fname)
        fig.savefig(save_path, dpi=300)
        plt.close(fig)
        print(f"Saved: {save_path}")
        
    # Comparison Plot (5mm vs 8mm)
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    ax.plot(res_5mm["df"]["V"], res_5mm["df"]["I"], '.-', color='#1f77b4', alpha=0.7, label='5 mm cavity micrometer')
    ax.plot(res_8mm["df"]["V"], res_8mm["df"]["I"], '.--', color='#d62728', alpha=0.7, label='8 mm cavity micrometer')
    ax.axvline(res_5mm["v_thresh_mid"], color='#1f77b4', linestyle=':', label='Threshold $V_0$ (5 mm)')
    ax.set_xlabel('Bias Voltage $V$ (V)')
    ax.set_ylabel('Current $I$ (mA)')
    ax.set_title('Gunn Diode I-V Characteristics Comparison (5 mm vs 8 mm)')
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend(loc='lower left', frameon=True)
    
    comp_path = os.path.join(PLOTS_DIR, "part_a_comparison.png")
    fig.savefig(comp_path, dpi=300)
    plt.close(fig)
    print(f"Saved: {comp_path}")


def run():
    print("\n" + "="*50)
    print("           PART (a): I-V CHARACTERISTICS")
    print("="*50)
    
    file_5mm = os.path.join(DATA_DIR, "a-i-IV-5mm.csv")
    file_8mm = os.path.join(DATA_DIR, "a-ii-IV-8mm.csv")
    
    res_5mm = analyze_dataset(file_5mm, "5mm")
    res_8mm = analyze_dataset(file_8mm, "8mm")
    
    for res in [res_5mm, res_8mm]:
        lbl = res["label"]
        print(f"\n--- Results for Cavity Setting: {lbl} ---")
        print(f"  Least Counts: Voltage = {res['lc_v']} V, Current = {res['lc_i']} mA")
        print(f"  Ohmic Dynamic Resistance: {format_val_err(res['r_ohmic'], res['r_ohmic_err'], 'Ohms')} (R^2 = {res['r2_ohmic']:.5f})")
        print(f"  Peak Current: {res['peak_i']} mA")
        print(f"  Threshold Voltage V_0: {format_val_err(res['v_thresh_mid'], res['v_thresh_err'], 'V')}")
        print(f"  Onset of NDM Drop: {res['v_thresh_drop']} V")
        print(f"  [Manual Eth = 36 kV/cm] Active Layer d: {format_val_err(res['d_um_manual'], res['d_um_err_manual'], 'microns')}")
        print(f"  [Manual Eth = 36 kV/cm] Predicted Transit Frequency: {format_val_err(res['nu_ghz_manual'], res['nu_ghz_err_manual'], 'GHz')}")
        print(f"  [Literature Eth = 3.2 kV/cm] Active Layer d: {format_val_err(res['d_um_lit'], res['d_um_err_lit'], 'microns')}")
        print(f"  [Literature Eth = 3.2 kV/cm] Predicted Transit Frequency: {format_val_err(res['nu_ghz_lit'], res['nu_ghz_err_lit'], 'GHz')}")
        print(f"  Domain Oscillation Regime: Mean I = {res['mean_i_osc']:.1f} ± {res['std_i_osc']:.1f} mA (Range: [{res['min_i_osc']}, {res['max_i_osc']}] mA)")
        print(f"  Sustained Oscillation NDR (V >= 5.5 V): {format_val_err(res['r_diff'], res['r_diff_err'], 'Ohms')} (R^2 = {res['r2_ndr']:.5f})")
        
    plot_part_a(res_5mm, res_8mm)
    return {"5mm": res_5mm, "8mm": res_8mm}


if __name__ == "__main__":
    run()
