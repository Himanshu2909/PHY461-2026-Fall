"""
part_c.py - Analysis of PIN Diode Modulator Characteristics (Part c)
Data:
  - c-PIN_Characteristics.csv
Calculations:
  - Least count detection for PIN bias, Vmax, and Vmin
  - Re-calculation of modulation depth m = (Vmax - Vmin) / (Vmax + Vmin) * 100%
  - Exact analytical error propagation delta m from delta Vmax and delta Vmin
  - Identification of peak modulation depth and saturation behavior
  - Publication-quality modulation depth vs bias plot
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from utils import DATA_DIR, PLOTS_DIR, detect_least_count, format_val_err


def analyze_part_c():
    file_c = os.path.join(DATA_DIR, "c-PIN_Characteristics.csv")
    df = pd.read_csv(file_c)
    
    col_bias = df.columns[0]
    col_vmax = df.columns[1]
    col_vmin = df.columns[2]
    col_m_raw = df.columns[3]
    
    lc_bias = detect_least_count(df[col_bias])
    lc_vmax = detect_least_count(df[col_vmax])
    lc_vmin = detect_least_count(df[col_vmin])
    lc_v = max(lc_vmax, lc_vmin)
    
    bias = df[col_bias].values
    vmax = df[col_vmax].values
    vmin = df[col_vmin].values
    m_raw = df[col_m_raw].values
    
    # Re-calculate modulation depth
    m_calc = ((vmax - vmin) / (vmax + vmin)) * 100.0
    
    # Exact analytical error propagation
    # dm/dVmax = 200 * Vmin / (Vmax + Vmin)^2
    # dm/dVmin = -200 * Vmax / (Vmax + Vmin)^2
    denom = (vmax + vmin) ** 2
    dm_dvmax = (200.0 * vmin) / denom
    dm_dvmin = (-200.0 * vmax) / denom
    m_err = np.sqrt((dm_dvmax * lc_v) ** 2 + (dm_dvmin * lc_v) ** 2)
    
    # Peak modulation depth
    peak_idx = np.argmax(m_calc)
    peak_m = m_calc[peak_idx]
    peak_m_err = m_err[peak_idx]
    peak_bias = bias[peak_idx]
    
    # High-bias plateau statistics (bias >= 1200 mV)
    plat_mask = bias >= 1200
    plat_m = m_calc[plat_mask]
    plat_err = m_err[plat_mask]
    mean_plat_m = np.mean(plat_m)
    std_plat_m = np.std(plat_m, ddof=1)
    
    results = {
        "lc_bias": lc_bias,
        "lc_v": lc_v,
        "bias": bias,
        "vmax": vmax,
        "vmin": vmin,
        "m_raw": m_raw,
        "m_calc": m_calc,
        "m_err": m_err,
        "peak_m": peak_m,
        "peak_m_err": peak_m_err,
        "peak_bias": peak_bias,
        "mean_plat_m": mean_plat_m,
        "std_plat_m": std_plat_m,
        "df": df
    }
    return results


def plot_part_c(res):
    bias = res["bias"]
    m_calc = res["m_calc"]
    m_err = res["m_err"]
    
    fig, ax = plt.subplots(figsize=(8, 5.5))
    
    # Plot experimental data with error bars
    ax.errorbar(bias, m_calc, yerr=m_err, xerr=res["lc_bias"], fmt='o-', color='#1f77b4',
                label='Calculated % Modulation $m$', capsize=4, markersize=6, zorder=4)
    
    # Highlight peak
    ax.scatter([res["peak_bias"]], [res["peak_m"]], color='#d62728', marker='*', s=180,
               zorder=5, label=f'Peak Modulation: {res["peak_m"]:.2f}% at {res["peak_bias"]} mV')
    
    # Plateau reference line
    ax.axhline(res["mean_plat_m"], color='#2ca02c', linestyle='--', alpha=0.8,
               label=f'Plateau Average ($V_{{PIN}} \\geq 1.2$ V): {res["mean_plat_m"]:.2f}%')
    
    ax.set_xlabel('PIN Diode Bias Voltage $V_{PIN}$ (mV)')
    ax.set_ylabel('Modulation Depth $m$ (%)')
    ax.set_title('PIN Modulator Characteristics (% Modulation vs Bias)')
    ax.set_xlim(150, 3050)
    ax.set_ylim(-0.5, 16.5)
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend(loc='upper right', frameon=True, fontsize=9.5)
    
    # Info box positioned in the clear lower-right quadrant
    info_text = (
        f"Peak: $m_{{max}} = {res['peak_m']:.2f} \\pm {res['peak_m_err']:.2f}\\%$\n"
        f"At Bias: ${res['peak_bias']}$ mV\n"
        f"Plateau Mean: ${res['mean_plat_m']:.2f} \\pm {res['std_plat_m']:.2f}\\%$\n"
        f"Instrument Least Counts:\n"
        f"  $\\Delta V_{{PIN}} = {res['lc_bias']}$ mV\n"
        f"  $\\Delta V_{{CRO}} = {res['lc_v']}$ mV"
    )
    ax.text(0.58, 0.08, info_text, transform=ax.transAxes, verticalalignment='bottom',
            fontsize=9.5, bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9, edgecolor='gray'))
    
    save_path = os.path.join(PLOTS_DIR, "part_c_pin_modulation.png")
    fig.savefig(save_path, dpi=300)
    plt.close(fig)
    print(f"Saved: {save_path}")


def run():
    print("\n" + "="*50)
    print("        PART (c): PIN MODULATOR CHARACTERISTICS")
    print("="*50)
    
    res = analyze_part_c()
    
    print(f"  Least Counts: PIN Bias = {res['lc_bias']} mV, CRO Voltages = {res['lc_v']} mV")
    print(f"  Peak Modulation Depth: {format_val_err(res['peak_m'], res['peak_m_err'], '%')} at {res['peak_bias']} mV")
    print(f"  High-Bias Plateau Mean Modulation: {res['mean_plat_m']:.2f} ± {res['std_plat_m']:.2f} %")
    
    # Print sample comparison table of raw vs recomputed
    print("\n  Sample Verification Table:")
    print(f"  {'Bias (mV)':<10} {'Vmax (mV)':<10} {'Vmin (mV)':<10} {'Raw m(%)':<10} {'Calc m(%)':<12} {'Error(%)':<10}")
    for i in range(len(res["bias"])):
        print(f"  {res['bias'][i]:<10.0f} {res['vmax'][i]:<10.1f} {res['vmin'][i]:<10.1f} {res['m_raw'][i]:<10.2f} {res['m_calc'][i]:<12.2f} {res['m_err'][i]:<10.2f}")
        
    plot_part_c(res)
    return res


if __name__ == "__main__":
    run()
