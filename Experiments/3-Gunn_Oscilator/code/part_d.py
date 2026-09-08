"""
part_d.py - Response Law of Detector Diode & Guide Wavelength Analysis (Part d)
Data:
  - Extremas: 0dB, 3dB, 6dB, 10dB
  - Fine scans: 0dB, 3dB, 6dB, 10dB
Calculations:
  - Least count detection for probe position and rectified current
  - Identification of standing-wave minimas (nulls)
  - Linear regression x_k = x_0 + k * (lambda_g / 2) to determine lambda_g +/- delta lambda_g
  - Successive differences statistics Delta x_min
  - Weighted average lambda_g across all 4 attenuations
  - Calculated microwave frequency f_guide +/- delta f_guide
  - Fine scan log-log regression: ln(I) = n * ln(sin beta x') + C to determine exponent n
  - Analysis of square-law response (n ~ 2 at low power) vs departure at high power
  - Publication-quality plots: standing waves, minima regression, log-log response, n vs attenuation
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from utils import DATA_DIR, PLOTS_DIR, detect_least_count, linear_fit_with_errors, weighted_average, format_val_err

# Physical constants
C_CMS = 2.99792458e10   # Speed of light in cm/s
A_WAVEGUIDE = 2.286     # WR-90 waveguide broad dimension in cm
LAMBDA_C = 2.0 * A_WAVEGUIDE  # Cutoff wavelength = 4.572 cm

EXTREMA_CONFIGS = [
    {"label": "0dB",  "file": "d-i-response-0dB-extremas.csv",   "scale": 1e-7, "color": "#1f77b4"},
    {"label": "3dB",  "file": "d-ii-response-3dB-extremas.csv",  "scale": 1e-8, "color": "#d62728"},
    {"label": "6dB",  "file": "d-iii-response-6dB-extremas.csv", "scale": 1e-9, "color": "#2ca02c"},
    {"label": "10dB", "file": "d-iv-response-10dB-extremas.csv", "scale": 1e-9, "color": "#9467bd"}
]

FINE_CONFIGS = [
    {"label": "0dB",  "file": "d-i-response-0dB-fine.csv",   "scale": 1e-7, "x_ref": 8.25},
    {"label": "3dB",  "file": "d-ii-response-3dB-fine.csv",  "scale": 1e-8, "x_ref": 8.30},
    {"label": "6dB",  "file": "d-iii-response-6dB-fine.csv", "scale": 1e-9, "x_ref": 8.19},
    {"label": "10dB", "file": "d-iv-response-10dB-fine.csv", "scale": 1e-9, "x_ref": 8.40}
]


def extract_minimas(df, label):
    """
    Isolates the true periodic standing wave minimas (nulls) from the extrema file.
    In 0dB, the experimenter recorded intermediate dips (abs_I ~ 25-33); the true nulls
    are at abs_I <= 20 (8.25, 10.20, 12.00, 13.80, 15.60 cm).
    In 3dB, 6dB, 10dB, the recorded tables consist of alternating minima and maxima,
    so even rows (indices 0, 2, 4, 6, 8) are the true standing-wave minima.
    """
    pos_col = df.columns[1]
    curr_col = df.columns[0]
    df_sorted = df.sort_values(by=pos_col).reset_index(drop=True)
    abs_i = df_sorted[curr_col].abs().values
    
    if label == "0dB":
        min_mask = abs_i <= 20
        minima_df = df_sorted[min_mask].copy()
    else:
        minima_df = df_sorted.iloc[::2].copy()
        
    return minima_df, df_sorted


def analyze_extrema():
    extrema_results = []
    
    for cfg in EXTREMA_CONFIGS:
        label = cfg["label"]
        file_path = os.path.join(DATA_DIR, cfg["file"])
        df = pd.read_csv(file_path)
        
        lc_curr = detect_least_count(df[df.columns[0]])
        lc_pos = detect_least_count(df[df.columns[1]])
        
        minima_df, full_df = extract_minimas(df, label)
        pos_col = df.columns[1]
        x_mins = minima_df[pos_col].values
        k = np.arange(len(x_mins))
        
        # Linear regression: x_k = x_0 + k * (lambda_g / 2)
        slope, slope_err, intercept, intercept_err, r2, _ = linear_fit_with_errors(
            k, x_mins, y_err=np.full_like(x_mins, lc_pos)
        )
        lambda_g = 2.0 * slope
        lambda_g_err = 2.0 * slope_err
        
        # Successive differences
        diffs = np.diff(x_mins)
        mean_diff = np.mean(diffs)
        err_diff = np.std(diffs, ddof=1) / np.sqrt(len(diffs)) if len(diffs) > 1 else lc_pos
        lambda_g_diff = 2.0 * mean_diff
        lambda_g_diff_err = 2.0 * err_diff
        
        # Calculate microwave frequency: f = c * sqrt(1/lambda_g^2 + 1/lambda_c^2)
        inv_lamg2 = 1.0 / (lambda_g ** 2)
        inv_lamc2 = 1.0 / (LAMBDA_C ** 2)
        f_hz = C_CMS * np.sqrt(inv_lamg2 + inv_lamc2)
        f_ghz = f_hz / 1.0e9
        
        # df/dlambda_g = - c / (lambda_g^3 * sqrt(1/lambda_g^2 + 1/lambda_c^2))
        df_dlamg = (C_CMS / (lambda_g ** 3)) / np.sqrt(inv_lamg2 + inv_lamc2)
        f_ghz_err = (df_dlamg * lambda_g_err) / 1.0e9
        
        res = {
            "label": label,
            "lc_curr": lc_curr,
            "lc_pos": lc_pos,
            "scale": cfg["scale"],
            "full_df": full_df,
            "minima_df": minima_df,
            "x_mins": x_mins,
            "k": k,
            "diffs": diffs,
            "slope": slope,
            "slope_err": slope_err,
            "r2": r2,
            "lambda_g": lambda_g,
            "lambda_g_err": lambda_g_err,
            "lambda_g_diff": lambda_g_diff,
            "lambda_g_diff_err": lambda_g_diff_err,
            "f_ghz": f_ghz,
            "f_ghz_err": f_ghz_err,
            "color": cfg["color"]
        }
        extrema_results.append(res)
        
    # Combined weighted average lambda_g
    lam_vals = [r["lambda_g"] for r in extrema_results]
    lam_errs = [r["lambda_g_err"] for r in extrema_results]
    comb_lam_g, comb_lam_g_err = weighted_average(lam_vals, lam_errs)
    
    # Combined frequency
    inv_lamg2 = 1.0 / (comb_lam_g ** 2)
    inv_lamc2 = 1.0 / (LAMBDA_C ** 2)
    comb_f_hz = C_CMS * np.sqrt(inv_lamg2 + inv_lamc2)
    comb_f_ghz = comb_f_hz / 1.0e9
    df_dlamg = (C_CMS / (comb_lam_g ** 3)) / np.sqrt(inv_lamg2 + inv_lamc2)
    comb_f_ghz_err = (df_dlamg * comb_lam_g_err) / 1.0e9
    
    summary = {
        "extrema_results": extrema_results,
        "combined_lambda_g": comb_lam_g,
        "combined_lambda_g_err": comb_lam_g_err,
        "combined_f_ghz": comb_f_ghz,
        "combined_f_ghz_err": comb_f_ghz_err
    }
    return summary


def analyze_fine(extrema_summary):
    ext_dict = {r["label"]: r for r in extrema_summary["extrema_results"]}
    fine_results = []
    
    for cfg in FINE_CONFIGS:
        label = cfg["label"]
        file_path = os.path.join(DATA_DIR, cfg["file"])
        df = pd.read_csv(file_path)
        
        pos_col = df.columns[1]
        curr_col = df.columns[0]
        lc_curr = detect_least_count(df[curr_col])
        lc_pos = detect_least_count(df[pos_col])
        
        df["abs_I"] = df[curr_col].abs()
        df = df.sort_values(by=pos_col).reset_index(drop=True)
        
        # Use lambda_g from this attenuation
        lam_g = ext_dict[label]["lambda_g"]
        lam_g_err = ext_dict[label]["lambda_g_err"]
        beta = 2.0 * np.pi / lam_g
        
        # Reference minimum
        x_ref = cfg["x_ref"]
        df["x_prime"] = (df[pos_col] - x_ref).abs()
        df["sin_bx"] = np.sin(beta * df["x_prime"])
        
        # Case A: Raw log-log fit: ln(I) = n * ln(sin beta x') + C
        mask_raw = (df["sin_bx"] > 0.05) & (df["abs_I"] > 0)
        log_sin_raw = np.log(df.loc[mask_raw, "sin_bx"].values)
        log_i_raw = np.log(df.loc[mask_raw, "abs_I"].values)
        
        slope_raw, err_raw, int_raw, _, r2_raw, fit_raw = linear_fit_with_errors(
            log_sin_raw, log_i_raw
        )
        
        # Case B: Baseline-subtracted log-log fit: ln(I - I_min) = n * ln(sin beta x') + C
        i_floor = df["abs_I"].min()
        mask_sub = (df["sin_bx"] > 0.05) & (df["abs_I"] > i_floor)
        
        if np.sum(mask_sub) >= 3:
            log_sin_sub = np.log(df.loc[mask_sub, "sin_bx"].values)
            log_i_sub = np.log(df.loc[mask_sub, "abs_I"].values - i_floor)
            slope_sub, err_sub, int_sub, _, r2_sub, fit_sub = linear_fit_with_errors(
                log_sin_sub, log_i_sub
            )
        else:
            slope_sub, err_sub, int_sub, r2_sub, fit_sub = None, None, None, None, None
            log_sin_sub, log_i_sub = None, None
            
        res = {
            "label": label,
            "df": df,
            "lc_curr": lc_curr,
            "lc_pos": lc_pos,
            "lam_g": lam_g,
            "x_ref": x_ref,
            "i_floor": i_floor,
            "mask_raw": mask_raw,
            "log_sin_raw": log_sin_raw,
            "log_i_raw": log_i_raw,
            "n_raw": slope_raw,
            "n_raw_err": err_raw,
            "r2_raw": r2_raw,
            "fit_raw": fit_raw,
            "mask_sub": mask_sub,
            "log_sin_sub": log_sin_sub,
            "log_i_sub": log_i_sub,
            "n_sub": slope_sub,
            "n_sub_err": err_sub,
            "r2_sub": r2_sub,
            "fit_sub": fit_sub
        }
        fine_results.append(res)
        
    return fine_results


def plot_part_d(ext_summary, fine_results):
    ext_res = ext_summary["extrema_results"]
    
    # Plot 1: 4-panel Standing Waves
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.5), sharex=True)
    axes = axes.flatten()
    
    legend_handles = []
    legend_labels = []
    
    for i, res in enumerate(ext_res):
        ax = axes[i]
        lbl = res["label"]
        df = res["full_df"]
        min_df = res["minima_df"]
        pos_col = df.columns[1]
        curr_col = df.columns[0]
        scale_label = curr_col.split('(')[1].split(')')[0]
        
        max_val = np.max(df[curr_col].abs())
        ax.set_ylim(0, max_val * 1.35)
        
        line1, = ax.plot(df[pos_col], df[curr_col].abs(), 'o-', color=res["color"], alpha=0.8,
                         label='Extrema sweep')
        scat1 = ax.scatter(min_df[pos_col], min_df[curr_col].abs(), color='black', marker='v', s=80,
                           zorder=5, label='Identified Minimas (Nulls)')
        
        if i == 0:
            legend_handles = [line1, scat1]
            legend_labels = ['Extrema Sweep Points', 'Identified Standing Wave Minimas (Nulls)']
        
        ax.set_ylabel(f'Rectified Current ({scale_label})')
        ax.set_title(f'Standing Wave Profile — Attenuation: {lbl}', fontsize=12, pad=8)
        ax.grid(True, linestyle='--', alpha=0.5)
        
        # Info box in the top margin (guaranteed empty due to 35% headroom)
        info = f"$\\lambda_g = {res['lambda_g']:.3f} \\pm {res['lambda_g_err']:.3f}$ cm   ($R^2 = {res['r2']:.5f}$)"
        ax.text(0.98, 0.94, info, transform=ax.transAxes, verticalalignment='top',
                horizontalalignment='right', fontsize=10,
                bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=0.9, edgecolor='gray'))
        
    axes[2].set_xlabel('Probe Position $x$ (cm)')
    axes[3].set_xlabel('Probe Position $x$ (cm)')
    
    # Unified figure legend at the top
    fig.legend(legend_handles, legend_labels, loc='upper center', bbox_to_anchor=(0.5, 0.98),
               ncol=2, frameon=True, fontsize=11)
    fig.suptitle('Standing Wave Extrema Profiles and Identified Minima across Attenuations',
                 fontsize=15, y=1.02)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    
    save_path1 = os.path.join(PLOTS_DIR, "part_d_standing_waves.png")
    fig.savefig(save_path1, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {save_path1}")
    
    # Plot 2: Minima Linear Regression
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    for res in ext_res:
        k = res["k"]
        x_mins = res["x_mins"]
        lbl = res["label"]
        ax.errorbar(k, x_mins, yerr=res["lc_pos"], fmt='o', color=res["color"],
                    label=f'{lbl} ($\\lambda_g = {res["lambda_g"]:.3f}$ cm)', capsize=4)
        k_fit = np.linspace(0, len(k)-1, 50)
        ax.plot(k_fit, res["slope"] * k_fit + (x_mins[0] - res["slope"]*0), '--', color=res["color"], alpha=0.7)
        
    ax.set_xlabel('Minima Order Index $k$')
    ax.set_ylabel('Minima Position $x_k$ (cm)')
    ax.set_title('Guide Wavelength Determination: Minima Position vs Fringe Index $k$')
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend(loc='lower right', frameon=True)
    
    comb_text = (
        f"Combined Weighted Mean $\\lambda_g$:\n"
        f"  $\\lambda_g = {ext_summary['combined_lambda_g']:.3f} \\pm {ext_summary['combined_lambda_g_err']:.3f}$ cm\n"
        f"Derived Frequency:\n"
        f"  $f_{{calc}} = {ext_summary['combined_f_ghz']:.3f} \\pm {ext_summary['combined_f_ghz_err']:.3f}$ GHz"
    )
    ax.text(0.05, 0.70, comb_text, transform=ax.transAxes,
            bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9, edgecolor='gray'))
    
    save_path2 = os.path.join(PLOTS_DIR, "part_d_minima_regression.png")
    fig.savefig(save_path2, dpi=300)
    plt.close(fig)
    print(f"Saved: {save_path2}")
    
    # Plot 3: 4-panel Log-Log Detector Response Curves
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    axes = axes.flatten()
    
    for i, res in enumerate(fine_results):
        ax = axes[i]
        lbl = res["label"]
        
        # Raw log-log
        ax.plot(res["log_sin_raw"], res["log_i_raw"], 'o', color='#1f77b4',
                label=f'Raw: $n = {res["n_raw"]:.3f} \\pm {res["n_raw_err"]:.3f}$')
        sin_grid = np.linspace(np.min(res["log_sin_raw"]), np.max(res["log_sin_raw"]), 50)
        ax.plot(sin_grid, res["fit_raw"](sin_grid), '-', color='#1f77b4', alpha=0.8)
        
        # Subtracted log-log if available
        if res["n_sub"] is not None:
            ax.plot(res["log_sin_sub"], res["log_i_sub"], 's', color='#d62728',
                    label=f'Baseline-sub: $n = {res["n_sub"]:.3f} \\pm {res["n_sub_err"]:.3f}$')
            sin_sub_grid = np.linspace(np.min(res["log_sin_sub"]), np.max(res["log_sin_sub"]), 50)
            ax.plot(sin_sub_grid, res["fit_sub"](sin_sub_grid), '--', color='#d62728', alpha=0.8)
            
        ax.set_xlabel(r'$\ln [\sin(\beta x^\prime)]$')
        ax.set_ylabel(r'$\ln (I)$ or $\ln (I - I_{min})$')
        ax.set_title(f'Detector Response Law: {lbl}')
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.legend(loc='lower right', frameon=True)
        
    fig.suptitle(r'Log-Log Detector Response Characterization ($\ln I = n \ln\sin\beta x^\prime + C$)', fontsize=15)
    save_path3 = os.path.join(PLOTS_DIR, "part_d_response_loglog.png")
    fig.savefig(save_path3, dpi=300)
    plt.close(fig)
    print(f"Saved: {save_path3}")
    
    # Plot 4: Exponent n vs Attenuation (Power Level)
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    att_vals = [0, 3, 6, 10]
    n_raw_vals = [r["n_raw"] for r in fine_results]
    n_raw_errs = [r["n_raw_err"] for r in fine_results]
    
    n_sub_vals = [r["n_sub"] for r in fine_results]
    n_sub_errs = [r["n_sub_err"] for r in fine_results]
    
    ax.errorbar(att_vals, n_sub_vals, yerr=n_sub_errs, fmt='s-', color='#d62728',
                label='Baseline-Subtracted Exponent $n_{sub}$', capsize=4, markersize=6)
    ax.errorbar(att_vals, n_raw_vals, yerr=n_raw_errs, fmt='o--', color='#1f77b4',
                label='Raw Current Exponent $n_{raw}$', capsize=4, markersize=6)
    
    ax.axhline(2.0, color='black', linestyle=':', label='Theoretical Square Law ($n = 2$)')
    ax.axhline(1.0, color='gray', linestyle='-.', label='Linear Detector Limit ($n = 1$)')
    
    ax.set_xlim(-0.8, 11.2)
    ax.set_ylim(-0.3, 5.8)
    ax.set_xlabel('RF Attenuation (dB) — Decreasing RF Power to the Right')
    ax.set_ylabel('Detector Response Exponent $n$')
    ax.set_title('Detector Response Exponent $n$ vs RF Power / Attenuation')
    ax.grid(True, linestyle='--', alpha=0.5)
    
    # Place legend at upper-left where y in [2.8, 5.5] and x in [0.2, 4.0] is empty except for (0, 4.6)
    ax.legend(loc='upper left', frameon=True, fontsize=9.5)
    
    # Place info box in the completely clear upper-right region (x in [6.5, 10.5], y in [3.0, 5.5])
    note_text = (
        "Square-law regime ($n \\approx 2$):\n"
        f"  At 10 dB: $n = {n_sub_vals[3]:.2f} \\pm {n_sub_errs[3]:.2f}$\n"
        "High power departure:\n"
        "  As power increases (lower attenuation),\n"
        "  detector departs from ideal square-law."
    )
    ax.text(0.97, 0.95, note_text, transform=ax.transAxes, verticalalignment='top',
            horizontalalignment='right', fontsize=9.5,
            bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9, edgecolor='gray'))
    
    save_path4 = os.path.join(PLOTS_DIR, "part_d_exponent_vs_power.png")
    fig.savefig(save_path4, dpi=300)
    plt.close(fig)
    print(f"Saved: {save_path4}")


def run():
    print("\n" + "="*50)
    print("  PART (d): DETECTOR RESPONSE LAW & GUIDE WAVELENGTH")
    print("="*50)
    
    ext_summary = analyze_extrema()
    ext_res = ext_summary["extrema_results"]
    
    print("\n--- Standing Wave Minima & Guide Wavelength Analysis ---")
    for r in ext_res:
        lbl = r["label"]
        print(f"\n  Attenuation: {lbl}")
        print(f"    Least Counts: Position = {r['lc_pos']} cm, Current = {r['lc_curr']} ({r['scale']} A)")
        print(f"    Minima Positions (cm): {r['x_mins']}")
        print(f"    Successive Spacing Delta x (cm): {r['diffs']}")
        print(f"    From Linear Fit: lambda_g = {format_val_err(r['lambda_g'], r['lambda_g_err'], 'cm')} (R^2 = {r['r2']:.5f})")
        print(f"    From Consecutive Diffs: lambda_g = {format_val_err(r['lambda_g_diff'], r['lambda_g_diff_err'], 'cm')}")
        print(f"    Calculated Frequency: {format_val_err(r['f_ghz'], r['f_ghz_err'], 'GHz')}")
        
    print("\n--- Combined Guide Wavelength & Frequency ---")
    print(f"  Weighted Average lambda_g: {format_val_err(ext_summary['combined_lambda_g'], ext_summary['combined_lambda_g_err'], 'cm')}")
    print(f"  Inferred Microwave Frequency: {format_val_err(ext_summary['combined_f_ghz'], ext_summary['combined_f_ghz_err'], 'GHz')}")
    
    fine_res = analyze_fine(ext_summary)
    
    print("\n--- Detector Response Exponent n from Fine Scans ---")
    for r in fine_res:
        lbl = r["label"]
        print(f"\n  Attenuation: {lbl}")
        print(f"    Raw Exponent: n = {format_val_err(r['n_raw'], r['n_raw_err'])} (R^2 = {r['r2_raw']:.5f})")
        if r["n_sub"] is not None:
            print(f"    Baseline-Subtracted Exponent: n = {format_val_err(r['n_sub'], r['n_sub_err'])} (R^2 = {r['r2_sub']:.5f})")
            
    plot_part_d(ext_summary, fine_res)
    return {"extrema_summary": ext_summary, "fine_results": fine_res}


if __name__ == "__main__":
    run()
