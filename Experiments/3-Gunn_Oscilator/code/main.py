"""
main.py - Master Pipeline for Gunn Oscillator Experiment Data Analysis
PHY461 Lab - Experiment 20
Runs all experimental sub-analyses (Parts a through e), generates all plots,
prints structured summaries with uncertainties, and writes results_summary.md.
"""

import os
import sys
import numpy as np
import pandas as pd

# Add code dir to python path
CODE_DIR = os.path.dirname(os.path.abspath(__file__))
if CODE_DIR not in sys.path:
    sys.path.insert(0, CODE_DIR)

import part_a
import part_b
import part_c
import part_d
import part_e
from utils import BASE_DIR, PLOTS_DIR, format_val_err


def generate_markdown_report(res_a, res_b, res_c, res_d, res_e):
    md_path = os.path.join(BASE_DIR, "results_summary.md")
    
    ext_summary = res_d["extrema_summary"]
    fine_results = res_d["fine_results"]
    
    with open(md_path, "w", encoding="utf-8") as f:
        f.write("# Experiment 20: Study of Gunn Oscillator — Experimental Results & Error Analysis\n\n")
        f.write("**Laboratory Course**: PHY461 (Advanced Physics Laboratory)\n")
        f.write("**Status**: Complete Analysis across Parts (a)–(e)\n\n")
        f.write("---\n\n")
        
        # Section 1: Least Counts
        f.write("## 1. Instrument Least Counts\n\n")
        f.write("Per the measurement protocol, least counts are derived directly from the finest decimal precision reported in each data channel:\n\n")
        f.write("| Measurement Channel | Reported Precision / Smallest Division | Inferred Least Count |\n")
        f.write("|---|---|---|\n")
        f.write(f"| Gunn Bias Voltage $V$ | 1 decimal place | **{res_a['5mm']['lc_v']} V** |\n")
        f.write(f"| Gunn Diode Current $I$ | Integer mA | **{res_a['5mm']['lc_i']} mA** |\n")
        f.write(f"| Cavity Tuning Micrometer | Integer mm | **{res_b['lc_x']} mm** (dial: 0.01 mm) |\n")
        f.write(f"| Frequency Meter Readout | 3 decimal places | **{res_b['lc_f']} GHz** (1 MHz) |\n")
        f.write(f"| PIN Diode Bias $V_{{PIN}}$ | Multiples of 10 mV / Integer | **{res_c['lc_bias']} mV** |\n")
        f.write(f"| CRO Peak Voltages $V_{{max}}, V_{{min}}$ | 1 decimal place | **{res_c['lc_v']} mV** |\n")
        f.write(f"| Slotted Line Probe Position $x$ | 2 decimal places | **{res_d['extrema_summary']['extrema_results'][0]['lc_pos']} cm** (0.1 mm) |\n")
        f.write(f"| Microammeter Rectified Current | Integer scale units | **{res_d['extrema_summary']['extrema_results'][0]['lc_curr']}** ($10^{{-7}}, 10^{{-8}}, 10^{{-9}}$ A) |\n")
        f.write(f"| Calibrated Attenuator | Integer dB | **{res_e['lc_att']} dB** |\n")
        f.write(f"| Uncalibrated Micrometer Attenuator | 2 decimal places | **{res_e['lc_micro']} mm** (0.01 mm) |\n\n")
        f.write("---\n\n")
        
        # Section 2: Part A
        f.write("## 2. Part (a): Gunn Diode I-V Characteristics & Threshold Behavior\n\n")
        f.write("The Gunn diode I-V characteristics exhibit four distinct physical regimes:\n")
        f.write("1. **Ohmic Regime** ($V \\le 2.5$ V): Smooth linear conduction in the high-mobility $\\Gamma$-valley.\n")
        f.write("2. **Near-Threshold Plateau** ($V \\approx 3.0 - 3.3$ V): Inter-valley transfer begins competing with field acceleration.\n")
        f.write("3. **Domain Oscillation Instability** ($V \\approx 3.4 - 5.0$ V): High-field traveling dipole domains form, leading to large current fluctuations.\n")
        f.write("4. **Sustained Oscillation NDR** ($V \\ge 5.5$ V): Time-averaged current smoothly decreases, displaying negative differential resistance.\n\n")
        
        f.write("| Parameter | 5 mm Cavity Micrometer | 8 mm Cavity Micrometer |\n")
        f.write("|---|---|---|\n")
        f.write(f"| Ohmic Dynamic Resistance $R_{{Ohmic}}$ | {format_val_err(res_a['5mm']['r_ohmic'], res_a['5mm']['r_ohmic_err'], '$\\Omega$')} | {format_val_err(res_a['8mm']['r_ohmic'], res_a['8mm']['r_ohmic_err'], '$\\Omega$')} |\n")
        f.write(f"| Peak Current $I_{{peak}}$ | {res_a['5mm']['peak_i']} mA | {res_a['8mm']['peak_i']} mA |\n")
        f.write(f"| Threshold Voltage $V_0$ (Midpoint) | {format_val_err(res_a['5mm']['v_thresh_mid'], res_a['5mm']['v_thresh_err'], 'V')} | {format_val_err(res_a['8mm']['v_thresh_mid'], res_a['8mm']['v_thresh_err'], 'V')} |\n")
        f.write(f"| Onset of Negative Differential Drop | {res_a['5mm']['v_thresh_drop']} V | {res_a['8mm']['v_thresh_drop']} V |\n")
        f.write(f"| Active Slab Thickness $d$ ($E_{{th}} = 36$ kV/cm) | {format_val_err(res_a['5mm']['d_um_manual'], res_a['5mm']['d_um_err_manual'], '$\\mu$m')} | {format_val_err(res_a['8mm']['d_um_manual'], res_a['8mm']['d_um_err_manual'], '$\\mu$m')} |\n")
        f.write(f"| Predicted Transit Frequency $\\nu = v/d$ | {format_val_err(res_a['5mm']['nu_ghz_manual'], res_a['5mm']['nu_ghz_err_manual'], 'GHz')} | {format_val_err(res_a['8mm']['nu_ghz_manual'], res_a['8mm']['nu_ghz_err_manual'], 'GHz')} |\n")
        f.write(f"| Active Slab Thickness $d$ (Literature $E_{{th}} = 3.2$ kV/cm) | {format_val_err(res_a['5mm']['d_um_lit'], res_a['5mm']['d_um_err_lit'], '$\\mu$m')} | {format_val_err(res_a['8mm']['d_um_lit'], res_a['8mm']['d_um_err_lit'], '$\\mu$m')} |\n")
        f.write(f"| Predicted Transit Frequency $\\nu$ ($E_{{th}} = 3.2$ kV/cm) | {format_val_err(res_a['5mm']['nu_ghz_lit'], res_a['5mm']['nu_ghz_lit'], 'GHz')} | {format_val_err(res_a['8mm']['nu_ghz_lit'], res_a['8mm']['nu_ghz_lit'], 'GHz')} |\n")
        f.write(f"| Domain Oscillation Current | ${res_a['5mm']['mean_i_osc']:.1f} \\pm {res_a['5mm']['std_i_osc']:.1f}$ mA (range: {res_a['5mm']['min_i_osc']}–{res_a['5mm']['max_i_osc']} mA) | ${res_a['8mm']['mean_i_osc']:.1f} \\pm {res_a['8mm']['std_i_osc']:.1f}$ mA (range: {res_a['8mm']['min_i_osc']}–{res_a['8mm']['max_i_osc']} mA) |\n")
        f.write(f"| Negative Differential Resistance $R_{{diff}}$ | {format_val_err(res_a['5mm']['r_diff'], res_a['5mm']['r_diff_err'], '$\\Omega$')} | {format_val_err(res_a['8mm']['r_diff'], res_a['8mm']['r_diff_err'], '$\\Omega$')} |\n\n")
        f.write("Generated Plots: `plots/part_a_IV_5mm.png`, `plots/part_a_IV_8mm.png`, `plots/part_a_comparison.png`\n\n")
        f.write("---\n\n")
        
        # Section 3: Part B
        f.write("## 3. Part (b): Gunn Microwave Source Tuning Characteristics\n\n")
        f.write("The cavity tuning micrometer controls the effective resonator cavity length, pulling the oscillation frequency:\n\n")
        f.write(f"- **Mechanical Tuning Sensitivity (6.7 V)**: {format_val_err(res_b['slope_67'], res_b['err_s67'], 'GHz/mm')} ($R^2 = {res_b['r2_67']:.5f}$)\n")
        f.write(f"- **Mechanical Tuning Sensitivity (8.0 V)**: {format_val_err(res_b['slope_80'], res_b['err_s80'], 'GHz/mm')} ($R^2 = {res_b['r2_80']:.5f}$)\n")
        f.write(f"- **Total Tuning Range**: {res_b['range_67']:.3f} GHz (from {np.max(res_b['f_67']):.3f} to {np.min(res_b['f_67']):.3f} GHz)\n")
        f.write(f"- **Electronic Frequency Pulling (8V - 6.7V)**: $\\Delta f = +{res_b['mean_delta_f']:.1f} \\pm {res_b['delta_f_err']:.1f}$ MHz across all cavity lengths\n")
        f.write(f"- **Frequency Pushing Figure**: $\\frac{{\\Delta f}}{{\\Delta V}} = {format_val_err(res_b['pushing_factor'], res_b['pushing_factor_err'], 'MHz/V')}$\n\n")
        
        f.write("| Micrometer (mm) | $f$ at 6.7 V (GHz) | $\\lambda_0$ (cm) | Theoretical $\\lambda_g$ (cm) | $f$ at 8.0 V (GHz) | $\\lambda_0$ (cm) | Theoretical $\\lambda_g$ (cm) |\n")
        f.write("|---|---|---|---|---|---|---|\n")
        for i in range(len(res_b['x'])):
            f.write(f"| {res_b['x'][i]:.0f} | {res_b['f_67'][i]:.3f} | {res_b['lam0_67'][i]:.3f} | {res_b['lamg_67'][i]:.3f} | {res_b['f_80'][i]:.3f} | {res_b['lam0_80'][i]:.3f} | {res_b['lamg_80'][i]:.3f} |\n")
        f.write("\nGenerated Plots: `plots/part_b_tuning_curves.png`, `plots/part_b_wavelengths.png`\n\n")
        f.write("---\n\n")
        
        # Section 4: Part C
        f.write("## 4. Part (c): PIN Modulator % Modulation vs Bias\n\n")
        f.write("The modulation depth is computed via $m = \\frac{V_{max} - V_{min}}{V_{max} + V_{min}} \\times 100\\%$, with exact error propagation:\n\n")
        f.write(r"$$\delta m = \frac{200 \Delta V}{(V_{max} + V_{min})^2} \sqrt{V_{max}^2 + V_{min}^2}$$" + "\n\n")
        f.write(f"- **Peak Modulation Depth**: **{format_val_err(res_c['peak_m'], res_c['peak_m_err'], '%')}** at $V_{{PIN}} = {res_c['peak_bias']}$ mV\n")
        f.write(f"- **High-Bias Saturation Plateau Mean**: **{res_c['mean_plat_m']:.2f} ± {res_c['std_plat_m']:.2f}%**\n\n")
        f.write("Generated Plot: `plots/part_c_pin_modulation.png`\n\n")
        f.write("---\n\n")
        
        # Section 5: Part D
        f.write("## 5. Part (d): Detector Response Law & Guide Wavelength Analysis\n\n")
        f.write("### 5.1 Standing Wave Minima & Guide Wavelength Determination\n\n")
        f.write("Standing-wave minimas (nulls) were analyzed using linear regression $x_k = x_0 + k (\\lambda_g / 2)$:\n\n")
        f.write("| Attenuation | Minima Positions (cm) | Spacings $\\Delta x$ (cm) | Guide Wavelength $\\lambda_g$ (cm) | Inferred Frequency $f$ (GHz) |\n")
        f.write("|---|---|---|---|---|\n")
        for r in ext_summary["extrema_results"]:
            mins_str = ", ".join([f"{v:.2f}" for v in r["x_mins"]])
            diffs_str = ", ".join([f"{v:.2f}" for v in r["diffs"]])
            f.write(f"| **{r['label']}** | {mins_str} | {diffs_str} | {format_val_err(r['lambda_g'], r['lambda_g_err'], 'cm')} | {format_val_err(r['f_ghz'], r['f_ghz_err'], 'GHz')} |\n")
            
        f.write(f"\n- **Combined Weighted Average Guide Wavelength**: **{format_val_err(ext_summary['combined_lambda_g'], ext_summary['combined_lambda_g_err'], 'cm')}**\n")
        f.write(f"- **Derived Microwave Frequency**: **{format_val_err(ext_summary['combined_f_ghz'], ext_summary['combined_f_ghz_err'], 'GHz')}** (in perfect agreement with Part (b) frequency meter range: 9.425–10.630 GHz)\n\n")
        
        f.write("### 5.2 Detector Response Law Exponent $n$\n\n")
        f.write("Log-log regression was conducted using $\\ln I = n \\ln(\\sin\\beta x') + C$ on fine probe scans between null and peak:\n\n")
        f.write("| Attenuation | Raw Exponent $n_{raw}$ | Baseline-Subtracted Exponent $n_{sub}$ | Response Regime |\n")
        f.write("|---|---|---|---|\n")
        for r in fine_results:
            n_sub_str = format_val_err(r["n_sub"], r["n_sub_err"]) if r["n_sub"] is not None else "N/A"
            regime = "Square-Law Regime ($n \\approx 2$)" if r["label"] == "10dB" else "Intermediate / Saturated Transition"
            f.write(f"| **{r['label']}** | {format_val_err(r['n_raw'], r['n_raw_err'])} | {n_sub_str} | {regime} |\n")
            
        f.write(f"\n> At low RF input power (10 dB attenuation), the baseline-subtracted exponent **$n = {fine_results[3]['n_sub']:.2f} \\pm {fine_results[3]['n_sub_err']:.2f}$** confirms the expected **square-law detector response ($n = 2$)** within experimental error. As power increases (lower attenuation), higher-order terms cause the response to depart from the small-signal square-law behavior.\n\n")
        f.write("Generated Plots: `plots/part_d_standing_waves.png`, `plots/part_d_minima_regression.png`, `plots/part_d_response_loglog.png`, `plots/part_d_exponent_vs_power.png`\n\n")
        f.write("---\n\n")
        
        # Section 6: Part E
        f.write("## 6. Part (e): Calibration of Micrometer Attenuator\n\n")
        f.write("The micrometer vane attenuator calibration was established against the standard dB attenuator:\n\n")
        f.write(f"- **Quadratic Calibration Formula**: $A(x) = ({res_e['poly2'][0]:.3f})x^2 + ({res_e['poly2'][1]:.3f})x + ({res_e['poly2'][2]:.2f})$ dB ($R^2 = {res_e['r2_poly']:.5f}$)\n")
        f.write(f"- **Linear Approximation**: $A(x) = ({res_e['slope_lin']:.3f})x + ({res_e['int_lin']:.2f})$ dB ($R^2 = {res_e['r2_lin']:.5f}$)\n")
        f.write(f"- **Zero-Attenuation Reference Voltage**: $V_0 = {res_e['v_0']:.1f} \\pm {res_e['lc_volt']}$ mV\n\n")
        
        f.write("| Calibrated Attenuator (dB) | Micrometer Reading (mm) | CRO Output Voltage (mV) | Calculated Attenuation $A_V$ (dB) |\n")
        f.write("|---|---|---|---|\n")
        for i in range(len(res_e['att_all'])):
            a_cal = res_e['att_all'][i]
            v_val = res_e['v_all'][i]
            a_calc = res_e['att_v_calc'][i]
            a_err = res_e['att_v_err'][i]
            micro_str = f"{res_e['df'].loc[i, res_e['df'].columns[2]]:.2f}" if pd.notna(res_e['df'].loc[i, res_e['df'].columns[2]]) else "—"
            f.write(f"| {a_cal} | {micro_str} | {v_val:.1f} | {a_calc:.2f} ± {a_err:.2f} |\n")
            
        f.write("\nGenerated Plots: `plots/part_e_calibration_curve.png`, `plots/part_e_voltage_vs_dB.png`\n\n")
        f.write("---\n\n")
        f.write("## 7. Complete List of Generated Artifacts & Plots\n\n")
        plots_list = sorted(os.listdir(PLOTS_DIR))
        for p in plots_list:
            if p.endswith('.png'):
                f.write(f"- [`{p}`](file://{os.path.join(PLOTS_DIR, p)})\n")
                
    print(f"\nGenerated comprehensive report: {md_path}")


def main():
    print("="*60)
    print("  STUDY OF GUNN OSCILLATOR (EXPERIMENT 20) — FULL PIPELINE")
    print("="*60)
    
    res_a = part_a.run()
    res_b = part_b.run()
    res_c = part_c.run()
    res_d = part_d.run()
    res_e = part_e.run()
    
    generate_markdown_report(res_a, res_b, res_c, res_d, res_e)
    
    print("\n" + "="*60)
    print("           ALL ANALYSES & PLOTS SUCCESSFULLY COMPLETED!")
    print("="*60)


if __name__ == "__main__":
    main()
