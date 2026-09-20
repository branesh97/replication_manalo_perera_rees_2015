Replication files – Macroeconometrics assignment ------------------------------------
Name: Branesh Prakash
Student ID: 35775246
-------------------------------------------------------------------------------------
1. Overview

These files reproduce the empirical results for replication and extension of: Manalo, J., Perera, D. and Rees, D. (2015), "Exchange rate movements and the Australian economy".

The analysis estimates a small-open-economy structural VAR for Australia. The baseline replication covers 1985Q1–2013Q2. Extensions consider 1985Q1–2019Q4 and 1985Q1–2026Q2. Additional specifications include pandemic controls and the VIX as a measure of global financial risk.

All MATLAB scripts use paths relative to the location of the script.
-------------------------------------------------------------------------------------
2. Software

The analysis was conducted in MATLAB. No additional MATLAB toolboxes are required beyond functions used in the supplied code.
-------------------------------------------------------------------------------------
3. Folder structure


    README.md
    baseline.m
    ext_sample_precovid.m
    ext_sample_full.m
    ext_sample_pandemic_controls.m
    ext_vix_original_sample.m
    ext_vix_precovid.m
    ext_vix_pandemic_controls.m
    make_main_irf_figure.m
    online_links.txt

    _func
        olsvar_soe.m
        calculate_IRF_FEVD.m
        bootstrap_after_bootstrap_soe.m

    data
        data_raw.xlsx
        data_sources.txt

    results
        baseline.mat
        ext_sample_1.precovid.mat
        ext_sample_2.full.mat
        ext_sample_3.pandemic_controls.mat
        results.txt

    figures
        [generated figure files]
-------------------------------------------------------------------------------------
3. Data

The MATLAB scripts all read: data/data_raw.xlsx, Sheet: Raw, Range: B2:H167
The observations are quarterly from 1985Q1 to 2026Q2.
Columns are:
1. US real GDP
2. Australian terms of trade
3. Australian real GDP
4. Australian trimmed mean inflation
5. Australian cash rate
6. Australian real trade-weighted exchange rate (RTWI)
7. CBOE Volatility Index (VIX)
VIX data are available from 1990Q1, so all VIX specifications use matched samples beginning in 1990Q1.

Detailed data sources, series identifiers, download links, and any processing performed before assembling data_raw.xlsx are documented in data_sources.txt.
-------------------------------------------------------------------------------------
4. Model specification

Specified in paper and each Matlab code file.

Additional notes: Confidence intervals use a Kilian-style bootstrap-after-bootstrap procedure adapted to the small-open-economy VAR. The first bootstrap stage estimates finite-sample bias in the VAR coefficients. The estimated bias is used to construct a bias-corrected VAR. If necessary, the bias correction is shrunk toward the original
estimate until the generating VAR is stable. The second bootstrap stage generates the distribution used for the reported impulse-response confidence intervals. Both stages use 1,000 bootstrap replications. A fixed random-number seed of 12345 is used for reproducibility. The small-open-economy restrictions and deterministic controls are
imposed in every bootstrap re-estimation. The reported confidence intervals are 85 percent and 95 percent percentile intervals.
-------------------------------------------------------------------------------------
5. Replication instructions
Run the scripts from the project directory in the following order:
1. baseline.m
2. ext_sample_precovid.m
3. ext_sample_full.m
4. ext_sample_pandemic_controls.m
The first four scripts generate the baseline and sample-extension results and save the corresponding MAT files in /results.
Then, if replicating the VIX extension, run:
5. ext_vix_original_sample.m
6. ext_vix_precovid.m
7. ext_vix_pandemic_controls.m
Finally run:
8. make_main_irf_figure.m
This reads the saved baseline, pre-COVID and preferred full-sample results and exports the combined main impulse-response figure to the figures directory.
-------------------------------------------------------------------------------------
6. Main script descriptions

baseline.m: Replicates the baseline 1985Q1–2013Q2 VAR, impulse responses and forecast-error variance decompositions, and compares the FEVD results with the published values.

ext_sample_precovid.m: Extends the baseline specification through 2019Q4 and computes bootstrap confidence intervals.

ext_sample_full.m: Extends the baseline specification through 2026Q2 without pandemic controls. This specification is used primarily to diagnose the influence of the pandemic period and unusual residual observations.

ext_sample_pandemic_controls.m: Estimates the preferred 1985Q1–2026Q2 extension with 2020Q2 and 2021Q3 pandemic pulse controls and bootstrap confidence intervals.

ext_vix_original_sample.m: Compares models with and without VIX over the matched 1990Q1–2013Q2 sample.

ext_vix_precovid.m: Compares models with and without VIX over 1990Q1–2019Q4 and computes bootstrap confidence intervals.

ext_vix_pandemic_controls.m: Compares models with and without VIX over 1990Q1–2026Q2 using the pandemic controls and computes bootstrap confidence intervals.

make_main_irf_figure.m: Constructs the combined impulse-response figure from saved result files, used in report.
-------------------------------------------------------------------------------------
7. Output files
The /results directory contains MAT files produced by the main scripts. results.txt is included as a convenience record of the numerical output obtained from the submitted code. It is not used as an input by any estimation script. The /figures directory contains figure files generated from the MATLAB results. Because the random-number seed is fixed, rerunning the bootstrap scripts should reproduce the submitted confidence intervals.