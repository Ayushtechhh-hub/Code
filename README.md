# Code
3.5 Estimation of Model Parameters Using Linear Regression of WVTR Data
WVTR-Based Diffusion Coefficient Estimation for Multilayer Coatings

Overview
This document presents a Python-based analytical framework for estimating effective diffusion coefficients of biopolymer and biopolymer –PLA multilayer coating systems using linear regression of experimentally measured WVTR data. The methodology is suitable for thesis-level research and continuum-scale transport modelling.

Objective
The objectives are to quantify moisture diffusion coefficients for biopolymer coatings at different coating builds, extract the permeable diffusion coefficient of PLA, compare results with an impermeable PLA limit, and visualize WVTR–thickness trends.

Physical Constants
Mw = 18.01528 g/mol (molecular weight of water)
Ci = 1.1 mol/m³ (inside concentration)
Co = 0.55 mol/m³ (outside concentration)
Hp = 62.5 µm (paper thickness)
D_PLA,imp = 1×10⁻¹⁵ m²/s (impermeable PLA)

Governing Equations

Linear Regression:
WVTR = mH + b

Diffusion Coefficient:
D = Mw (Ci − Co) / |m|

Impermeable PLA Limit:
WVTR = Mw (Ci − Co) / (Hp/Dp + Hc/Dc)

Where,  Dc is the Diffusion  coefficient  of coating
        Hc is the  Thickness of coating
Outputs
The script generates individual regression plots for each coating regime and a combined comparison plot. Extracted diffusion coefficients are printed in results-ready format.

Scientific Relevance
The framework provides experimentally derived transport parameters suitable for COMSOL simulations and multilayer barrier analysis.

Execution
Run the script using:
python wvtr_diffusivity_analysis.py

