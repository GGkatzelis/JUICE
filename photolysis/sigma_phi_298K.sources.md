# Sources and validation for `sigma_phi_298K.csv`

| | |
|---|---|
| Table | `photolysis/sigma_phi_298K.csv`: sigma(lambda, 298 K) x phi(lambda, 298 K, 1 atm) in cm^2, 1-nm grid 200.0-750.0 nm inclusive (551 rows), one column per MCM v3.3.1 J (numeric order J1 ... J61). Empty cells = no data (J15, J16, J17, J61). `0` = outside the source's wavelength coverage or a genuine zero of the data. |
| Built by | `photolysis/build_sigma_phi.py` (importable without side effects) |
| Rebuild | `python -I photolysis/build_sigma_phi.py <path-to-tuv-x-checkout> --validate` with the project's Python environment (needs numpy, netCDF4; runs in ~3 s; `--validate` prints the spot checks and the bracket table below) |
| Data source | NCAR **tuv-x** (https://github.com/NCAR/tuv-x), commit `097b54cceb8edaaac67ffea520b43233728ecb69` (2026-10-19), licensed **Apache-2.0** (Copyright University Corporation for Atmospheric Research). The netCDF data files are TUV 5.4 data sets converted by "TUV Data Xformer 1.0"; the original TUV file name is quoted for each entry below. |
| Built on | 2026-11-02 |
| Conditions | T = 298 K, p = 1013.25 hPa, M = 2.46e19 cm^-3 for every temperature- or pressure-dependent parametrisation. The tuv-x `tint` rule is applied to multi-temperature tables (clamp T to the tabulated range, linear interpolation between rows). |
| Gridding | Linear interpolation of each source onto the 1-nm grid, zero outside the source coverage unless tuv-x's `photolysis_rate_constants.json` configures a `lower/upper extrapolation` (noted per J). Exception: the 0.01-nm O3 data (and the solar spectrum used in the test) are averaged over each 1-nm bin [lambda-0.5, lambda+0.5) instead of point-sampled. |
| MCM reference | https://cri.york.ac.uk/MCM/rates/photolysis (reaction/products, l,m,n, MCM data sources -> `photolysis/mcm_photolysis.csv`). The page header says v3.3.1 but the download/J21 note still refer to v3.2; the numbers are identical to `mcm_converter.py` lines 88-113. |

The intended use is J_i / J_NO2 = integral(sigma_phi_i F) / integral(sigma_phi_J4 F) for a lamp spectrum F(lambda) (class-A / "bottom-up" method of Carter, SAPRC and F0AM). Columns are *channel* values: J1 + J2 = sigma(O3), J31 + J32 + J33 = sigma(CHOCHO) x phi_total, etc.

## Per-J provenance

Abbreviations: xs = cross section, qy = quantum yield, cov = wavelength coverage of the source (zero outside unless stated), `tuv-x:` = reaction name in `data/photolysis_rate_constants.json`.

| J | MCM reaction | tuv-x reaction | sigma source (original TUV file / reference) | phi source | cov / caveats |
|---|---|---|---|---|---|
| J1 | O3 -> O(1D) + O2 | `O3+hv->O2+O(1D)` | `data/cross_sections/O3_1.nc` = DATAE1NEW/O3/1995Malicet_1998Brion_O3_295.txt (Malicet et al. 1995, Brion et al. 1998; 295 K, 0.01 nm), bin-averaged | `src/quantum_yields/o3-o2_o1d.F90`: Matsumi et al. (2002) JPL expression at 298 K: 0.90 for lambda <= 305, three-term formula 305-328, 0.08 for 328-340, 0 above 340 nm | xs cov 195-830 nm (full Chappuis band). MCM used Molina & Molina 1986 xs; T dependence of xs (295 vs 298 K) ignored |
| J2 | O3 -> O(3P) + O2 | `O3+hv->O2+O(3P)` | as J1 | 1 - phi(O1D) (`o3-o2_o3p.F90`) | non-zero to 750 nm (Chappuis) |
| J3 | H2O2 -> 2 OH | `H2O2+hv->OH+OH` | `H2O2_1.nc` = DATAJ1NEW/ABS/H2O2_jpl94_Kahan.abs outside 260-350 nm; inside 260-350 nm the JPL temperature polynomial of `src/cross_sections/h2o2-oh_oh.F90` (Nicovich & Wine 1988 form, chi(298 K) = 0.9859) | 1.0 (tuv-x constant; MCM: 1.0 > 230 nm) | cov 190-410 nm |
| J4 | NO2 -> NO + O(3P) | `NO2+hv->NO+O(3P)` | `NO2_1.nc` ("NO2 tint": rows at 220 and 294 K; 298 K clamps to the 294 K row). The netCDF header gives no file name (`../YYYY`); the values match the JPL-11 recommendation (Vandaele et al. 1998: 1.05e-20 at 254.8 nm, 6.44e-19 at 400 nm) | `data/quantum_yields/NO2_1.nc` = DATAJ1/YLD/NO2_jpl11.yld, 298 K row (248 K row unused); lower extrapolation `boundary` (phi = 1 below 300 nm) | xs cov 242-660 nm, qy data 300-422 nm. **sigma_phi = 0 above 422 nm is data** (JPL-11 yields: 0.02 at 420, 0.01 at 422, table ends) not truncation of the xs |
| J5 | NO3 -> NO + O2 | `NO3+hv->NO+O2` | `NO3_1.nc` = DATAJ1NEW/ABS/NO3_jpl11.abs | `NO3-NO+O2_1.nc` = DATAJ1NEW/YLD/NO3_jpl2011.qy1, 298 K row (also 230, 190 K) | xs cov 403-691 nm, qy 585-640 nm and 0 at 585 -> J5 non-zero only 586-640 nm (JPL-11 / Johnston et al. 1996) |
| J6 | NO3 -> NO2 + O(3P) | `NO3+hv->NO2+O(3P)` | as J5 | `NO3-NO2+O(3P)_1.nc` = NO3_jpl2011.qy2, 298 K row; configured lower extrapolation constant 1.0 (phi = 1 below 585 nm) | qy 0 above 640 nm (data end) |
| J7 | HONO -> OH + NO | `HNO2+hv->OH+NO` | `HNO2_1.nc` = DATAJ1NEW/ABS/HONO_jpl11.abs (JPL-11, Stutz/Bongartz based) | 1.0 | cov 184-396 nm. MCM used Bongartz et al. 1991/1994 (higher by ~20-30 % band-averaged; see bracket test) |
| J8 | HNO3 -> OH + NO2 | `HNO3+hv->OH+NO2` | `HNO3_1.nc` = DATAJ1NEW/ABS/HNO3_burk.abs (Burkholder et al. 1993): col 1 sigma(298 K), col 2 B; sigma = col1 exp(B (T-298)) (`hno3-oh_no2.F90`) -> col 1 at 298 K | 1.0 | cov 186-350 nm |
| J11 | HCHO -> H + HCO | `CH2O+hv->H+HCO` | `CH2O_1.nc` = DATAJ1NEW/CH2O/CH2O_jpl11.abs (Meller & Moortgat 2000): col 1 sigma(298), col 2 dsigma/dT (`ch2o.F90`) -> col 1 | `data/quantum_yields/CH2O_1.nc` = DATAJ1/CH2O/CH2O_jpl11.yld col 1 (H+HCO); lower extrapolation `boundary` (0.31 below 250 nm) | xs cov 226-375 nm, qy 250-361 nm (0 above 361) |
| J12 | HCHO -> H2 + CO | `CH2O+hv->H2+CO` | as J11 | CH2O_jpl11.yld col 2 (H2+CO at 1 atm) with the 330-360 nm pressure/temperature correction of `src/quantum_yields/ch2o.F90` evaluated at 298 K, M = 2.46e19 (reduces to col 2 at 300 K, 2.45e19) | as J11 |
| J13 | CH3CHO -> CH3 + HCO | `CH3CHO+hv->CH3+HCO` | `CH3CHO_1.nc` = DATAJ1NEW/CH3CHO/CH3CHO_jpl11.abs (Martinez et al. 1992) | `data/quantum_yields/CH3CHO_1.nc` = DATAJ1/CH3CHO/CH3CHO_iup.yld (IUPAC). **Column meaning from the data and `ch3cho-ch3_hco.F90`: col 1 = phi(CH4+CO), col 2 = phi(CH3+HCO) at 1 atm** (0.46/0.31 at 260 nm, 0/0.43 at 300 nm); Stern-Volmer form of the F90 at M = 2.46e19 (reduces to col 2 at 2.465e19) | xs cov 202-360 nm; qy data 260-330 nm, **no extrapolation configured in tuv-x -> 0 below 260 nm** (truncation, not data; irrelevant above 280 nm) |
| J14 | C2H5CHO -> C2H5 + HCO | `C2H5CHO+hv->C2H5+HCO` | `C2H5CHO_1.nc` = DATAJ1NEW/C2H5CHO/C2H5CHO_iup.abs | `C2H5CHO_1.nc` (qy) = DATAJ1NEW/C2H5CHO/C2H5CHO_iup.yld (IUPAC, 1 atm) with the Stern-Volmer scaling of `c2h5cho.F90` at M = 2.46e19 | xs cov 202-365 nm; qy data 294-340 nm, **no extrapolation configured -> 0 below 294 nm** (truncation; IUPAC yield there is ~0.9). Affects a 254/280-nm lamp, not UV-A lamps |
| J15 | n-C3H7CHO -> n-C3H7 + HCO | none | **no n-butanal data in tuv-x -> empty** | MCM: phi 0.21 (Tadic 2001) | JUICE fallback: scale from J14 (same Martinez 1992 family) |
| J16 | n-C3H7CHO -> C2H4 + CH3CHO | none | **empty** | MCM: phi 0.10 | as J15 |
| J17 | i-C3H7CHO -> i-C3H7 + HCO (the MCM page prints "n-C4H9 + HCO", a typo) | none | **no i-butanal data in tuv-x -> empty** | MCM: Chen et al. 2002 | fallback J14 |
| J18 | MACR -> CH2=CCH3 + HCO | `CH2=C(CH3)CHO+hv->Products` (lumped) | `Methacrolein_1.nc` = DATAJ1NEW/ABS/Methacrolein_jpl11.abs (Gierczak et al. 1997) | **MCM channel value phi = 1.95e-3** (Pinho et al. 2005 optimisation); tuv-x uses a lumped 0.01 | cov 250-395 nm |
| J19 | MACR -> CH2=C(CH3)CO + H | as J18 | as J18 | phi = 1.95e-3 (MCM: same as J18) | as J18 |
| J20 | C5HPALD1 -> CH3C(CHO)=CHCH2O + OH | none for HPALD | **proxy, as defined by the MCM itself:** sigma(MACR, Gierczak 1997) = `Methacrolein_1.nc` | phi = 1.0 (MCM assumption after Peeters et al. 2009) | Not HPALD data; the MCM l (7.60e-4) equals l(J18)/1.95e-3 = 7.60e-4, confirming this definition. Treat as a surrogate |
| J21 | CH3COCH3 -> CH3CO + CH3 | `CH3COCH3+hv->CH3CO+CH3` | `CH3COCH3_1.nc` = DATAJ1NEW/CH3COCH3/CH3COCH3_jpl11.abs (Gierczak et al. 1998): sigma = c1 (1 + T(c2 + T(c3 + T c4))), T clamped to 235-298 K (`src/cross_sections/acetone-ch3co_ch3.F90`) | `src/quantum_yields/acetone-ch3co_ch3.F90` = Blitz et al. (2004) pressure/temperature parametrisation, branch CH3CO, T = 298 K, M = 2.46e19; defaults 0.95 below 279 nm, 0 above 327 nm. phi(290/300/310/320) = 0.314/0.213/0.053/0.009 | cov 215-349 nm; the CO + 2 CH3 channel (Blitz "CO" branch) is not an MCM channel and is excluded |
| J22 | MEK -> CH3CO + C2H5 | `CH3COCH2CH3+hv->CH3CO+CH2CH3` | `CH3COCH2CH3_1.nc` = DATAJ1NEW/ABS/Martinez.abs (Martinez et al. 1992) | **phi = 0.400** = tuv-x `ch3coch2ch3-ch3co_ch2ch3.F90`, 1/(0.96 + 2.22e-3 p_torr) with p_torr = 760 M/2.69e19 = 695 (IUPAC pressure-dependent yield). The MCM page quotes 0.16 (Pinho et al. 2005), but the MCM l,m,n reproduce the bracket test only with ~0.35-0.40 (with 0.16: 0.56x / 0.33x of MCM; with 0.40: 1.41x / 0.83x). Switch `PHI_MEK_SOURCE = "mcm_page"` in the script, or multiply the column by 0.4, to get the page value | cov 202-355 nm |
| J23 | MVK -> CH3CH=CH2 + CO | `MVK+hv->Products` (lumped) | `MVK_1.nc` = DATAJ1NEW/ABS/MVK_jpl11.abs (Gierczak et al. 1997) | total phi = min(1, exp(-0.055 (lambda-308)) / (5.5 + 9.2e-19 M)) (`src/quantum_yields/mvk.F90`, Gierczak et al. 1997 / IUPAC) at M = 2.46e19, **split 50/50** between J23 and J24 (MCM: identical l,m,n). phi_total(300/320/340/360) = 0.055/0.018/0.006/0.002 | cov 250-395 nm; pressure dependent (1 atm) |
| J24 | MVK -> CH3CO + CH2=CH | as J23 | as J23 | 0.5 x total | as J23 |
| J31 | GLYOX -> H2 + 2 CO | `CHOCHO+hv->H2+CO+CO` | `CHOCHO_1.nc` = DATAJ1NEW/CHOCHO/glyoxal_jpl11.abs (Volkamer et al. 2005) | `CHOCHO-H2_CO_CO_1.nc` = glyoxal_jpl11.qy2 (JPL-11 1-atm channel yields built on Tadic et al. 2006 / Feierabend et al. 2009); `boundary` below 225 nm | xs cov 250-526 nm; qy 225-420 nm, 0 above 420 nm and this channel is 0 above ~355 nm. **Channel partition differs from the MCM, see bracket test** |
| J32 | GLYOX -> HCHO + CO | `CHOCHO+hv->CH2O+CO` | as J31 | `CHOCHO-CH2O_CO_1.nc` = glyoxal_jpl11.qy3 | as J31 |
| J33 | GLYOX -> 2 HCO | `CHOCHO+hv->HCO+HCO` | as J31 | `CHOCHO-HCO_HCO_1.nc` = glyoxal_jpl11.qy1 | as J31 |
| J34 | MGLYOX -> CH3CO + HCO | `CH3COCHO+hv->CH3CO+HCO` | `CH3COCHO_1.nc` = DATAJ1NEW/CH3COCHO/CH3COCHO_jpl11.abs | `src/quantum_yields/ch3cocho.F90` (JPL-11 form of Chen et al. 2000): phi0 = clip(1 - (lambda-380)/60, 0, 1); for 380-440 nm phi = phi0/(phi0 + 1.36e8 exp(-8793/lambda) x 760 M/2.456e19). phi(380/400/420/440) = 0.098/0.022/0.004/0 (1 atm) | cov 200-493 nm; strongly pressure dependent |
| J35 | BIACET -> 2 CH3CO | `CH3COCOCH3+hv->Products` | `CH3COCOCH3_1.nc` = DATAJ1NEW/CH3COCOCH3/biacetyl_horowitz.abs (Horowitz et al. 2001) | 0.158 constant (Plum et al. 1983) - identical in MCM and tuv-x | cov 206-492 nm |
| J41 | CH3OOH -> CH3O + OH | `CH3OOH+hv->CH3O+OH` | `CH3OOH_1.nc` = DATAJ1NEW/CH3OOH/CH3OOH_jpl11.abs (Vaghjiani & Ravishankara 1989) | 1.0 | cov 210-405 nm |
| J51 | CH3ONO2 -> CH3O + NO2 | `CH3ONO2+hv->CH3O+NO2` | `CH3ONO2_1.nc` = DATAJ1NEW/RONO2/CH3ONO2_jpl11.abs: col 1 sigma(298), col 2 B; sigma = col1 exp(B (T-298)) (`ch3ono2-ch3o_no2.F90`) -> col 1 | 1.0 | cov 190-344 nm |
| J52 | C2H5ONO2 -> C2H5O + NO2 | `C2H5ONO2+hv->C2H5O+NO2` | `C2H5ONO2_1.nc` = DATAJ1NEW/RONO2/C2H5ONO2_iup2006.abs (IUPAC 2006 average, incl. Talukdar 1997); same exp(B dT) form -> col 1 | 1.0 | cov 185-340 nm. tuv-x also has `CH3CH2ONO2_1.nc` (Talukdar only, 236-350 nm), not used |
| J53 | n-C3H7ONO2 -> n-C3H7O + NO2 | `nC3H7ONO2+hv->C3H7O+NO2` | `nC3H7ONO2_1.nc` = DATAJ1NEW/RONO2/nC3H7ONO2_iup2006.abs | 1.0 | cov 185-340 nm |
| J54 | i-C3H7ONO2 -> CH3C(O.)CH3 + NO2 | `CH3CHONO2CH3+hv->CH3CHOCH3+NO2` | `CH3CHONO2CH3_1.nc` = DATAJ1NEW/RONO2/isopropyl_talukdar.abs (Talukdar et al. 1997): col 1 sigma(298), col 2 B (`rono2.F90`) -> col 1 | 1.0 | cov 236-360 nm |
| J55 | t-C4H9ONO2 -> t-C4H9O + NO2 | `C(CH3)3(ONO2)+hv->C(CH3)(O.)+NO2` | **analytic fit** `src/cross_sections/t_butyl_nitrate.F90`: sigma = exp(-115.5 + 0.5307 lambda - 0.993e-3 lambda^2), valid 270-330 nm (Roberts & Fajer 1989) | 1.0 | 0 outside 270-330 nm by construction (tail beyond 330 nm missing) |
| J56 | NOA (CH3C(O)CH2ONO2) -> CH3C(O)CH2O + NO2 | `nitro_acetone+hv->Products` | **analytic fit** `src/cross_sections/nitroxy_acetone.F90`: sigma = exp(-156.8 + 0.7834 lambda - 1.365e-3 lambda^2), valid 284-335 nm (Roberts & Fajer 1989 / Barnes et al. 1993) | **MCM phi = 0.9** (Mueller et al. 2014); tuv-x assumes 1.0 | 0 outside 284-335 nm by construction |
| J57 | - | - | not listed on the MCM photolysis page and not in `mcm_converter.py`; no column | | |
| J61 | **unidentified** (l,m,n = 7.537e-4, 0.499, 0.266) | - | **empty.** Not on the MCM photolysis page; present with the same numbers in `mcm_converter.py` and in AtChem2 `mcm/photolysis-rates_v3.3.1`; used by **no** reaction in the full MCM v3.3.1 extract `MCM_tested_scenarios/IDL Input/mcm_full/chemistry/MCM_all.fac.rate` | | |

## Validation

### (a) Spot values (from `--validate`)

| quantity | table | expectation given | comment |
|---|---|---|---|
| sigma_phi(NO2) 254 nm | 1.037e-20 cm^2 | ~1.05e-20 (phi = 1) | OK (JPL-11 1.05e-20 at 254.8 nm, interpolated to 254.0) |
| sigma_phi(NO2) maximum | 6.30e-19 cm^2 at 390 nm; 1.20e-20 at 420, 5.9e-21 at 422, 0 at 423 nm | 5-6e-19 near 380-400, 0 above ~422 | OK; the 294-K Vandaele peak is 6.3-6.4e-19 |
| sigma(O3) 305 nm, phi(O1D) | sigma = 1.992e-19, phi = 0.900, product 1.79e-19 cm^2 | "~5e-19 x 0.9" | phi is right; sigma(O3, 305 nm) = 1.99e-19 is the Malicet/Brion (and JPL) value - the 5e-19 figure belongs to a shorter wavelength on the Hartley slope |
| sigma(H2O2) 254 nm | 7.02e-20 cm^2 | ~7e-20 | OK |
| sigma(HONO) 354 / 368 nm | 4.87e-19 / 4.36e-19 cm^2 | ~2e-19 | the JPL-11 (Stutz/Bongartz) table at ~1 nm resolution peaks at ~4.9e-19; 2e-19 corresponds to 5-10 nm band averages. Not an error |
| HCHO | sigma_phi(H+HCO) 305 nm = 3.60e-20, sigma_phi(H2+CO) 330 nm = 2.54e-20 | | consistent with JPL-11 |

### (b) Clear-sky bracket test, SZA = 30 deg

Extraterrestrial spectrum: `data/profiles/solar/sao2010.solref.converted` (Chance & Kurucz 2010, W m^-2 nm^-1 at 0.01 nm) converted to photons cm^-2 s^-1 nm^-1 and bin-averaged to 1 nm. Direct beam only, Beer-Lambert with sec(30 deg) = 1.1547: case (i) ozone only, 300 DU = 8.06e18 molec cm^-2 with sigma(O3, 295 K) from the table (J1 + J2); case (ii) ozone + Rayleigh tau_R = 0.00877 lambda_um^-4.05. MCM: J_i/J4 = l_i cos^m_i(30) exp(-n_i sec 30) / (same for J4). Case (i) has too much UV-B (no scattering), case (ii) too little (no diffuse light), so the MCM value should lie between or near them.

| J | MCM J_i/J4 | case i (O3 only) | case ii (O3+Rayleigh) | (i)/MCM | (ii)/MCM | verdict |
|---|---|---|---|---|---|---|
| J1 | 3.308e-03 | 4.003e-03 | 2.166e-03 | 1.21 | 0.65 | bracketed |
| J2 | 5.047e-02 | 4.460e-02 | 6.615e-02 | 0.88 | 1.31 | bracketed |
| J3 | 8.226e-04 | 8.628e-04 | 5.934e-04 | 1.05 | 0.72 | bracketed |
| J5 | 2.591e+00 | 2.305e+00 | 3.772e+00 | 0.89 | 1.46 | bracketed |
| J6 | 1.790e+01 | 1.861e+01 | 2.919e+01 | 1.04 | 1.63 | bracketed |
| J7 | 2.210e-01 | 1.640e-01 | 1.507e-01 | 0.74 | 0.68 | ~25 % low in both: JPL-11 HONO xs vs MCM's Bongartz 1991/94 (level difference, same shape) |
| J8 | 6.623e-05 | 2.443e-04 | 4.444e-05 | 3.69 | 0.67 | bracketed (deep-UV absorber: wide bracket) |
| J11 | 3.349e-03 | 3.830e-03 | 2.447e-03 | 1.14 | 0.73 | bracketed |
| J12 | 5.332e-03 | 5.770e-03 | 4.187e-03 | 1.08 | 0.79 | bracketed |
| J13 | 4.619e-04 | 5.761e-04 | 3.298e-04 | 1.25 | 0.71 | bracketed |
| J14 | 1.976e-03 | 1.987e-03 | 1.186e-03 | 1.01 | 0.60 | bracketed |
| J18 | 1.201e-04 | 1.093e-04 | 8.831e-05 | 0.91 | 0.74 | near (MCM phi applied) |
| J19 | 1.201e-04 | 1.093e-04 | 8.831e-05 | 0.91 | 0.74 | near |
| J20 | 6.158e-02 | 5.603e-02 | 4.529e-02 | 0.91 | 0.74 | near (MACR proxy, phi = 1) |
| J21 | 5.636e-05 | 5.315e-05 | 2.840e-05 | 0.94 | 0.50 | near/bracketed (Blitz 2004 in both) |
| J22 | 3.884e-04 | 5.468e-04 | 3.240e-04 | 1.41 | 0.83 | bracketed with phi = 0.400; with the page's 0.16 it would be 0.56 / 0.33 |
| J23 | 1.969e-04 | 2.393e-04 | 1.679e-04 | 1.22 | 0.85 | bracketed |
| J24 | 1.969e-04 | 2.393e-04 | 1.679e-04 | 1.22 | 0.85 | bracketed |
| J31 | 6.446e-03 | 1.881e-03 | 1.262e-03 | 0.29 | 0.20 | **genuine discrepancy** (see below) |
| J32 | 9.718e-04 | 3.357e-03 | 2.384e-03 | 3.45 | 2.45 | **genuine discrepancy** |
| J33 | 2.925e-03 | 7.969e-03 | 7.703e-03 | 2.72 | 2.63 | **genuine discrepancy** |
| J34 | 1.427e-02 | 1.479e-02 | 1.346e-02 | 1.04 | 0.94 | bracketed |
| J35 | 3.074e-02 | 2.961e-02 | 3.690e-02 | 0.96 | 1.20 | bracketed |
| J41 | 6.080e-04 | 6.121e-04 | 4.396e-04 | 1.01 | 0.72 | bracketed |
| J51 | 1.127e-04 | 5.076e-04 | 6.497e-05 | 4.50 | 0.58 | bracketed (UV-B absorber: wide bracket) |
| J52 | 1.311e-04 | 5.860e-04 | 1.084e-04 | 4.47 | 0.83 | bracketed |
| J53 | 1.734e-04 | 6.767e-04 | 1.335e-04 | 3.90 | 0.77 | bracketed |
| J54 | 2.932e-04 | 3.003e-04 | 1.807e-04 | 1.02 | 0.62 | bracketed |
| J55 | 8.356e-04 | 7.276e-04 | 4.380e-04 | 0.87 | 0.52 | near (fit truncated at 330 nm) |
| J56 | 3.110e-03 | 2.877e-03 | 1.785e-03 | 0.92 | 0.57 | near (phi 0.9 applied) |

Reference values quoted in the task (J1/J4 3.3e-3, J7/J4 0.22, J3/J4 8.2e-4, J11/J4 3.3e-3, J12/J4 5.3e-3, J13/J4 4.6e-4, J21/J4 5.6e-5, J41/J4 6.1e-4, J5/J4 2.6, J6/J4 17.9) are reproduced by the MCM column.

**Glyoxal (J31-J33) - documented, not "fixed".** The *total* glyoxal photolysis agrees: (J31+J32+J33)/J4 = 1.034e-2 (MCM) vs 1.321e-2 (case i) and 1.135e-2 (case ii), i.e. bracketed within 1.1-1.3x. Only the channel partition differs. JPL-11 (tuv-x) at 1 atm: HCO+HCO is the main channel at 300-420 nm, H2+2CO is zero above ~355 nm, and all yields are zero above 420 nm, so the strong 440-460 nm band contributes nothing. The MCM (Tadic et al. 2006 as implemented in MCM v3.2/3.3.1) assigns the dominant, weakly SZA-dependent share (m = 0.13, n = 0.201, i.e. visible-light driven) to J31 (H2+2CO) and a UV-type dependence (m = 0.644) to J33. No mapping or unit error is involved (the three tuv-x channel files were checked against the qy1/qy2/qy3 column meaning in `photolysis_rate_constants.json` and against the values: qy1 peaks at 0.41 near 345 nm, qy2 falls from 0.65 at 225 nm to 0 at 365 nm, qy3 peaks at 0.46 near 305 nm). If JUICE wants MCM-consistent channel ratios for a lamp, use the sum of the three columns as sigma x phi_total and apply the MCM channel fractions at SZA 30 deg (J31 : J32 : J33 = 0.623 : 0.094 : 0.283), noting that the lamp-weighted total will still differ from the MCM for lamps emitting mainly above 420 nm.

**MEK (J22).** The MCM page quotes phi = 0.16 (Pinho et al. 2005) but the l,m,n of J22 correspond to phi ~ 0.35-0.40 with the Martinez et al. (1992) cross section; the table therefore uses the IUPAC pressure-dependent yield as implemented in tuv-x (0.400 at 1 atm). This is a known inconsistency between the MCM web table and its parameters.

**HONO (J7).** 25-30 % below the MCM in both cases. The shape of the JPL-11 and Bongartz spectra is the same; the difference is the absolute level of the data sets the two compilations chose. Left as is (it is the JPL/IUPAC recommendation).

### (c) Runtime

`build_sigma_phi.py <tuv-x> --validate` completes in about 3 s (netCDF4 reads + numpy), well under the 2-minute limit.

## Known limitations

* Temperature dependence of all cross sections is frozen at 298 K (NO2 at 294 K, O3 at 295 K); pressure dependence at 1 atm (HCHO H2+CO, CH3CHO, C2H5CHO, acetone, MVK, methylglyoxal, MEK). Chambers operating at reduced pressure need the F90 parametrisations re-evaluated (all are in the script with `M_AIR`, `T_K` as module constants).
* tuv-x's default *zero* extrapolation outside quantum-yield tables truncates CH3CHO (< 260 nm) and C2H5CHO (< 294 nm); the HCHO and glyoxal yields use the configured `boundary` extrapolation downwards.
* J15, J16, J17, J61 are empty; J20 is a MACR surrogate by MCM definition; J18/J19/J56 use MCM constant yields; J22 uses the IUPAC/tuv-x yield rather than the page's 0.16.
