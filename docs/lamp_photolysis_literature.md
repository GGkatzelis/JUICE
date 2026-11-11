# Photolysis frequencies of lamp-lit simulation chambers — literature survey

Compiled 2026-11-02 for JUICE chamber profiles (`setups/chambers/<id>.chamber.json`).
Scope: published j(NO2) and ratios j(X)/j(NO2) for artificial light sources
(blacklights, 254/185 nm Hg lamps, Xe/Ar arcs, LEDs), plus clear-sky reference values.

How to read the labels used throughout:

| Label | Meaning |
|---|---|
| **[A]** | measured by chemical actinometry (NO2 photo-stationary state or NO2 decay in N2) |
| **[S]** | derived from a measured lamp spectrum × σ × φ (spectroradiometer, scaled or absolute) |
| **[C]** | computed in the source paper from a model (TUV, MCM, KinSim) |
| **[R]** | computed by me in this review from quoted numbers — not in the source |
| **[?]** | secondary citation, snippet, or value I could not verify in the primary text |

Units: all j in s⁻¹ unless written `min⁻¹`. 1 min⁻¹ = 1.667e-2 s⁻¹. Photon fluxes in photons cm⁻² s⁻¹.

Notation: j(O1D) = j(O3 → O(1D)+O2); HCHO_r = HCHO → H+HCO; HCHO_m = HCHO → H2+CO;
NO3_a = NO3 → NO+O2; NO3_b = NO3 → NO2+O(3P). MCM J-numbers: J1=O1D, J3=H2O2, J4=NO2, J5/J6=NO3_a/b,
J7=HONO, J11/J12=HCHO_r/m, J13=CH3CHO, J21=acetone, J41=CH3OOH.

---

## 1. Summary table

Ratios are relative to j(NO2) of the same light source. Where a cell says "sun-like" the source
only demonstrated that j(O1D)/j(NO2) matches sunlight; other ratios then follow the clear-sky row.

| Light source (typical lamps) | j(NO2) full power, s⁻¹ | j(O1D)/j(NO2) | j(HONO)/j(NO2) | j(H2O2)/j(NO2) | j(HCHO_r)/j(NO2) | j(HCHO_m)/j(NO2) | other ratios | basis / sources |
|---|---|---|---|---|---|---|---|---|
| **Clear sky, SZA ≈ 30° (Jülich noon, late June)** | 8–9e-3 outside; 6–7e-3 inside SAPHIR | 3.3e-3 (MCM, SZA 30°) [C]; 3.0–3.5e-3 SAPHIR site [R from Bohn & Zilken Fig. 14] | 0.22 | 8.2e-4 | 3.3e-3 | 5.3e-3 | CH3CHO 4.6e-4; acetone 5.6e-5; CH3OOH 6.1e-4; NO3_a 2.6; NO3_b 17.9 | §2 (MCM l,m,n; Bohn & Zilken 2005; Bohn et al. 2005) |
| **350 nm blacklight** (Sylvania F40/350BL, GE F40BL, UVA-340/350; 310–400 nm, peak 350–355 nm) | 1.4e-3 (128 lamps, 10 m³ [?]) … 2.0e-3 (NIES, 276 lamps) … 2.2e-3 (ILMARI, 47 lamps) … 4.3e-3 (UCR EPA, 80% power) … 4.7e-3 (GTEC) … 6e-3 (Aerodyne OFR, 4 lamps) … 2.5e-2 (old Caltech 28 m³, 300 lamps [?]) | **no actinometric value found.** Crude estimate [R]: 0.5–3e-3, extremely sensitive to the <315 nm tail of the phosphor; Carter: "abnormally low" vs sunlight for 300–340 nm-poor lamps | ≈0.3–0.5 [R, from σ overlap; Platt gets 0.27 for a 368 nm lamp] | 5.3e-3 (NIES, 350 nm/40 nm FWHM) [S]; ≈3e-4 implied by Caltech H2O2 experiments [R?]; use 1e-3 ×/÷ 4 | ≈2–4e-3 [R] | ≈4–6e-3 [R] | NO3 ≈ 0 (no emission >410 nm); CH3CHO ≈ 2–4e-4; acetone ≈ 2–4e-5; CH3OOH ≈ 1e-3 [all R] | §3.4–3.9, §3.13–3.15 |
| **365 nm "BLB" / actinic blacklight** (Sylvania F36/T8 BLB, GE F40BLB, Philips TL-K 40W/05 (peak ≈365 nm); 340–400 nm, peak 365–369 nm) | 3.7e-3 / 5.8e-3 / 1.03e-2 for 54/108/160 lamps (ILMARI) [A]; 8.2e-3 (GIG-CAS, 135 × 60 W, peak 369 nm) [A]; 1.9e-3 (HIRAC, 24 × TL-K 40W/05) [A] | ≈0 (no emission <340 nm) — Carter/Leskinen; use ≤1e-4 | ≈0.3 [R] | "only to a very small degree" (Leskinen); use ≤2e-4 | ≈1e-3 [R] | ≈3e-3 [R] | NO3 ≈ 0 | §3.7, §3.14 |
| **368 nm UVA solarium lamp** (Philips Cleo Performance 100 W, 310–430 nm; PSI mobile chamber, 40 lamps) | 8.0 ± 0.7e-3 [A] | **5.4e-2 [S, Platt Table 5] — suspiciously high, see §3.11** | 0.27 [S] | not given; ≈1e-3 [R from CH3OOH] | 1.9e-3 [S] | 4.7e-3 [S] | CH3CHO 1.4e-4; CH3OOH 5.4e-4; glyoxal 8.9e-4; MVK 4.3e-4; acetone 1.3e-5; CH3ONO2 3.6e-4; PAN 5e-5; NO3_a 0; NO3_b 0.56 [all S] | Platt et al. 2013 (§3.11) |
| **254 nm low-pressure Hg, 185 nm blocked (OFR254)** (BHK 82-9304-03 with quartz sleeve; Light Sources GPH436T5L) | per 1e15 photons cm⁻² s⁻¹ at 254 nm: j(NO2) = 1.05e-5 from the 254 nm line alone [R]; actinometry gives ≈13× more because of the weak 313/365/405/436 nm Hg lines (Lambe 2019) | ≈1000 (254 nm line only) / ≈75 (with minor-line j(NO2)) [R] | σ(HONO,254) not verified; order 1 [R?] | 6.7 / 0.5 [R] | ≈0.3–0.7 [R, σ ≈ 0.5–1e-20, φ_r ≈ 0.3] | ≈0.5–1 [R] | j(O1D) = 1.0e-2 per 1e15; j(H2O2) = 7e-5; j(CH3OOH) = 3.3e-5; j(acetone) ≈ 3e-5; j(HNO3) = 1.2e-5; j(NO3) ≈ 0 | §3.2, §4 |
| **185 + 254 nm Hg (OFR185)** (BHK 82-9304-03 / GPH436T5VH, Suprasil) | as OFR254; F185/F254 ≈ 0.01 (BHK, Li 2015) to 0.10 (GPH436T5VH, Rowe 2020) | as OFR254 | — | — | — | — | j(H2O) = 7.1e-20·F185; j(O2) ≈ 1e-20·F185 (path dependent); j(N2O) = 1.4e-19·F185; j(O3,185) ≈ 6e-19·F185 | §3.2 |
| **Xe arc, 9 mm Pyrex filtered** (CESAM, 3 × 4 kW; MAC, 2 × 6 kW XBO + 112 halogen) | 4.2 ± 0.2e-3 (CESAM) [A]; 2.25 ± 0.4e-3 (MAC) [A]; 2.0e-3 (PSI 27 m³, Xe + blacklights) [?] | 3.3e-3 (CESAM: 1.4e-5/4.2e-3) [A/S]; ratio j(NO2)/j(O1D) = 298 vs 260 for TUV ground level | sun-like (use clear-sky row) | sun-like | sun-like | sun-like | NO3: sun-like unless visible is filtered | §3.8–3.9 |
| **Ar arc solar simulator** (UCR EPA, 200 kW) | 4.2–4.3e-3 (0.25–0.26 min⁻¹ at 80% power) [A]; 0.09–0.26 min⁻¹ used in Warren 2008 | sun-like per Carter; table of ratios not retrievable (UCR site down) | sun-like | sun-like | sun-like | sun-like | — | §3.4 |
| **LED solar simulator** (AIDA: 3000 LEDs, 305–528 nm) | 1.58 ± 0.21e-3 [A] | not reported; depends on 305–320 nm LED share | — | — | — | — | 2,3-pentanedione 1.03e-4 [A] | Vallon et al. 2022 (§3.12) |

---

## 2. Clear-sky reference

### 2.1 MCM photolysis parameterisation (Jenkin et al. 1997; Saunders et al. 2003)

The MCM website (https://cri.york.ac.uk/MCM/rates/photolysis) tabulates `J = l·cos(χ)^m·exp(−n·sec χ)`
(χ = solar zenith angle) for 34 photolysis processes. The parameters are stated to be consistent with
Saunders et al. (2003), ACP 3, 161–180, DOI 10.5194/acp-3-161-2003 (MCM v3 protocol), which itself
builds on Jenkin et al. (1997), Atmos. Environ. 31, 81–104, DOI 10.1016/S1352-2310(96)00105-7.
They represent **clear-sky, surface, low-albedo, fixed-ozone-column** conditions from a two-stream
radiative transfer calculation — not any particular site or season. They are therefore directly usable
as "ratio to j(NO2) under clear sky", and they are what JUICE should use as the solar reference.

Parameters (as listed on the MCM page, 2026-11-02):

| J | process | l | m | n |
|---|---|---|---|---|
| J1 | O3 → O(1D) | 6.073e-5 | 1.743 | 0.474 |
| J2 | O3 → O(3P) | 4.775e-4 | 0.298 | 0.080 |
| J3 | H2O2 → 2 OH | 1.041e-5 | 0.723 | 0.279 |
| J4 | NO2 → NO + O(3P) | 1.165e-2 | 0.244 | 0.267 |
| J5 | NO3 → NO + O2 | 2.485e-2 | 0.168 | 0.108 |
| J6 | NO3 → NO2 + O(3P) | 1.747e-1 | 0.155 | 0.125 |
| J7 | HONO → OH + NO | 2.644e-3 | 0.261 | 0.288 |
| J8 | HNO3 → OH + NO2 | 9.312e-7 | 1.230 | 0.307 |
| J11 | HCHO → H + HCO | 4.642e-5 | 0.762 | 0.353 |
| J12 | HCHO → H2 + CO | 6.853e-5 | 0.477 | 0.323 |
| J13 | CH3CHO → CH3 + HCO | 7.344e-6 | 1.202 | 0.417 |
| J21 | CH3COCH3 → CH3CO + CH3 | 7.992e-7 | 1.578 | 0.271 |
| J22 | MEK | 5.804e-6 | 1.092 | 0.377 |
| J41 | CH3OOH → CH3O + OH | 7.649e-6 | 0.682 | 0.279 |

(The page also lists J14–J20, J23, J24, J31–J35, J51–J56; see the site.)

**Computed from these parameters [R]** (absolute J in s⁻¹ and ratio to J4):

| SZA | J4 (NO2) | J1/J4 (O1D) | J7/J4 (HONO) | J3/J4 (H2O2) | J11/J4 (HCHO_r) | J12/J4 (HCHO_m) | J13/J4 (CH3CHO) | J21/J4 (acetone) | J41/J4 (CH3OOH) | J5/J4 (NO3_a) | J6/J4 (NO3_b) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0° | 8.92e-3 | 4.24e-3 | 0.222 | 8.83e-4 | 3.66e-3 | 5.56e-3 | 5.43e-4 | 6.83e-5 | 6.49e-4 | 2.50 | 17.3 |
| 20° | 8.64e-3 | 3.81e-3 | 0.222 | 8.56e-4 | 3.52e-3 | 5.46e-3 | 5.06e-4 | 6.29e-5 | 6.31e-4 | 2.54 | 17.5 |
| **30°** | **8.26e-3** | **3.31e-3** | **0.221** | **8.23e-4** | **3.35e-3** | **5.33e-3** | **4.62e-4** | **5.64e-5** | **6.08e-4** | **2.59** | **17.9** |
| 40° | 7.70e-3 | 2.67e-3 | 0.220 | 7.74e-4 | 3.10e-3 | 5.14e-3 | 4.02e-4 | 4.78e-5 | 5.75e-4 | 2.68 | 18.5 |
| 50° | 6.90e-3 | 1.95e-3 | 0.218 | 7.10e-4 | 2.77e-3 | 4.86e-3 | 3.27e-4 | 3.78e-5 | 5.31e-4 | 2.83 | 19.5 |
| 60° | 5.77e-3 | 1.22e-3 | 0.215 | 6.26e-4 | 2.34e-3 | 4.48e-3 | 2.40e-4 | 2.70e-5 | 4.73e-4 | 3.09 | 21.2 |
| 70° | 4.11e-3 | 5.70e-4 | 0.210 | 5.16e-4 | 1.78e-3 | 3.89e-3 | 1.46e-4 | 1.62e-5 | 3.96e-4 | 3.68 | 25.0 |

Jülich (50.9° N) summer-solstice noon has SZA ≈ 27.5°, so the SZA 30° row is the right "SAPHIR summer noon" reference.
Note how strongly J1/J4 falls with SZA (factor 7 from 0° to 70°) while J7/J4 is flat — this is the physical
reason the j(O1D)/j(NO2) ratio is the diagnostic quantity for any light source.

### 2.2 SAPHIR / Jülich measurements

- **Bohn, B. and Zilken, H.** (2005): Model-aided radiometric determination of photolysis frequencies in a sunlit atmosphere simulation chamber, *Atmos. Chem. Phys.* 5, 191–206, DOI 10.5194/acp-5-191-2005.
  Abstract: "Under overcast conditions we obtain fractions of 0.74 and 0.67 of external values for photolysis frequencies j(NO2) and j(O1D), respectively. On a clear sky summer day these values are time-dependent within ranges 0.65–0.86 and 0.60–0.73." Text: on 28 July 2002 the inside/outside ratios "vary between 0.65 and 0.86 with a mean of 0.78 for j(NO2) and between 0.60 and 0.73 with a mean of 0.68 for j(O1D)." Fig. 14 (external, clear sky, 28 July 2002, SZA_min 32°): j(NO2) peaks at ≈ 9e-3 s⁻¹ (axis to 10e-3) and j(O1D) at ≈ 2.5–3e-5 s⁻¹ [values read off the figure, R]. Direct sun contributes up to 50 % of j(NO2) and 30 % of j(O1D) at noon.
  → **Outside: j(O1D)/j(NO2) ≈ 3.0–3.5e-3; inside SAPHIR ≈ 0.68/0.78 × that ≈ 2.6–3.0e-3** [R]. This agrees with the MCM SZA 30° value (3.3e-3).
- **Bohn, B., Rohrer, F., Brauers, T., Wahner, A.** (2005): Actinometric measurements of NO2 photolysis frequencies in the atmosphere simulation chamber SAPHIR, *Atmos. Chem. Phys.* 5, 493–503, DOI 10.5194/acp-5-493-2005. NO2 actinometry with 20–50 ppb NO2 in the dark chamber, photo-stationary and time-dependent (FACSIMILE) evaluation; calibration factor for the radiometry-model approach h = 1.34 ± 0.10, overall accuracy of chamber j(NO2) 13 %; NO + O3 rate constant (1.93 ± 0.12)e-14 cm³ s⁻¹ at 298 K re-determined. Fig. 2/9: external spectroradiometer j(NO2) on 28 July 2002 ≈ 8–9e-3 s⁻¹ at noon; 10 March 2002 ≈ 5–6e-3; 10 Dec 2002 ≈ 2e-3.
- **Bohn, B., et al.** (2008): Photolysis frequency measurement techniques: results of a comparison within the ACCENT project, *Atmos. Chem. Phys.* 8, 5373–5391, DOI 10.5194/acp-8-5373-2008 (Jülich intercomparison; j(NO2) spectroradiometers agree within ±3 %).
- Related method papers: Hofzumahaus, A., et al. (1999), *Appl. Opt.* 38, 4443 (solar actinic flux spectroradiometry; j(O1D), j(NO2), j(HCHO)); Shetter, R. E. and Müller, M. (1999), *J. Geophys. Res.* 104, 5647 (ratios at Mauna Loa); Kraus, A., Hofzumahaus, A. (1998), *J. Atmos. Chem.* 31, 161 (O3, NO2, HCHO, H2O2, HONO by UV spectroradiometry).

### 2.3 Other clear-sky anchors quoted in chamber papers

- Wang et al. 2011 (CESAM): TUV/NCAR, equator, 21 June noon: j(NO2) = 1.27e-2, j(O1D) = 4.8e-5 (ratio 3.8e-3); ground-level ratio j(NO2)/j(O1D) = 260.
- Shao et al. 2022 (MAC): spectroradiometer, Manchester, clear-sky summer: j(NO2) = 7e-3.
- Vallon et al. 2022 (AIDA): Quick TUV, Karlsruhe, max clear winter / summer day: j(NO2) = 3.4e-3 / 6.0e-3.
- Wang et al. 2014 (GIG-CAS): TUV, Guangzhou, 22 June noon: j(NO2) = 0.58 min⁻¹ (9.7e-3).
- Glowacki et al. 2007 (HIRAC) citing Carter et al. (1995): a chamber light source "is able to produce j(NO2) values on the order of 5e-3 s⁻¹", i.e. solar-noon-like.

---

## 3. Per-chamber / per-group blocks

### 3.1 CU Boulder — Jimenez / Ziemann environmental chamber facility

- **Krechmer, J. E., Pagonis, D., Ziemann, P. J., Jimenez, J. L.** (2016): Quantification of gas-wall partitioning in Teflon environmental chambers using rapid bursts of low-volatility oxidized species generated in situ, *Environ. Sci. Technol.* 50, 5757–5765, DOI 10.1021/acs.est.6b00606 (8 m³ chamber; full text paywalled — lamp details and j(NO2) not retrieved).
- **Krechmer, J. E., Day, D. A., Ziemann, P. J., Jimenez, J. L.** (2017): Direct measurements of gas/particle partitioning and mass accommodation coefficients in environmental chambers, *Environ. Sci. Technol.* 51, 11867–11875, DOI 10.1021/acs.est.7b02144 [DOI from memory — verify].
- **ICARUS** chamber record "1-White (1W)", Jimenez Lab (https://icarus.ucdavis.edu/chamber/534): 20 m³ FEP, UV lamps **F40/350BL/ECO (40 W)** and **FR48T12/BL-350/VHO/182 (115 W)**, UV peak 350 nm (FWHM 20 nm), visible T5HO 54 W tubes (490 nm); light-flux characterisations #38/#39/#40/#42 (2021) exist but their j values are not shown in the page text. Linked experiment sets: Krechmer ES&T 2017 "fast burst photooxidation", Liu 2019.
- **Result: no published j(NO2) or ratios found for this chamber in open sources.** The lamp type is the standard 350 nm blacklight; use the generic 350 nm profile (§5) and the ICARUS spectrum files if an absolute value is needed.
- Zhang et al. 2014 PNAS (vapor wall loss) is a **Caltech** chamber paper (see §3.3), not CU.

### 3.2 CU Boulder / Jimenez group and partners — oxidation flow reactors with low-pressure Hg lamps

Primary sources:

- **Li, R., Palm, B. B., Ortega, A. M., Hlywiak, J., Hu, W., Peng, Z., Day, D. A., Knote, C., Brune, W. H., de Gouw, J. A., Jimenez, J. L.** (2015): Modeling the radical chemistry in an oxidation flow reactor: radical formation and recycling, sensitivities, and the OH exposure estimation equation, *J. Phys. Chem. A* 119, 4418–4432, DOI 10.1021/jp509534k. Characterised the emission of BHK 82-9304-03 lamps vs. voltage; F185/F254 "around 1 % with a slightly decreasing trend with lamp setting" (as summarised in Peng & Jimenez 2020); 254 nm flux "80–250 times larger than at 185 nm" (Peng et al. 2016).
- **Peng, Z., Day, D. A., Stark, H., Li, R., Lee-Taylor, J., Palm, B. B., Brune, W. H., Jimenez, J. L.** (2015): HOx radical chemistry in oxidation flow reactors with low-pressure mercury lamps systematically examined by modeling, *Atmos. Meas. Tech.* 8, 4863–4890, DOI 10.5194/amt-8-4863-2015. Sect. 2: "The range of photon fluxes studied at 185 and 254 nm are 10^11–10^14 and 4.2e13–8.5e15 photons cm⁻² s⁻¹". Table 1 pairs: L = 1e11 / 4.2e13; M = 1e13 / 1.4e15; H = 1e14 / 8.5e15 (185 / 254 nm). Cross sections from JPL, uncertainty via Hébrard et al.; "most of the uncertainty is contributed by photolysis rates of O3, O2, and H2O". Supplement Fig. S18: 254 nm emission vs. lamp voltage (no emission below 17.025 V).
- **Peng, Z., Day, D. A., Ortega, A. M., Palm, B. B., Hu, W., Stark, H., Li, R., Tsigaridis, K., Brune, W. H., Jimenez, J. L.** (2016): Non-OH chemistry in oxidation flow reactors for the study of atmospheric chemistry systematically examined by modeling, *Atmos. Chem. Phys.* 16, 4283–4305, DOI 10.5194/acp-16-4283-2016. Same flux grid (L: 1e11/4.17e13, M: 1e13/1.45e15, H: 1e14/8.51e15; one field case 7.9e11/2.04e14). Tables S1–S3 of the supplement list σ(185 nm) and σ(254 nm) for many VOCs (unit quantum yield assumed unless stated). Photolysis at 254 nm is "a concern mainly" for conjugated species (σ > 1e-18 cm²).
- **Peng, Z., Jimenez, J. L.** (2017): Modeling of the chemistry in oxidation flow reactors with high initial NO, *Atmos. Chem. Phys.* 17, 11991–12010, DOI 10.5194/acp-17-11991-2017 (same flux grid; NOy chemistry with HONO/NO2/HNO3).
- **Peng, Z., Jimenez, J. L.** (2020): Radical chemistry in oxidation flow reactors for atmospheric chemistry research, *Chem. Soc. Rev.* 49, 2570–2616, DOI 10.1039/C9CS00766K. Review of all OFR modes (OFR185, OFR254-X, -iNO, -iN2O, OFR369-i(iPrONO)); UV grid as above; "UVA (369 nm) light intensity should be very strong (on the orders of 10^15–16 photons cm⁻² s⁻¹)" for nitrite photolysis; OHexp estimation equations (Table in Sect. 6); KinSim mechanism files for all modes in the ESI. KinSim (Peng & Jimenez 2019, *J. Chem. Educ.* 96, 806) runs in Igor Pro; the OFR mechanism and lamp-flux inputs are distributed with the review and on the PAM users' wiki (https://sites.google.com/site/pamwiki/ — lamp spectra and calibration notes are kept there, not retrieved here).
- **Lambe, A. T., et al.** (2011): Characterization of aerosol photooxidation flow reactors: heterogeneous oxidation, secondary organic aerosol formation and cloud condensation nuclei activity measurements, *Atmos. Meas. Tech.* 4, 445–461, DOI 10.5194/amt-4-445-2011. PAM: "Four mercury lamps (BHK Inc.) with peak emission intensity at λ = 254 nm are mounted in Teflon-coated quartz cylinders"; TPOT: ozone-free Pen-Ray lamp, 254 nm; lamp voltage 0–110 V. Rowe et al. 2020 cite this work for **I254,max = (3.5 ± 0.7)e15 photons cm⁻² s⁻¹ for two lamps at maximum output**.
- **Lambe, A. T., Krechmer, J. E., Peng, Z., Casar, J. R., Carrasquillo, A. J., Raff, J. D., Jimenez, J. L., Worsnop, D. R.** (2019): HOx and NOx production in oxidation flow reactors via photolysis of isopropyl nitrite, isopropyl nitrite-d7, and 1,3-propyl dinitrite at λ = 254, 350, and 369 nm, *Atmos. Meas. Tech.* 12, 299–311, DOI 10.5194/amt-12-299-2019. **The most useful lamp-calibration paper for JUICE:**
  - Lamps: 4 × GPH436T5L (254 nm, Light Sources Inc.), 4 × F436T5/BL/4P-350 (350 nm, Aerodyne), 4 × F436T5/BLC/4P-369 (369 nm, Aerodyne); emission spectra from the manufacturer in Fig. S1.
  - Table 1 cross sections (Atkinson et al. 2004), averaged over each lamp spectrum: **σ(NO2) = 1.05e-20 (254 nm), 4.70e-19 (350 nm), 5.60e-19 (369 nm) cm²**.
  - NO2 actinometry in the OFR (τ = 98 s) gave "maximum I254 = 8.6e16, I350 = 6.3e15, and I369 = 6.5e15 photons cm⁻² s⁻¹. While I350 and I369 values were in agreement with values calculated from lamp manufacturer specifications (I350 = 5.8e15 and I369 = 6.2e15) within uncertainties, I254 obtained from our calibration was ≈13 times larger than expected. We hypothesize that this discrepancy was due to the presence of additional minor mercury lines (e.g., λ ≈ 313, 365, 405) that induce NO2 photolysis … Thus, we instead assume maximum I254 = 6.5e15 photons cm⁻² s⁻¹ based on manufacturer specifications."
  - "In our OFR, jNO2,max = 0.36 min⁻¹ at λ = 350 nm" (= 6.0e-3 s⁻¹ with 4 lamps). A blacklight chamber they compare with had "jNO2 = 0.27 min⁻¹, assumed 350 nm wavelength".
  - **Consequence for a 254 nm reactor:** j(NO2) measured by actinometry is dominated by the minor UVA/visible Hg lines, not by the 254 nm line. For j(O1D)/j(NO2) one must therefore decide whether "j(NO2)" means the 254 nm-line value (σ·F254) or the actinometric value (≈13× larger). See §4.3 for both.
- **Rowe, J. P., Lambe, A. T., Brune, W. H.** (2020): Technical note: Effect of varying the λ = 185 and 254 nm photon flux ratio on radical generation in oxidation flow reactors, *Atmos. Chem. Phys.* 20, 13417–13424, DOI 10.5194/acp-20-13417-2020 (not AMT). Lamp type A = GPH436T5VH/4P (Light Sources Inc.), type 214 quartz, "a standard component of the Aerodyne PAM OFR"; types C–E (GPH436T5L/VH/4P 90/10, 96/4, 98.5/1.5) fuse 185-blocking and 185-transmitting quartz. Model inputs: "I254,max = (3.5 ± 0.7)e15 photons cm⁻² s⁻¹ for two lamps operated at maximum output (Lambe et al. 2011)"; "I185,max : I254,max ≈ 0.10 for lamp types A and B"; OHexp equations valid for I185:I254 ≈ 0.001 to 0.1. Measured lamp spectra are not published as data files in this note (Fig. S2 shows relative output).
- Other field/lab OFR papers from the group (Palm et al. 2016 ACP 16, 2943; Ortega et al. 2016 ACP 16, 7411; Hu et al. 2016 ACP 16, 11563) use the same flux ranges and the Li/Peng estimation equations; none adds lamp-spectrum data beyond the above.

Computed OFR photolysis frequencies [R] — see §4.3 for the table and the cross sections used.

### 3.3 Caltech (Seinfeld / Flagan / Wennberg / Ng-era)

- **Cocker, D. R., Flagan, R. C., Seinfeld, J. H.** (2001): State-of-the-art chamber facility for studying atmospheric aerosol chemistry, *Environ. Sci. Technol.* 35, 2594–2601, DOI 10.1021/es0019169. Dual 28 m³ FEP chambers, 300 blacklights (same type as Carter's). A search snippet of the PDF (intra.engr.ucr.edu/~dcocker/JA7.pdf, host unreachable from here) states "the measured k1 value in the chamber at full light intensity is 1.5 min⁻¹", and Leskinen et al. (2015) quote the same number ("1.5 min⁻¹ in the 28 m³ chamber at the California Institute of Technology (Cocker et al., 2001)"). **1.5 min⁻¹ = 2.5e-2 s⁻¹, i.e. ≈3× solar noon [?]** — treat as unverified; it is plausible for 300 lamps around 28 m³ but most later Caltech work ran at reduced intensity.
- ICARUS light-flux characterisations #18–#21 (2013-12-16, "Seinfeld chambers"): "UV black lamps (350 nm centered), 320 total lights", spectroradiometer irradiance converted to actinic flux and "scaled using H2O2 photolysis"; j fields blank; spectrum CSVs downloadable (JHSflux50_Dec2013.csv, JHSflux_onebulb_Dec2013.csv).
- **Zhang, X., Cappa, C. D., Jathar, S. H., McVay, R. C., Ensberg, J. J., Kleeman, M. J., Seinfeld, J. H.** (2014): Influence of vapor wall loss in laboratory chambers on yields of secondary organic aerosol, *PNAS* 111, 5802–5807, DOI 10.1073/pnas.1404727111 — current 24 m³ dual chambers with 350 nm blacklights; SI not retrievable (publisher block); no j value collected.
- **Huang, Y., Coggon, M. M., Zhao, R., Lignell, H., Bauer, M. U., Flagan, R. C., Seinfeld, J. H.** (2017): The Caltech Photooxidation Flow Tube reactor, *Atmos. Meas. Tech.* 10, 839–867, DOI 10.5194/amt-10-839-2017: UVA blacklights centred at 350 nm, NO2 actinometry **j(NO2) = 5.36e-3 s⁻¹** (flow tube, not the chamber).
- Current Ng-group page (nggroup.caltech.edu): blacklights 300–400 nm peak 354 nm plus "sunshine" fluorescent 300–900 nm; no j values.
- Lambe et al. 2011 note that the Caltech chamber OH source is "UV photolysis of H2O2 (λ = 350 nm)". Typical Caltech H2O2 experiments (≈2 ppm H2O2 → [OH] ≈ 2e6 cm⁻³ with OH reactivity ≈ 100 s⁻¹ dominated by H2O2 itself) imply j(H2O2) ≈ 2e-6 s⁻¹ [R, inference], i.e. j(H2O2)/j(NO2) ≈ 3e-4 if j(NO2) ≈ 6e-3 — an order of magnitude below the NIES estimate (§3.15). This spread is the reason for the large uncertainty attached to j(H2O2) for blacklights in §5.

### 3.4 UC Riverside CE-CERT (Carter, Cocker)

- **Carter, W. P. L., Cocker, D. R., Fitz, D. R., Malkina, I. L., Bumiller, K., Sauer, C. G., Pisano, J. T., Bufalino, C., Song, C.** (2005): A new environmental chamber for evaluation of gas-phase chemical mechanisms and secondary aerosol formation, *Atmos. Environ.* 39, 7768–7788, DOI 10.1016/j.atmosenv.2005.08.038. Two 90 m³ FEP reactors; **200 kW argon arc solar simulator or blacklights**; Fig. 2 compares arc and blacklight spectra "normalized to give the same NO2 photolysis rate". UCR progress reports give "light intensity with lamp at 80 % recommended maximum gives NO2 photolysis rate of 0.25–0.26 min⁻¹" (= 4.2–4.3e-3 s⁻¹). Shao et al. (2022, Table 2) list UCR j(NO2) = 3.83e-3 citing this paper. **The report's table of photolysis rates of other species relative to NO2 for blacklights vs. arc vs. sunlight (the item requested) could not be retrieved: intra.cert.ucr.edu / intra.engr.ucr.edu time out from this network and the Wayback Machine has no copy.** Carter's own characterisation of blacklights (referee comment on Leskinen et al., AMTD 7, C5010, 2014): "Figure 2 in Carter et al (2005) shows the spectrum of the type of blacklights usually used in environmental chambers … Although light sources should not have significant intensity below the 300 nm cutoff, it is important that they have a spectrum similar to sunlight in the important 300–400 nm region"; lamps with "no intensity in the UV region from 300 to about 340 nm … will result in abnormally low photolysis rates relative to atmospheric for other photolyses that are potentially important radical sources, such as photolysis of O3 forming O1D … and photolysis of aldehydes".
- **Warren, B., Song, C., Cocker, D. R.** (2008): Light intensity and light source influence on secondary organic aerosol formation for the m-xylene/NOx photooxidation system, *Environ. Sci. Technol.* 42, 5461–5466, DOI 10.1021/es702985n [DOI from memory — verify]: "80 115-W black lights and a variable voltage 200 kW argon arc lamp"; k1 varied 0.09–0.26 min⁻¹; at k1 = 0.09 min⁻¹ the arc lamp gave ≈20 % more SOA than blacklights.
- **Carter, W. P. L., Luo, D., Malkina, I. L.** (1997): Environmental chamber studies for development of an updated photochemical mechanism for VOC reactivity assessment, Final report to CARB contract 92-345 (https://ww2.arb.ca.gov/sites/default/files/classic/research/apr/past/92-345.pdf). Blacklight DTC and xenon-arc CTC chambers; Cl2/n-butane actinometry converted to k1 via the measured spectral distributions: **k1 = 2.89·k(Cl2) for the blacklight chamber and k1 = 5.15·k(Cl2) for the xenon arc chamber** (Eqs. V, VI) — a published spectrum-weighted ratio for two lamp types; blacklight runs had k1 ≈ 0.19 min⁻¹ (DTC); quartz-tube actinometry correction factor 1.66.
- **Carter, W. P. L.** (2010): Development of the SAPRC-07 chemical mechanism and updated ozone reactivity scales, report to CARB contracts 03-318/06-408/07-730 (https://ww2.arb.ca.gov/sites/default/files/classic/research/apr/past/03-318.pdf; also *Atmos. Environ.* 44, 5324, 2010). Table 19 lists the light source of every UCR chamber (EC/XTC/CTC: xenon arc; ITC/ETC/DTC: blacklights; EPA: argon arc or blacklights; TVA: blacklights + sunlamps). The companion spreadsheet **Photks.xls** (SAPRC file page) holds the photolysis-rate calculations incl. the blacklight and arc spectra — not retrievable from this network.

### 3.5 Carnegie Mellon (Donahue / Robinson / Presto / Hennigan)

- **Hennigan, C. J., et al.** (2011): Chemical and physical transformations of organic aerosol from the photo-oxidation of open biomass burning emissions in an environmental chamber, *Atmos. Chem. Phys.* 11, 7669–7686, DOI 10.5194/acp-11-7669-2011: "UV lights (GE model 10526 blacklights)" and natural sunlight; no j(NO2) given.
- ICARUS lists CMU CAPS as "9 m³ Teflon chamber illuminated with UVA radiation". Presto et al. 2005 (ES&T 39, 7046) and later CMU papers state j(NO2) only occasionally; **no value retrieved**. Springer guide Table 1.4: Carnegie Mellon, 16 m³, fluorescent tubes.

### 3.6 Georgia Tech Environmental Chamber (GTEC, Ng group)

- **Boyd, C. M., Sanchez, J., Xu, L., Eugene, A. J., Nah, T., Tuet, W. Y., Guzman, M. I., Ng, N. L.** (2015): Secondary organic aerosol formation from the β-pinene+NO3 system: effect of humidity and peroxy radical fate, *Atmos. Chem. Phys.* 15, 7497–7522, DOI 10.5194/acp-15-7497-2015: "The chambers are surrounded by black lights … output predominantly in the ultraviolet between 300 and 400 nm, with a maximum at 354 nm. The black lights are supplemented by natural sunshine fluorescent lights … wavelengths between 300 and 900 nm. **The jNO2 of the chamber facility is 0.28 min⁻¹ when all of the black lights are turned on**" (= 4.7e-3 s⁻¹).
- ICARUS light-flux characterisations #35/#37 (2020-07-28, GTEC Near/Far): black lights on left and right walls, UV-VIS spectrometer 170–840 nm averaged over four directions; "Photolysis frequency of NO2 at maximum lighting: 485.86 Hz" — almost certainly a unit slip for 4.86e-3 s⁻¹ (consistent with 0.28 min⁻¹) [R]; H2O2/O3 fields blank; spectrum file GT_UV_SL_20200728.csv.

### 3.7 ILMARI, University of Eastern Finland, Kuopio

- **Leskinen, A., Yli-Pirilä, P., Kuuspalo, K., Sippula, O., Jalava, P., Hirvonen, M.-R., Jokiniemi, J., Virtanen, A., Komppula, M., Lehtinen, K. E. J.** (2015): Characterization and testing of a new environmental chamber, *Atmos. Meas. Tech.* 8, 2267–2278, DOI 10.5194/amt-8-2267-2015. 29 m³ FEP. Two lamp sets: **Sylvania F36/T8 BLB 36 W ("365 nm lamps", 160 lamps, 18 switch groups)** and **Sylvania F 40W/350 BL ("350 nm lamps", 47 lamps)**. Total UV irradiance with all 160 lamps 29.7 W m⁻². NO2 actinometry (PSS, k = 1.8e-14): "The NO2 photolysis rate for the 365 nm lamps at 54, 108 and 160 lamps was 0.0037, 0.0058 and **0.0103 s⁻¹** (0.22, 0.35, 0.62 min⁻¹) … For the 350 nm lamps, the NO2 photolysis rate with all lamps on was **0.0022 s⁻¹** (0.132 min⁻¹)." "For the 365 and 350 nm lamps, there is no irradiation in the [<310 nm] region", so glyoxal etc. photolysis is suppressed; "in our first setup, the 365 nm blacklights would have been able to dissociate H2O2 into OH radicals only in a very small degree" (author reply). The paper's comparison list: PSI 27 m³ 0.12 min⁻¹ (Paulsen 2005), UCR 90 m³ 0.19 min⁻¹, Caltech 28 m³ 1.5 min⁻¹, PSI mobile 8.0e-3 s⁻¹.

### 3.8 Manchester Aerosol Chamber (MAC)

- **Shao, Y., Wang, Y., Du, M., Voliotis, A., Alfarra, M. R., O'Meara, S. P., Turner, S. F., McFiggans, G.** (2022): Characterisation of the Manchester Aerosol Chamber facility, *Atmos. Meas. Tech.* 15, 539–559, DOI 10.5194/amt-15-539-2022. 18 m³ FEP; "Two 6 kW arc xenon lamps (XBO 6000 W/HSLA OFR, Osram)" with 4 mm quartz filters "to filter out unwanted UV light" and "a bank of 112 halogen bulbs" for 290–800 nm. "j(NO2) as derived from our steady-state actinometry experiments was comparable … with that directly measured from the integrated absorption across the measured wavelengths (**2.25 ± 0.4 vs. 1.5e-3 s⁻¹**, respectively)"; total actinic flux in the MAC ≈ 1/3.5 of Manchester June midday clear sky (7e-3). Table 2 compiles j(NO2) of other chambers: KNU 2.83e-3 (Babar et al. 2016), GIG-CAS 8.17e-3 (Wang et al. 2014), UCR 3.83e-3 (Carter et al. 2005), AIOFM-CAS 3.17/3.50e-3 (Hu et al. 2014). No j(O1D) or other ratios.

### 3.9 CESAM, LISA Créteil

- **Wang, J., Doussin, J. F., Perrier, S., Perraudin, E., Katrib, Y., Pangui, E., Picquet-Varrault, B.** (2011): Design of a new multi-phase experimental simulation chamber for atmospheric photosmog, aerosol and cloud chemistry research, *Atmos. Meas. Tech.* 4, 2465–2494, DOI 10.5194/amt-4-2465-2011. 4.2 m³ stainless steel; **three 4 kW high-pressure xenon arc lamps** behind **9 mm Pyrex** filters (chosen so the spectrum matches ground-level sun; 2 mm or no filter mimic 45 / 60 km altitude). Spectroradiometer (LI-COR LI-1800) single point at the bottom centre: **j(NO2) = 4.5e-3, j(O1D) = 1.5e-5**; NO2 actinometry with wall-chemistry model: **j(NO2) = (4.2 ± 0.2)e-3** — "comparable to the 24 h-average of j(NO2) in the atmosphere". Abstract: j(NO2) = 4.2e-3 and j(O1D) = 1.4e-5. **Table 3, j(NO2)/j(O1D): 9 mm Pyrex 298 vs. solar ground level 260; 2 mm Pyrex 3.2 vs. 45 km 3.3; no filter 2.0 vs. 60 km 2.0.** Light flux varies by a factor 5 across the reactor (photodiode tomography). For modelling, "other photolysis frequencies used in the MCM were scaled from their atmospheric values using j(NO2) as a reference value" — i.e. MCM ratios were assumed, consistent with the measured j(O1D)/j(NO2) = 3.3e-3 [R].
- EUROCHAMP guide (§3.16) uses a CESAM N2-actinometry example: j(NO2) = (2.3 ± 0.2)e-3 (initial-slope) vs (2.8 ± 0.1)e-3 (full Eq. 2.3.5) — a different lamp configuration/date than 2011.

### 3.10 Sunlit chambers (reference only): EUPHORE, HELIOS, SAPHIR, UNC

- **Becker, K. H.** (ed.) (1996): The European Photoreactor EUPHORE — Final report, EC contract EV5V-CT92-0059, Wuppertal. Two 200 m³ FEP hemispheres, Valencia; photolysis frequencies from spectroradiometry/filter radiometers inside the chamber.
- **Zádor, J., Turányi, T., Wirtz, K., Pilling, M. J.** (2006): Measurement and investigation of chamber radical sources in the European Photoreactor (EUPHORE), *J. Atmos. Chem.* 55, 147–166, DOI 10.1007/s10874-006-9033-y. Paywalled; abstract gives HONO and HCHO wall-source parameterisations scaled to j(NO2) (W(HCHO) = c·j(NO2)·exp(−5686 K/T), c = 3.1e17 cm⁻³), no j tables retrieved. Bloss et al. (2005, ACP 5, 623) and Heitlin/Wirtz EUPHORE studies scale all photolysis rates with measured j(NO2)/calculated j(NO2) — i.e. MCM ratios again.
- **Ren, Y., Grosselin, B., Daële, V., Mellouki, A.** (2017): Investigation of the reaction of ozone with isoprene, methacrolein and methyl vinyl ketone using the HELIOS chamber, *Faraday Discuss.* 200, 289–311, DOI 10.1039/C7FD00014F. HELIOS: 90 m³ FEP hemisphere, Orléans, >90 % solar transmission; an AMT characterisation paper was not found.
- SAPHIR: §2.2.

### 3.11 PSI PACS-C3 and the PSI mobile chamber

- **Paulsen, D., Dommen, J., Kalberer, M., Prévôt, A. S. H., Richter, R., Sax, M., Steinbacher, M., Weingartner, E., Baltensperger, U.** (2005): Secondary organic aerosol formation by irradiation of 1,3,5-trimethylbenzene–NOx–H2O in a new reaction chamber for atmospheric chemistry and physics, *Environ. Sci. Technol.* 39, 2668–2678, DOI 10.1021/es0489137. 27 m³ FEP; four 4 kW xenon arc lamps plus blacklights; j(NO2) = 0.12 min⁻¹ = 2.0e-3 s⁻¹ (as quoted by Leskinen 2015 and Vallon 2022) [?].
- **Platt, S. M., El Haddad, I., Zardini, A. A., Clairotte, M., Astorga, C., Wolf, R., Slowik, J. G., Temime-Roussel, B., Marchand, N., Ježek, I., Drinovec, L., Močnik, G., Möhler, O., Richter, R., Barmet, P., Bianchi, F., Baltensperger, U., Prévôt, A. S. H.** (2013): Secondary organic aerosol formation from gasoline vehicle emissions in a new mobile environmental reaction chamber, *Atmos. Chem. Phys.* 13, 9141–9158, DOI 10.5194/acp-13-9141-2013. 9 m³ (later 7 m³) FEP, **40 × 100 W Philips Cleo Performance solarium lamps**, emission "peak at 368 nm … narrow window between 310 and 430 nm" (Fig. 6). NO2 actinometry at 25 °C: **j(NO2) = (8.0 ± 0.7)e-3 s⁻¹** (5 % sampling-line correction). Table 5 "Predicted gas phase photolysis rates" from the measured spectrum scaled to j(NO2) with IUPAC σ and φ (Eq. 4), 25 °C / −7 °C, in 1e-5 s⁻¹:

  | process | 25 °C | −7 °C | ratio to j(NO2) at 25 °C [R] |
  |---|---|---|---|
  | NO2 → NO + O(3P) | 797 | 259 | 1 |
  | O3 → O(1D) | 43.1 | 17.8 | **5.4e-2 (!)** |
  | HONO → OH + NO | 214 | 68.4 | 0.27 |
  | NO3 → NO + O2 | 0.00 | 0.00 | 0 |
  | NO3 → NO2 + O(3P) | 444 | 135 | 0.56 |
  | HCHO → HCO + H | 1.50 | 0.53 | 1.9e-3 |
  | HCHO → CO + H2 | 3.74 | 1.35 | 4.7e-3 |
  | CH3CHO → products | 0.11 | — | 1.4e-4 |
  | CH3OOH → products | 0.43 | — | 5.4e-4 |
  | (CHO)2 → products | 0.71 | — | 8.9e-4 |
  | MVK → products | 0.34 | — | 4.3e-4 |
  | CH3C(O)CH3 → products | 0.01 | — | 1.3e-5 |
  | CH3ONO2 → products | 0.29 | — | 3.6e-4 |
  | PAN → products | 0.04 | — | 5e-5 |

  Lamp output at −7 °C is 0.28 of that at 25 °C, so all j scale by 0.30–0.41 (Eq. 5 gives the temperature fit).
  **Caveat on the O3 row:** a 5.4 % O1D/NO2 ratio is 13× the clear-sky value (4.2e-3 at SZA 0) and physically implausible for a lamp with no emission below 310 nm, where σ(O3)·φ(O1D) is ≪ σ(NO2); the other rows (HCHO, CH3OOH, acetone, CH3CHO all *below* solar ratios) show the expected UV-poor behaviour. The O3 entry is likely a tabulation or quantum-yield error in the source; do not use it. The remaining rows are the best published spectrum-based ratio set for a UVA fluorescent lamp and are adopted (with adjustment for the 350 nm phosphor) in §5.

### 3.12 AIDA LED solar simulator (KIT)

- **Vallon, M., Gao, L., Jiang, F., Krumm, B., Nadolny, J., Song, J., Leisner, T., Saathoff, H.** (2022): LED-based solar simulator to study photochemistry over a wide temperature range in the large simulation chamber AIDA, *Atmos. Meas. Tech.* 15, 1795–1810, DOI 10.5194/amt-15-1795-2022. 3000 LEDs, 11 wavelengths 305–528 nm (FWHM 10–15 nm) on top of the 85 m³ chamber, 1400 W optical; NO2 actinometry (k = 2.07e-12·exp(−1400/T)): **j(NO2) = (1.58 ± 0.21)e-3 s⁻¹**, 213–295 K; ≈1/3 of Karlsruhe clear summer noon; vertical gradient (centre ≈ 45 % of top). 2,3-pentanedione j = (1.03 ± 0.15)e-4 (400–460 nm). No j(O1D). Comparison values quoted: PSI 27 m³ 2.0e-3 (Xe + blacklights), CESAM 4.2e-3. Related LED chambers: CLOUD (CERN), SAPHIR-PLUS (FZJ), UBC ATMOX (Tsai et al., EGUsphere 2025-3041).

### 3.13 HIRAC, University of Leeds

- **Glowacki, D. R., Goddard, A., Hemavibool, K., Malkin, T. L., Commane, R., Anderson, F., Bloss, W. J., Heard, D. E., Ingham, T., Pilling, M. J., Seakins, P. W.** (2007): Design of and initial results from a Highly Instrumented Reactor for Atmospheric Chemistry (HIRAC), *Atmos. Chem. Phys.* 7, 5371–5390, DOI 10.5194/acp-7-5371-2007. 2.25 m³ stainless steel; **24 × Philips TLK40W/05 actinic UV blacklamps** in 8 quartz tubes (output spectrum Fig. 2); NO2 PSS actinometry with 100–200 ppb NO2: **j(NO2) = (1.93 ± 0.10)e-3 s⁻¹**; analytical lamp-radiation-field model (ray tracing, reflectivity 0.55). Lamps are interchangeable: Philips TL40W/12 RS (UV-B, 275–390 nm) are used when O3 → O(1D) or HCHO photolysis is wanted (Winiberg et al. 2015, AMT 8, 523–540, DOI 10.5194/amt-8-523-2015), and 254 nm germicidal lamps are available — but no j(O1D) actinometry value was found in these two papers.

### 3.14 GIG-CAS, Guangzhou (365-nm-type blacklamps)

- **Wang, X., Liu, T., Bernard, F., Ding, X., Wen, S., Zhang, Y., Zhang, Z., He, Q., Lü, S., Chen, J., Saunders, S., Yu, J.** (2014): Design and characterization of a smog chamber for studying gas-phase chemical mechanisms and aerosol formation, *Atmos. Meas. Tech.* 7, 301–313, DOI 10.5194/amt-7-301-2014. 30 m³ FEP; **135 black lamps (1.2 m, 60 W)** "irradiation over the range from 340 to 400 nm with peak intensity at 369 nm"; NO2 PSS actinometry: **j(NO2) = 0.49 ± 0.01 min⁻¹ (8.2e-3 s⁻¹)** at full intensity vs. 0.58 min⁻¹ TUV Guangzhou 22 June noon; four lamp groups give 0–0.49 min⁻¹; a Metcon j(NO2) filter radiometer monitors drift. HONO wall source 0.005·j(NO2) ppb s⁻¹.

### 3.15 NIES, Tsukuba (350 nm blacklights) — the one published blacklight j(H2O2)/j(NO2)

- **Sato, K., et al.** (2011): author reply, *Atmos. Chem. Phys. Discuss.* 11, C5678–C5682 (paper: Sato, K., Nakao, S., Clark, C. H., Qi, L., Cocker, D. R., *Atmos. Chem. Phys.* 11, 7301–7317, DOI 10.5194/acp-11-7301-2011): "Black light spectrum is peaked at 350 nm with 40 nm FWHM. The H2O2 photolysis rate when 276 lights are used was evaluated to be **6.4e-4 min⁻¹** from the NO2 photolysis rate, the black light spectrum, and the H2O2 and NO2 absorption spectra (DeMore et al., 1997)"; standard conditions j(NO2) = 0.12 min⁻¹. → **j(H2O2)/j(NO2) = 5.3e-3 [S/R]**, 6× the clear-sky ratio, because the 350 nm phosphor has no 400–420 nm output where NO2 (but not H2O2) absorbs. The sensitivity to the <320 nm tail is large (compare §3.3).

### 3.16 EUROCHAMP-2020 guide and the ICARUS archive

- **Doussin, J.-F., Fuchs, H., Kiendler-Scharr, A., Seakins, P., Wenger, J.** (eds.) (2023): *A Practical Guide to Atmospheric Simulation Chambers*, Springer, open access, DOI 10.1007/978-3-031-22277-1. Chapter 2, **Alfarra, M. R., Camredon, M., Cazaunau, M., Doussin, J.-F., Fuchs, H., Jorga, S., McFiggans, G., Newland, M., Pandis, S. N., Rickard, A. R., Saathoff, H.**: Physical and chemical characterization of the chamber, pp. 73–111. Sect. 2.3: j(X) = ∫σφF dλ with σ, φ from Atkinson et al. 2004 / Burkholder et al. 2020 / Keller-Rudek et al. 2013; chemical actinometry "has mostly been confined to the determination of photolysis frequencies j(NO2) and j(O1D)"; protocols for j(NO2) actinometry in air (PSS, Leighton ratio; Eq. 2.3.3; SAPHIR example) and in N2 (NO2 decay; Eq. 2.3.4/2.3.5; CESAM example 2.3–2.8e-3 s⁻¹); recommendation to combine spectroradiometry (shape) with actinometry (absolute scale) and to record lamp spectra with the data. Chapter 1 Table 1.4 inventories ≈50 chambers with their light source (fluorescent tubes dominate; Xe arc: CESAM, MAC, PACS-C3, NIES, UCR; LED: AIDA, CLOUD, SAPHIR-PLUS; sun: EUPHORE, HELIOS, SAPHIR, UNC, CRAES). **The guide contains no table of j(X)/j(NO2) ratios per lamp type.** Doussin's closing chapter: blacklights remain "the most common light source … even if the atmospheric relevance of their emission spectrum can be questioned".
- **ICARUS** (https://icarus.ucdavis.edu): chamber records give lamp models and link "light flux characterization" entries with downloadable spectra (CSV irradiance or flux vs. wavelength). Entries checked: Jimenez 1W (#38–#42, no j shown), Seinfeld (#18–#21, 320 × 350 nm lamps, j blank), GTEC (#35/#37, "485.86 Hz"), Hildebrandt Ruiz LHR (#24, 240 UVA blacklights, Apogee PS-200 irradiance, j blank), Tran big bag (#27, 10–20 % top/bottom inhomogeneity). The archive is the best place to get *actual* lamp spectra for the US chambers; j values are rarely filled in.

### 3.17 Chambers for which nothing usable was found

- **CHARME, Dunkerque (ULCO):** no characterisation paper located (searched standard + extended); only a thesis record ("characterization of the new atmospheric simulation chamber CHARME … γ-terpinene ozonolysis"). The "Lamkaddam et al. 2019" lead is LISA/CESAM work (Lamkaddam did the CESAM n-dodecane study), not CHARME. No Xe-arc j ratios collected.
- **QUAREC, Wuppertal:** Barnes, I., Becker, K. H., Mihalopoulos, N. (1994), *J. Atmos. Chem.* 18, 267–289 — 1080 L quartz reactor with 32 Philips TL05 40 W superactinic (320–480 nm, max 360 nm) and 32 Philips TUV 40 W (254 nm) lamps; j(NO2) not found. A Wuppertal flow-reactor paper quotes J(NO2) = 0.027 s⁻¹ for 6 × Philips TL/05 lamps [?].
- **UT Austin (Hildebrandt Ruiz):** ICARUS LHR chamber 8.92 m³, 240 UVA blacklights + visible fluorescent; no j value published that I could find.
- **UC Davis (Kleeman / Cappa / Clark):** not found in ICARUS listings under those names; nothing collected.
- **NCAR:** no in-house chamber paper by "Zhang et al. 2015 AMT" exists to my knowledge (the 2015 Zhang AMT paper on highly reactive VOCs is a measurement paper). A search snippet mentions an "NCAR chamber j(NO2) ≈ 1.27e-3" [?] without a traceable source.
- **KNU (Korea), AIOFM-CAS, Purdue 10 m³ (128 lamps, 1.4e-3 s⁻¹) [?]**: only via Shao's Table 2 or search snippets; primary texts not retrieved.
- **CSIRO (Hynes et al., Environ. Chem. 2018, DOI 10.1071/EN18005):** compares two blacklight types (350 nm vs "blue" BLB) and reports j(NO2) drift over ≈150 experiments for the 350 nm lamps; page returned 404, values not collected.

---

## 4. How to compute J<n> for an arbitrary lamp spectrum

### 4.1 The integral and its practical form

`j(X) = ∫ σ_X(λ,T) · φ_X(λ,T,p) · F(λ) dλ`, with F the **spectral actinic flux density** (photons cm⁻² s⁻¹ nm⁻¹, integrated over 4π sr), σ in cm², φ dimensionless.

Practical chamber recipe (Platt et al. 2013 Eq. 4; Alfarra et al. 2023 Sect. 2.3; Bohn & Zilken 2005):

1. Measure or obtain the **relative** lamp spectrum I(λ) (spectroradiometer or manufacturer curve). Irradiance (W m⁻² nm⁻¹, cosine-weighted) is not actinic flux; for an all-around lamp array convert either by a geometric factor (Wang 2011 combine up- and down-looking spectra with pseudo-zenith angles) or simply treat I(λ) as relative and fix the scale in step 3.
2. Convert to photons: F(λ) ∝ I(λ)·λ (λ in nm, I in W). Resample σ·φ onto the lamp grid (0.5–1 nm; Hg lines need ≤0.1 nm or line-by-line treatment).
3. Fix the absolute scale with NO2 actinometry: `j(X) = j(NO2)_act · ∫σ_Xφ_X F dλ / ∫σ_NO2φ_NO2 F dλ`. For OFRs with Hg lamps, use the measured I254 (O3 or SO2/CO tracer, Lambe 2019; Rowe 2020) rather than NO2 actinometry (see 4.3).
4. Apply temperature: lamp output (fluorescent tubes lose ≈70 % between 25 and −7 °C, Platt Eq. 5) and σ, φ temperature dependence (O3→O1D φ tail; acetone φ pressure/T).
5. Record spatial inhomogeneity (factor 5 in CESAM, 10–20 % top/bottom in bag chambers) and lamp aging (CSIRO drift; GIG-CAS radiometer).

### 4.2 Data sources

- **Cross sections / quantum yields:** JPL Publication 19-5, Burkholder, J. B., et al. (2020), *Chemical Kinetics and Photochemical Data for Use in Atmospheric Studies*, https://jpldataeval.jpl.nasa.gov; IUPAC Task Group (Atkinson et al. 2004 ACP 4, 1461; 2006 ACP 6, 3625; Mellouki et al. 2021), https://iupac.aeris-data.fr; **MPI-Mainz UV/VIS Spectral Atlas** (Keller-Rudek, H., Moortgat, G. K., Sander, R., Sörensen, R., 2013, ESSD 5, 365, DOI 10.5194/essd-5-365-2013), https://uv-vis-spectral-atlas-mainz.org — files such as `NO2_JPL-2010(2011)_294K_240-662.5nm(rec).txt`, `H2O2_JPL-2010(2011)_298K_190-350nm(rec).txt`, `HONO_JPL-2010(2011)_298K_184-396nm(rec).txt`, `O3_Brion(1985)_298K_253.65-334.15nm(rec).txt`, `CH3OOH_JPL-2010(2011)_298K_210-405nm(rec).txt` (the site serves the data only through its own pages; deep links returned the index from this network). The MCM photolysis page links the σ/φ spreadsheets used for each J.
- **Solar reference:** TUV (Madronich, S., Flocke, S., 1999, in *Handbook of Environmental Chemistry* 2L, 1–26; Quick TUV calculator https://www.acom.ucar.edu/Models/TUV/Interactive_TUV/) or the MCM l,m,n fits (§2.1).
- **Lamp spectra:** Lambe et al. 2019 Fig. S1 (GPH436T5L 254 nm, F436T5/BL-350, F436T5/BLC-369 manufacturer curves); Rowe et al. 2020 Fig. S2 (relative 185/254 output); Li et al. 2015 (BHK 82-9304-03 vs voltage); Carter et al. 2005 Fig. 2 (UCR blacklight and Ar arc); Leskinen et al. 2015 Fig. 2 and reply Fig. R1 (Sylvania 350BL vs BLB); Platt et al. 2013 Fig. 6 (Cleo Performance); Glowacki et al. 2007 Fig. 2 (TLK40W/05); Wang et al. 2011 Fig. 4 (Xe + 0/2/9 mm Pyrex vs TUV); Vallon et al. 2022 (LED set); Boyd et al. 2015 (GTEC); ICARUS CSV files (GTEC 2020, Seinfeld 2013, LHR 2018); PAM wiki (OFR lamps). Manufacturer data sheets: Light Sources Inc. (GPH series), Sylvania/Osram (F40/350BL, F36T8/BLB), Philips (TL-K 40W/05, TL 40W/12, TUV), Q-Lab (UVA-340 — a "sunlight-matching" 295–365 nm fluorescent tube, not covered by any chamber paper found).

### 4.3 Low-pressure Hg lamps: numbers for the 185/254 nm case [R]

Cross sections at 253.65 nm (298 K): σ(O3) = 1.15e-17 (Hartley peak; JPL/Brion), φ(O1D) = 0.90; σ(NO2) = 1.05e-20 (IUPAC, as used by Lambe et al. 2019); σ(HNO3) = 1.18e-20 (JPL-2010, Mainz atlas value at 254 nm — verified); the following are **recalled JPL/IUPAC values, verify before use:** σ(H2O2) ≈ 7.0e-20 (φ = 1); σ(CH3OOH) ≈ 3.3e-20 (φ = 1); σ(acetone) ≈ 3.0e-20 (φ ≈ 1 at 254 nm, 1 atm); σ(CH3CHO) ≈ 1.2e-20 (φ_rad ≈ 0.3–0.5); σ(HCHO) ≈ 0.5–1e-20 (φ_r ≈ 0.3, φ_m ≈ 0.5 at 250–280 nm); σ(HONO,254) not verified (weak; HONO's strong bands are <220 nm and 300–400 nm); σ(NO3,254) ≈ 0. At 184.9 nm: σ(H2O) = 7.1e-20 (φ = 1); σ(O2) ≈ 1e-20 effective (Schumann–Runge/Herzberg, path and p dependent); σ(N2O) = 1.4e-19; σ(O3) ≈ 6e-19.

Per **F254 = 1e15 photons cm⁻² s⁻¹** (mid-range of Peng et al.; Aerodyne PAM full power ≈ 3.5–6.5e15):

| process | j (s⁻¹) | j/j(NO2) with j(NO2) = σ254·F254 = 1.05e-5 | j/j(NO2) with actinometric j(NO2) ≈ 13 × 1.05e-5 = 1.4e-4 (Lambe 2019 minor lines) |
|---|---|---|---|
| O3 → O(1D) | 1.04e-2 | ≈ 990 | ≈ 75 |
| O3 → O(3P) (φ = 0.1) | 1.2e-3 | ≈ 110 | ≈ 8 |
| H2O2 | 7.0e-5 | 6.7 | 0.5 |
| CH3OOH | 3.3e-5 | 3.1 | 0.24 |
| acetone | ≈ 3e-5 | ≈ 3 | ≈ 0.2 |
| HNO3 | 1.2e-5 | 1.1 | 0.09 |
| HCHO (total) | ≈ 0.5–1e-5 | ≈ 0.5–1 | ≈ 0.05 |
| CH3CHO | ≈ 5e-6 | ≈ 0.5 | ≈ 0.04 |
| NO2 | 1.05e-5 | 1 | (1) |
| NO3 | ≈ 0 | 0 | 0 |

For OFR185 add j(H2O) = 7.1e-20·F185 (= 7.1e-7 s⁻¹ at F185 = 1e13) and j(O2) ≈ 1e-20·F185; F185/F254 is lamp specific (≈0.01 BHK with Li et al. calibration; ≈0.1 GPH436T5VH per Rowe et al.).
The take-away for a profile: **in a 254 nm reactor the j(O1D)/j(NO2) ratio is 10⁴–10⁵ × solar and depends on how j(NO2) is defined**; JUICE should carry F254 (and F185) explicitly and compute j from σ·F rather than scale from j(NO2).

### 4.4 Using the MCM table inside JUICE

For a sunlit reference, `ratio_X(χ) = J_X(χ)/J4(χ)` from the l,m,n fits; for any lamp, store either (a) F(λ) and compute, or (b) a ratio set like §5. Do **not** scale lamp chambers by MCM ratios unless the source is a Pyrex-filtered Xe/Ar arc with demonstrated j(O1D)/j(NO2) ≈ 3e-3 (CESAM) — this is exactly what CESAM and EUPHORE modellers do, and it is wrong for blacklights by one to two orders of magnitude for O1D, H2O2 and carbonyls (and in the opposite direction for HONO and NO3).

---

## 5. Recommended default ratio sets for JSON profiles

Values are j(X)/j(NO2) unless stated. "u" is the multiplicative 1σ-ish uncertainty factor I would attach
(e.g. u = 3 means ×/÷ 3). Entries marked `est` are my estimates from the overlap arguments in §3/§4, not literature values.

```json
{
  "clear_sky_sza30": {
    "source": "MCM l,m,n (Saunders et al. 2003) at SZA 30 deg; consistent with SAPHIR 28 Jul 2002 (Bohn & Zilken 2005)",
    "jNO2_abs_s-1": 8.3e-3, "jNO2_u": 1.15,
    "O1D": 3.3e-3, "O3P": 5.0e-2, "HONO": 0.22, "H2O2": 8.2e-4,
    "HCHO_r": 3.3e-3, "HCHO_m": 5.3e-3, "CH3CHO": 4.6e-4, "acetone": 5.6e-5,
    "CH3OOH": 6.1e-4, "NO3_a": 2.6, "NO3_b": 17.9, "HNO3": 6.6e-5,
    "ratio_u": 1.2,
    "note": "inside SAPHIR multiply jNO2 by 0.65-0.86 and O1D ratio by ~0.87 (0.68/0.78)"
  },
  "generic_blacklight_350": {
    "source": "Sylvania F40/350BL, GE F40BL, TL-K 40W/05, UVA-350 class; 310-400 nm, peak 350-355 nm",
    "jNO2_abs_s-1_range": [1.4e-3, 1.0e-2], "jNO2_typical": 4.5e-3,
    "O1D": 1.0e-3, "O1D_u": 3, "O1D_note": "est; no actinometric value in the literature; governed by the 300-315 nm phosphor tail - measure or convolve the real spectrum",
    "HONO": 0.40, "HONO_u": 1.5,
    "H2O2": 1.0e-3, "H2O2_u": 4, "H2O2_note": "NIES 350/40 nm lamps 5.3e-3 (Sato 2011); Caltech H2O2 experiments imply ~3e-4",
    "HCHO_r": 3.0e-3, "HCHO_r_u": 1.7, "HCHO_m": 5.0e-3, "HCHO_m_u": 1.5,
    "CH3CHO": 3.0e-4, "CH3CHO_u": 2, "acetone": 3.0e-5, "acetone_u": 3,
    "CH3OOH": 1.0e-3, "CH3OOH_u": 3,
    "NO3_a": 0.0, "NO3_b": 0.05, "NO3_u": 10, "NO3_note": "no emission beyond ~410 nm; NO3 photolysis negligible",
    "HNO3": 5e-5, "HNO3_u": 5
  },
  "generic_blacklight_365": {
    "source": "Sylvania F36/T8 BLB, GE F40BLB, Aerodyne F436T5/BLC-369, GIG-CAS lamps; 340-400 nm, peak 365-369 nm; ILMARI 1.03e-2 (160 lamps), GIG-CAS 8.2e-3 (135 lamps)",
    "jNO2_abs_s-1_range": [2e-3, 1.0e-2],
    "O1D": 1.0e-4, "O1D_u": 10, "HONO": 0.30, "HONO_u": 1.5,
    "H2O2": 1.5e-4, "H2O2_u": 5, "HCHO_r": 1.0e-3, "HCHO_r_u": 2, "HCHO_m": 3.0e-3, "HCHO_m_u": 2,
    "CH3CHO": 1.0e-4, "CH3CHO_u": 3, "acetone": 1.0e-5, "acetone_u": 5, "CH3OOH": 2.0e-4, "CH3OOH_u": 4,
    "NO3_a": 0.0, "NO3_b": 0.05, "NO3_u": 10,
    "note": "Platt et al. 2013 (368 nm Cleo lamp, 310-430 nm) measured set: HONO 0.27, HCHO_r 1.9e-3, HCHO_m 4.7e-3, CH3CHO 1.4e-4, CH3OOH 5.4e-4, acetone 1.3e-5, NO3_b 0.56 (that lamp reaches 430 nm); its O1D entry (5.4e-2) is not credible"
  },
  "generic_hg254": {
    "source": "OFR254: GPH436T5L / BHK 82-9304-03 with 185-blocking quartz; Peng et al. 2015/2016, Lambe et al. 2019",
    "F254_photons_cm-2_s-1_range": [4.2e13, 8.5e15], "F254_full_power_2_lamps": 3.5e15, "F254_full_power_4_lamps": 6.5e15,
    "mode": "compute j = sigma(254 nm) * F254; do not scale from jNO2",
    "sigma_cm2_254nm": {"O3": 1.15e-17, "O3_phiO1D": 0.9, "NO2": 1.05e-20, "H2O2": 7.0e-20, "CH3OOH": 3.3e-20, "acetone": 3.0e-20, "HNO3": 1.18e-20, "CH3CHO": 1.2e-20, "HCHO": 0.7e-20, "HONO": null},
    "jNO2_minor_line_factor": 13, "jNO2_minor_line_u": 2,
    "O1D_over_jNO2_254only": 1.0e3, "O1D_over_jNO2_actinometric": 75, "ratio_u": 2
  },
  "generic_hg185_254": {
    "source": "OFR185: GPH436T5VH (F185/F254 ~0.10, Rowe 2020) or BHK 82-9304-03 (F185/F254 ~0.01, Li 2015)",
    "F185_over_F254": 0.01, "F185_over_F254_range": [0.001, 0.1],
    "sigma_cm2_185nm": {"H2O": 7.1e-20, "O2_effective": 1.0e-20, "N2O": 1.4e-19, "O3": 6e-19},
    "inherit": "generic_hg254"
  },
  "generic_xe_arc_pyrex": {
    "source": "CESAM 3x4 kW Xe + 9 mm Pyrex (Wang et al. 2011); MAC 2x6 kW XBO + halogen (Shao et al. 2022); UCR 200 kW Ar arc (Carter et al. 2005)",
    "jNO2_abs_s-1_range": [2.0e-3, 4.5e-3],
    "O1D": 3.3e-3, "O1D_u": 1.3,
    "others": "use clear_sky_sza30 ratios with u = 1.5 (HONO, HCHO, H2O2, CH3OOH, CH3CHO, acetone); NO3 ratios depend on visible filtering/halogen fill - u = 2",
    "note": "unfiltered or 2 mm Pyrex Xe gives O1D ratios of 0.5 and 0.3 (stratospheric); never use an unfiltered arc profile for the boundary layer"
  },
  "generic_led_solar": {
    "source": "AIDA LED simulator (Vallon et al. 2022): jNO2 1.58e-3; CLOUD, SAPHIR-PLUS, UBC ATMOX similar concept",
    "jNO2_abs_s-1": 1.6e-3, "jNO2_u": 1.3,
    "O1D": null, "note": "ratios depend entirely on the 305-320 nm LED share; must be computed from the installed LED spectrum (manufacturer curves + measured weights)"
  }
}
```

Rules of thumb behind these defaults:

- A blacklight shifts weight from the NO2-dominated 380–420 nm band to 330–370 nm. Relative to sunlight this **raises** j(HONO)/j(NO2) (×1.5–2), keeps j(HCHO)/j(NO2) within a factor ≈1.5, **lowers** j(CH3CHO), j(acetone), j(CH3OOH)/j(NO2) by 2–5 (368 nm lamp: by 4–5), and makes j(NO3) ≈ 0. j(O1D) and j(H2O2) are set by the <320 nm tail and vary by more than an order of magnitude between a 350 nm and a 365 nm phosphor — this is the single largest uncertainty and should be flagged in any JUICE profile that lacks a measured spectrum.
- A Pyrex-filtered Xe/Ar arc reproduces solar ratios to within the measurement precision of j(O1D)/j(NO2) (CESAM: 298 vs 260); scale MCM.
- A 254 nm Hg reactor cannot be described by ratios to j(NO2) at all; carry photon fluxes.

---

## 6. Gaps (what was not found)

- No chamber paper reports an **actinometric j(O1D)** for blacklights; all blacklight O1D ratios here are computed estimates.
- Carter et al. 2005 / SAPRC Photks.xls ratio tables for UCR blacklights and the Ar arc: not retrievable (UCR hosts unreachable, no archive copy).
- CU Boulder chamber j(NO2), Caltech current-chamber j(NO2), CMU, UT Austin, UC Davis, CHARME, QUAREC: no published values located.
- Zádor et al. 2006 EUPHORE j tables: paywalled.
- MPI-Mainz atlas deep links served only the index page from this network; H2O2, CH3OOH, acetone, CH3CHO, HCHO, HONO cross sections at 254 nm are therefore recalled values marked for verification.
