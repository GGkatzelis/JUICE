# How chamber studies specified photolysis frequencies for MCM / SAPRC — historical practice

Compiled 2026-11-02 for JUICE chamber profiles (`chambers/<id>.chamber.json`). Companion to
`docs/lamp_photolysis_literature.md` (which collects the *numbers*; this file collects the *methods*).

Question answered here: when people modelled chamber experiments with the Master Chemical Mechanism
(MCM) or SAPRC, **how did they fix j(NO2), how did they get every other J, did they use the MCM
`J = l·cos^m χ·exp(−n·sec χ)` parameterisation at some effective zenith angle, and where did the chamber
auxiliary mechanism come from?**

Labels used throughout:

| Label | Meaning |
|---|---|
| **[V]** | verbatim quote checked against the full text (PDF extracted with `pdftotext`) |
| **[S]** | quote taken from a search-engine snippet or a secondary source; not checked against the full text |
| **[?]** | not retrieved / inferred; treat as unverified |

Method classes (used in the table and in §2):

| Class | Definition |
|---|---|
| **A** | lamp spectrum × σ(λ)·φ(λ), integrated, with the absolute scale anchored on NO2 actinometry (or an equivalent absolute measurement) |
| **B** | solar (MCM / radiative-transfer) parameterisation for the chamber's location, date and real SZA, with all J scaled by j(NO2)_measured / j(NO2)_calculated (sometimes with a second, species-dependent factor) |
| **C** | ratios j(X)/j(NO2) imported from a reference chamber or a canonical lamp spectrum (not measured in the chamber in question) |
| **D** | J measured directly in the chamber (spectroradiometer or actinometry) for the species in question |
| **E** | J (or an effective quantum yield) fitted to the experiment |

Important nuance found while reading: almost every real study is a *combination* — typically **D for
j(NO2) (and sometimes a few more), then A or B for everything else**. The table therefore has two
separate columns.

---

## 1. Summary table

| Study | Chamber, light source | Mechanism, solver | How j(NO2) was fixed | How all other J were specified | Chamber auxiliary mechanism: source | DOI / source |
|---|---|---|---|---|---|---|
| Carter 2000, *Documentation of the SAPRC-99 chemical mechanism*, §"Light Characterization" | SAPRC/CE-CERT indoor chambers ETC, DTC (blacklights), ITC, EC, XTC, CTC (Xe arc), OTC (outdoor) | SAPRC-99 | **D**: periodic quartz-tube NO2 actinometry; Cl2–n-butane actinometry for Xe-arc chambers | **A**: "NO2 photolysis rate … and the relative spectral distribution of the light source"; blacklight spectrum assumed constant = the Carter et al. (1995d) recommended spectrum (**C-flavoured A**); Xe-arc spectrum measured around the time of each run; outdoor: solar model fitted to UV + TSR radiometers | Carter's own: `Walls + hν → HONO` (RN), RS, NO2W, N25I/N25S, O3W; parameters per "characterization set" from CO–NOx / n-butane–NOx / pure-air runs | Report to CARB, contracts 92-329 & 95-308 (UCR, May 2000) [V] |
| Carter et al. 2005 AE 39, 7768 (report `chamrpt.pdf`, 2004) | UCR EPA chamber, 2 × 90 m³; Ar arc solar simulator and bank of 115 W blacklights | SAPRC-99 | **D**: in-reactor quartz-tube NO2 actinometry (0.26 min⁻¹ arc; blacklight trend from QSL sensor calibrated against in-chamber actinometry) | **A**: "calculated using the measured NO2 photolysis rates, the relative measured spectral distributions for the light sources, and the absorption cross sections and quantum yields" (LI-COR LI-1800 spectroradiometer) | Carter's own: HONO off-gassing parameter RN (ratio to j(NO2)), NOx off-gassing, HCHO off-gassing; compared with SAPHIR LOPAP data of Rohrer et al. | 10.1016/j.atmosenv.2005.08.040; `intra.cert.ucr.edu/~carter/pubs/chamrpt.pdf` [V] |
| Carter 2010, SAPRC-07 report / AE 44, 5324 | all UCR chambers above | SAPRC-07 | **D** (as above) | **A**: "all photolysis rates are normalized to measured NO2 photolysis rates" | Carter's own, RN-I parameters re-fitted per chamber (e.g. EPA 0.0075 ppb, ITC 0.0165, DTC 0.0115, XTC 0.0055) | 10.1016/j.atmosenv.2010.01.026 (paper); `intra.cert.ucr.edu/~carter/SAPRC/saprc07.pdf` (report) [V] |
| Carter, Luo, Malkina, Pierce 1995 (ARB 92-345) | DTC (blacklight) and CTC (Xe arc) | SAPRC-90/-93 | **D**: quartz-tube NO2 actinometry (blacklight); Cl2–n-butane actinometry (Xe arc); NO/NO2/O3 steady-state judged unreliable | **A** (light-source-specific photolysis rates; "differences between blacklights and solar lighting can be accounted for in model simulations") | Carter's own | ARB report 92-345 (`ww2.arb.ca.gov/.../92-345.pdf`) [V] |
| Pinho, Pio, Jenkin 2005 AE 39, 1303 (isoprene, MACR, MVK) | SAPRC ETC/DTC/ITC/EC/XTC/CTC datasets; blacklights (DTC/ETC/ITC), Xe arc (XTC/CTC) | MCM v3 (FACSIMILE) | **D** (Carter's actinometry; snippet quotes 0.21 min⁻¹ for one chamber) | **A**, Carter's method: "calculated using the NO2 photolysis rates (experimentally determined …), the relative spectral distributions for the light sources, the absorption cross sections and quantum yields" [S] | Carter's: "chamber-dependent auxiliary mechanism" (Carter et al. 1995a), "updated by Pinho et al. (2005)" (per Metzger 2008 [V]; Pinho EGU 2006 abstract [V]) | 10.1016/j.atmosenv.2004.11.014 (full text not retrieved) |
| Pinho, Pio, Carter, Jenkin 2006 J. Atmos. Chem. 55, 55 (alkenes); Pinho, Pio, Green, Jenkin 2007 J. Atmos. Chem. 57, 171 (α-, β-pinene) | same SAPRC chambers | MCM v3 / v3.1 | as 2005 [?] | as 2005 [?] (Carter is co-author of the 2006 paper) | Carter's [?] | 10.1007/s10874-006-9025-y; 10.1007/s10874-007-9071-0 (not retrieved) |
| Hynes, Angove, Saunders, Haverd, Azzi 2005 AE 39, 7251 (propene, 1-butene) | CSIRO indoor chamber; "mercury blacklight blue lamps emitting roughly 350–390 nm" [S] | MCM v3.1 | **D**: measured NO2 photolysis rate [S] | constant in time; "the photolysis rates in MCM had to be adjusted accordingly … based on measured NO2 photolysis rates" (Ginnebaugh et al. 2010 SI) [V of secondary]; method (A vs B) not determinable from secondary text; sensitivity to "key photolysis rates" studied (**E-ish**) | own wall-reaction table (Teflon walls; provided to Ginnebaugh et al.), Carter-style | 10.1016/j.atmosenv.2005.09.005 (full text not retrieved) |
| Wagner et al. 2003 ACP 3, 89 (toluene) | EUPHORE (outdoor, sunlight) | MCM v3 | **D**: spectroradiometer (Bentham DM300) | **D** for O3, NO2, HCHO, HONO; **B** for all others: "boundary layer photolysis rates were calculated for the chamber location and experiment date (Saunders 2002; Jenkin 2002)" × one factor = average measured/calculated deviation for O3, NO2, HCHO, HONO | EUPHORE: k(NO2→HONO wall) = (1.15 ± 0.45)e-5 s⁻¹; dilution | 10.5194/acp-3-89-2003 [V] |
| Bloss et al. 2005a ACP 5, 623 (MCMv3/3.1 aromatics evaluation) | EUPHORE | MCM v3, v3.1 | **D**: filter radiometer J(NO2) (two instruments) | **B** with two factors: f1 = J(NO2)_meas / J(NO2)_clear-sky-calc (time-varying, applied to all J); f2 = species-dependent wall/cloud transmission factor from the 22 Oct 1997 dataset (O3, HCHO, HONO measured); average f2 for all other species | EUPHORE characterisation + literature: Table 2 (NO2 → 0.5 HONO + 0.5 wHNO3 at 1.15e-5 s⁻¹; light-induced NO2 off-gassing; N2O5 hydrolysis; O3 and HNO3 wall loss; initial HONO 0.5 ppb); limits from UNC/SAPRC values of **Carter & Lurmann (1991)** and Tuazon et al. (1983) | 10.5194/acp-5-623-2005 [V] |
| Bloss et al. 2005b ACP 5, 641 (MCMv3.1 aromatics development) | EUPHORE | MCM v3.1 | **D** (as 2005a) | **B** (as 2005a: "calculated rates were used with scaling factors …"); α-dicarbonyl J set to measured 0.0057·J(NO2) etc. (**D** for those) | same as 2005a | 10.5194/acp-5-641-2005 [V] |
| Zádor, Turányi, Wirtz, Pilling 2006 J. Atmos. Chem. 55, 147 (EUPHORE radical sources) | EUPHORE | MCM | **D** [?] | **B** [?] (EUPHORE standard practice) | own: empirical HONO (and HCHO) wall-source expressions as f(j(NO2), T, RH) [S]; Eq. (1) for W(HONO) quoted by Metzger 2008 [V] | 10.1007/s10874-006-9033-y (not retrieved) |
| Metzger et al. 2008 ACP 8, 6453 (1,3,5-TMB) | PSI 27 m³ chamber, 4 × 4 kW Xe arc; also SAPRC CTC/DTC datasets | MCM v3.1 (+ updates) | **D**: NO2 actinometry (Paulsen et al. 2005), ±12 %; trend tracked by J(NO2) filter radiometer | **A**: "determined from the measured light spectrum and using the absorption cross sections and quantum yield data summarized by Jenkin et al. (1997) for the MCM" (Bentham DMc 150 FC spectrograph) | own characterisation (12 clean-air, 2 NOx-air, 4 CO-air, 3 CO–NOx-air runs); N2O5 wall rates "adopted from (Carter and Lurmann, 1991) and (Bloss et al., 2005a)"; for SAPRC chambers: Carter et al. (1995a) updated by Pinho et al. (2005) | 10.5194/acp-8-6453-2008 [V] |
| Rickard et al. 2010 AE 44, 5423 (TMB, PSI, MCM) / Wyche et al. 2009 ACP 9, 635 | PSI chamber, same campaign | MCM v3.1 | as Metzger 2008 [?] | as Metzger 2008 [?] | as Metzger 2008 [?] | Wyche 10.5194/acp-9-635-2009 [V, no J method stated]; Rickard not retrieved |
| Paulsen et al. 2005 ES&T 39, 2668 (PSI chamber characterisation) | PSI 27 m³, 4 Xe arc lamps (16 kW) | — | **D**: NO2 actinometry (cited by Metzger) | — | — | 10.1021/es0489137 (not retrieved) |
| Platt et al. 2013 ACP 13, 9141 (PSI mobile chamber) | 40 × Philips Cleo Performance 100 W, peak 368 nm | MCM-independent J table (no MCM run reported) | **D**: NO/NO2/O3 photostationary state → (8.0 ± 0.7)e-3 s⁻¹ | **A**: Ji = JNO2 · ∫σiφi I dλ / ∫σNO2 φNO2 I dλ with IUPAC σ, φ and the measured (relative) lamp spectrum; Table 5 | HONO injection as OH source; no auxiliary mechanism reported | 10.5194/acp-13-9141-2013 [V] |
| Glowacki et al. 2007 ACP 7, 5371 (HIRAC design) | HIRAC, Philips TL-K 40W/05 blacklamps in 8 quartz tubes | — (MCM used later) | **D**: NO2 actinometry (NO2 in air) → (1.93 ± 0.10)e-3 s⁻¹ integrated; spatial map by ray tracing (OptisWorks) with blacklamp F(λ) from Carter et al. (2005, 1995) | **A** in principle: "JNO2 = F(x)∫σ φ F(λ) dλ"; blacklamp spectrum taken from Carter | not discussed | 10.5194/acp-7-5371-2007 [V] |
| Malkin et al. 2010 ACP 10, 1441 (isoprene ozonolysis, HIRAC) | HIRAC, dark | MCM v3.1 | n/a (dark) | n/a | none mentioned | 10.5194/acp-10-1441-2010 [V] |
| Winiberg et al. 2015 AMT 8, 523 (HIRAC FAGE calibration) | HIRAC blacklamps | — | instrument paper; no chamber J specification | — | — | 10.5194/amt-8-523-2015 [V] |
| Winiberg et al. 2016 ACP 16, 4023 (CH3C(O)O2 + HO2, HIRAC) | HIRAC, Philips TL-D 36W/BLB (350–400 nm), 2/4/8 lamps | Kintecus model (MCM-based rates) | **D**: j(NO2) = (2.4 ± 0.8)e-3 s⁻¹ (8 lamps) | **A** for j(Cl2) estimate from IUPAC σ (3.8e-4 s⁻¹) then **E**: "j(Cl2) has been determined for each set of experimental data using the measured Cl-atom-induced decay"; measured lamp-intensity time profile used as a model constraint | none (short experiments) | 10.5194/acp-16-4023-2016 [V] |
| Wang et al. 2011 AMT 4, 2465 (CESAM design) | CESAM 4.2 m³ stainless steel, 3 × 4 kW Xe arc, Pyrex filter | MCM v3 (propene) + auxiliary; FACSIMILE | **D** by fit (**E**): NOx-air photostationary-state runs fitted with a 36-reaction model → (4.2 ± 0.2)e-3 s⁻¹; cross-checked with LI-COR LI-1800 spectroradiometer (4.5e-3) | **A** for NO2 and O3→O(1D) (irradiance → actinic flux with pseudo-zenith geometry, σ from Voigt/Molina, φ from DeMore); **B** for everything else in the MCM run: "the other photolysis frequencies used in the MCM were scaled from their atmospheric values using JNO2 as a reference value" | own: NO2 wall uptake/re-release, `NO2(wall) → NO + HONO` dark and photo-enhanced (k = 1.4e-2 s⁻¹ × 0.51/0.10), O3 wall loss 3e-4 s⁻¹; "added without any tuning"; alternatives from Carter 2000/2006, Bloss 2005, Carter 2005, Hynes 2005, Metzger 2008 were tested | 10.5194/amt-4-2465-2011 [V] |
| Picquet-Varrault et al. 2020 ACP 20, 487 (carbonyl nitrates, CESAM) | CESAM Xe arc | — (kinetic model of Wang 2011) | **D**: dedicated NO2-in-N2 runs fitted → 2.2e-3 s⁻¹ (lamps aged) | **D** for the nitrates themselves; reported as j(nitrate)/j(NO2) "to allow for the calculation of jnitrate under various sunlight conditions" | own (Wang 2011) | 10.5194/acp-20-487-2020 [V] |
| Fuchs et al. 2014 ACP 14, 7895 (MACR; method identical to Fuchs 2013 Nat. Geosci. isoprene) | SAPHIR (outdoor, sunlight) | MCM v3.2 | **D**: spectroradiometer (actinic flux) | **D** for NO2, HCHO, O3, MACR, HONO; **B** for the rest: "first calculated by a MCM (version 3.1) function for clear-sky conditions and then scaled by the ratio of the measured and calculated NO2 photolysis frequencies" | SAPHIR own: HONO, HCHO, CH3CHO sources as f(radiation, RH, T) (Rohrer et al. 2005; Karl et al. 2006), scaled to the zero-air phase; background OH reactivity; photolytic RO2 source | 10.5194/acp-14-7895-2014 [V]; Fuchs 2013: 10.1038/ngeo1964 (methods not retrieved) |
| Jenkin et al. 2015 ACP 15, 11433 (MCM v3.3.1 isoprene) | no chamber runs in the paper; boundary-layer box model | MCM v3.3.1 | — | MCM parameterisation, "clear sky conditions … latitude of 0° at equinox" (ambient illustration, not a chamber) | — | 10.5194/acp-15-11433-2015 [V] |
| Schwantes et al. 2017 ACP 17, 3453 (toluene/cresol, Caltech) | Caltech 24 m³ Teflon, "UV lights centered ~350 nm" | MCM v3.3.1-based kinetic model | "calculated" j(NO2) = 4.4e-3 s⁻¹ (basis not stated in main text; Caltech practice = actinometry-anchored lamp spectrum [?]) | **A** [?]: "calculated photolysis rates for NO2, H2O2, and NO3 are 4.4e-3, 3.2e-6 and 2.3e-3 s⁻¹"; nitrophenol J from literature "corrected for the difference in the photolysis of NO2 between the Caltech chamber and the atmosphere" (**C**-type scaling for one species) | none beyond first-order wall losses | 10.5194/acp-17-3453-2017 [V] |
| Boyd et al. 2015 ACP 15, 7497 (β-pinene + NO3, GTEC) | GTEC 12 m³ Teflon, blacklights (300–900 nm emission range quoted), j(NO2) = 0.28 min⁻¹ all lamps | "simple kinetic model based on MCMv3.2" for dark NO3 chemistry | **D** (actinometry, value quoted) | n/a (dark experiments) | wall loss only | 10.5194/acp-15-7497-2015 [V] |
| Nah et al. 2016 ACP 16, 9361 (GTEC vapour wall loss) | GTEC | no MCM | — | — | — | 10.5194/acp-16-9361-2016 [V, not relevant to J] |
| Peng et al. 2015 AMT 8, 4863 (OFR185/OFR254 model; KinSim lineage) | Aerodyne PAM OFR, low-pressure Hg lamps (185 + 254 nm) | plug-flow kinetic model (JPL/IUPAC σ) | n/a (no NO2) | **A** with *line* fluxes: F185 = 1e11–1e14, F254 = 4.2e13–8.5e15 photons cm⁻² s⁻¹ "based on the characterization of Li et al. (2015)"; J = σ·F per line; flux uncertainty flagged as dominant near lamp threshold | n/a | 10.5194/amt-8-4863-2015 [V] |
| Sommariva et al. 2020 GMD 13, 169 (AtChem2) | tool | MCM | user: constrained (`photolysisConstrained.config`) or constant (`photolysisConstant.config`) | **B** by default: MCM l,m,n at the real SZA × τ (chamber-wall transmission factor) × JFAC = j_meas/j_calc (usually NO2); **or** constant J per reaction ("The typical scenario for constant photolysis rates is a lamp or solar simulator") | user | 10.5194/gmd-13-169-2020 [V]; AtChem2 manual [V] |
| Wolfe et al. 2016 GMD 9, 3309 (F0AM) | tool | MCM or others | user: direct J input, or global scaling factor "taken as the ratio of an observed photolysis frequency to its model-calculated value" | three options: **"bottom-up"** = user actinic-flux spectrum × IUPAC/JPL σφ ("most useful when simulating photochemical chamber experiments with non-solar light sources") = **A**; **"MCM"** l,m,n (**B**); **"hybrid"** TUV look-up tables; chamber example §3.1 uses bottom-up with a 350 nm UV-bulb spectrum from Crounse et al. (2011) | user (first-order wall loss) | 10.5194/gmd-9-3309-2016 [V] |
| Leskinen et al. 2015 AMT 8, 2267 (ILMARI) with referee comment by W. Carter | ILMARI, 365 nm "BLB" lamps then 350 nm lamps | — | NO2 actinometry | — | — | 10.5194/amt-8-2267-2015; referee comment amtd-7-C5010-2015 [V] |
| Carter, Atkinson, Winer, Pitts 1982 Int. J. Chem. Kinet. 14, 1071 "Experimental investigation of chamber-dependent radical sources" (the user's "Carter 1984") | SAPRC evacuable chamber | — | — | — | origin of the chamber radical source concept used in all Carter auxiliary mechanisms | 10.1002/kin.550141003 (crossref-verified citation; text not retrieved) |
| Carter & Lurmann 1991 AE 25A, 2771 | UNC outdoor + SAPRC chambers | SAPRC-90 | — | — | origin of the UNC/SAPRC wall parameters quoted as limits by Bloss 2005a and Metzger 2008 | not retrieved (cited via Bloss 2005a, Metzger 2008, Carter 2000) |

Not found in this corpus: a paper stating "the MCM photolysis parameterisation was used with a fixed solar
zenith angle to mimic the lamps". The closest things are (i) the EUPHORE/SAPHIR practice (real SZA plus
JFAC-type scaling), (ii) CESAM 2011 (MCM atmospheric J scaled by j(NO2), justified by the realism of the
Xe spectrum), and (iii) AtChem2's default, which does exactly that if a user gives a lamp-lit chamber a
latitude/date instead of a `photolysisConstant.config` file.

---

## 2. Classification of methods actually used (counts)

Counting *studies* (not tools), with the dominant method for the "all other J" column; studies marked
[?] counted where the inference is strong (Pinho 2006/2007 follow Pinho 2005; Rickard 2010 follows
Metzger 2008). Mixed cases are listed in both classes where appropriate.

| Class | Studies | Count |
|---|---|---|
| **A** spectrum × σφ anchored on NO2 actinometry | Carter 1995 (ARB), Carter 2000 (SAPRC-99), Carter 2005 (UCR EPA), Carter 2010 (SAPRC-07), Pinho 2005 [S], Pinho 2006/2007 [?], Metzger 2008, Rickard 2010 [?], Platt 2013, Glowacki 2007 (j(NO2) map), Winiberg 2016 (j(Cl2) prior), Wang 2011 (NO2, O3 only), Schwantes 2017 [?], Peng 2015 (Hg lines) | **14** (10 firmly verified) |
| **B** solar parameterisation at the real SZA × j(NO2)_meas/j(NO2)_calc (± species factor) | Wagner 2003, Bloss 2005a, Bloss 2005b, Zádor 2006 [?], Fuchs 2013/2014, Wang 2011 (all J except NO2, O3) | **6** |
| **C** ratios from a reference chamber / canonical spectrum | Carter 2000 blacklight runs ("spectrum is assumed to be constant and the spectrum recommended by Carter et al (1995d) is used"); Glowacki 2007 (blacklamp F(λ) from Carter); Schwantes 2017 nitrophenol (literature J × j(NO2) ratio) | **3**, all as sub-steps of A |
| **D** measured directly for the species in question | j(NO2) in essentially every study (actinometry or radiometer); additionally O3, HCHO, HONO at EUPHORE (Wagner, Bloss); NO2, HCHO, O3, MACR, HONO at SAPHIR (Fuchs); α-dicarbonyls at EUPHORE (Bloss 2005b); carbonyl nitrates at CESAM (Picquet-Varrault 2020) | j(NO2): **all**; other species: **4** chambers |
| **E** fitted | Wang 2011 (j(NO2) by model fit to NOx-air runs — strictly an actinometry evaluation); Winiberg 2016 (j(Cl2) from reactant decay); Hynes 2005 (sensitivity on "key photolysis rates") [S]; Carter 2000/2010 fit *quantum yields* (MEK, DCB1-3, glyoxal radical yield) to chamber data rather than J | **2 (+ fitted φ in SAPRC)** |

Reading of the counts:

1. **Indoor-lamp chambers modelled with SAPRC or MCM used class A** whenever the group had a
   spectroradiometer (UCR, PSI, Caltech/F0AM, HIRAC, CESAM for NO2/O3). This is also the method
   embedded in the SAPRC modelling software that Pinho/Jenkin used for the MCM evaluations.
2. **Class B is a sunlight-chamber method** (EUPHORE, SAPHIR). The one lamp chamber that used it (CESAM
   2011) explicitly justified it by the solar-likeness of the filtered Xe spectrum, and later CESAM work
   (2020) reports j(X)/j(NO2) ratios measured in the chamber instead.
3. Nobody in this corpus ran the MCM l,m,n at a fictitious SZA for blacklights. The two public MCM box
   models differ: AtChem2 offers "constant J per reaction" for lamps (and otherwise applies SZA × JFAC),
   F0AM offers "bottom-up" (class A) and recommends it for non-solar light.
4. The chamber auxiliary mechanism is **Carter's in lineage everywhere**: UCR/SAPRC chambers (Carter
   1982 → 1991 → 1995 → 2000 → 2005, updated by Pinho 2005 for the MCM runs), EUPHORE (own
   characterisation, with Carter & Lurmann 1991 UNC/SAPRC values as bounds), PSI (own characterisation,
   N2O5 rates from Carter & Lurmann 1991 / Bloss 2005a), CESAM (own, after testing Carter 2000/2006 and
   others), SAPHIR (own, Rohrer 2005 — which Carter 2005/2010 in turn cite as the direct confirmation of
   his HONO-off-gassing representation).

---

## 3. What the mechanism-development papers did (verbatim)

These define "what SAPRC and the MCM consider correct" for chamber photolysis.

### 3.1 SAPRC — Carter 2000, *Documentation of the SAPRC-99 Chemical Mechanism*, Section V.B.1 "Light Characterization" [V]

> "Light characterization requires specification of both the intensity and the spectrum of the light
> source used in the experiments. As discussed by Carter et al (1995c,d) for indoor chamber runs, this is
> determined by the NO2 photolysis rate (usually derived from results of NO2 actinometry experiments),
> and the relative spectral distribution of the light source. For blacklight chambers the spectrum is
> assumed to be constant and the spectrum recommended by Carter et al (1995d) is used. For Xenon arc
> chambers, the spectrum tends to vary with time, and the spectrum used for modeling is based on
> measurements made during or around the time of the experiment as discussed by Carter et al (1995d).
> For the outdoor chamber, the spectrum was calculated as a function of solar zenith angle as described
> by Carter et al. (1995c)."

> "For the blacklight chambers and the EC, the NO2 photolysis rates were determined by carrying out
> periodic NO2 actinometry experiments, with the values assigned to individual runs being based on the
> trends or averages of the measurements carried out around the time of the experiment (Carter et al,
> 1995d, 1997a). For the XTC and CTC, the relative trends in light intensity with time were determined
> primarily using absolute spectral intensity measurements, whose data were placed on an absolute basis
> using occasional in-chamber steady-state or Cl2 actinometry measurements (Carter et al, 1995d, 1997a)."

Section V.B.2 "Representation of Chamber Wall Effects" [V]:

> "The chamber wall effects that were represented in the simulations of these experiments were the
> chamber radical source (Carter et al, 1982), NOx offgasing, heterogeneous formation of HONO from NO2,
> N2O5 hydrolysis, O3 dark decay, and background effects causing excess NO to NO2 conversions. These by
> the following psuedo-reactions: Walls + hν → HONO (RN); Walls + NO2 + hν → 0.5 {HONO + Wall-NOx} (RS);
> Walls + NO2 → yHONO HONO + (1-yHONO) Wall-NOx (NO2W); N2O5 + walls → 2 Wall-NOx (N25I);
> N2O5 + H2O + walls → 2 Wall-NOx (N25S); Walls + O3 → (loss of O3) (O3W)."

> "The formation of HONO from the walls (Reaction RN) us used to represent both the chamber radical
> source and NOx offgasing, because the HONO so formed would photolyze rapidly to form both OH radicals
> and NO … Previously, the radical source and NOx offgasing were represented as separate processes …
> However, in most cases the rates of these two reactions tended to be the same to within the uncertainty
> of the determination, suggesting that they may be due to the same process."

### 3.2 SAPRC — Carter et al. 2004/2005 UCR EPA chamber report, "Light Characterization" [V]

> "Photolysis rates used when modeling chamber experiments are calculated using the measured NO2
> photolysis rates, the relative measured spectral distributions for the light sources, and the
> absorption cross sections and quantum yields for NO2 and the other photolysis reactions in the chemical
> mechanism being evaluated. Therefore the measured NO2 photolysis rates serve as the measure of the
> absolute light intensity, and the relative spectral distributions of the light sources serve as the
> means to calculate the other photolysis rates relative to that for NO2. The precisions of the
> photolysis rates so derived are determined primarily by the precision of the NO2 actinometry
> measurement."

> "… they were all assumed to have the same NO2 photolysis rate for modeling purposes, which was
> 0.26 min-1, the average of the in-chamber actinometry results. The other photolysis rates were
> calculated using this, the assigned spectral distribution (shown on Figure 3) and the absorption
> cross-sections and quantum yields given with the mechanism. The relative spectrum of the arc light
> source was measured using a LI-COR LI-1800 spectroradiometer …"

> "Blacklights tend to decrease in intensity in time as they are used (Carter et al, 1995b), and this
> needs to be taken into account when assigning photolysis rates when modeling blacklight experiments."
> (In-reactor quartz-tube actinometry is "the standard to derive absolute photolysis rates within the
> chamber"; the 'standard position' tube reads high by 1/0.698; a QSL spherical irradiance sensor,
> calibrated against in-chamber actinometry, tracks the decline against "blacklight run count".)

Chamber radical source [V]:

> "It has been known for some time that environmental chamber experiments could not be modeled
> consistently unless some sources of radicals attributed to chamber effects is assumed (e.g., Carter et
> al, 1982; Carter and Lurmann, 1991; Carter, 2000). The most sensitive experiments to this effect are
> NOx-air irradiations of compounds, such as CO or alkanes …"
> "For this chamber, assuming HONO offgasing at a similar magnitude as the apparent NOx offgasing rate …
> is usually sufficient to account for most of the chamber-dependent radical source, though results of
> some of the experiments are somewhat better simulated if a small amount (100 ppt or less) of HONO is
> also assumed to be initially present."
> "Direct evidence for this comes from the data of Rohrer et al (2004), who used sensitive long path
> absorption photometer (LOPAP) instrument to detect ppt levels of HONO emitted from the walls during
> irradiations in the large outdoor SAPHIR chamber …"

### 3.3 SAPRC-07 — Carter 2010 report [V]

> "The ~7% increase in the NO2 photolysis rate under atmospheric conditions has no effect on NO2
> photolysis rates used in chamber simulations because all photolysis rates are normalized to measured
> NO2 photolysis rates. However, rates of other photolysis reactions in chamber simulations will
> decrease accordingly."

> "This is represented in the chamber model by the parameter RN, which is the rate of HONO offgasing
> relative to the light intensity as measured by the NO2 photolysis rate. Since HONO has not been
> measured directly in any of the chambers used for mechanism evaluation … the HONO offgasing rate
> parameter has to be determined by adjusting the parameter so the model calculations can simulate
> results of the appropriate characterization experiments. The most sensitive experiments are the CO -
> NOx and alkane (primarily n-butane) - NOx experiments …"
> Table 22: "RN-I (ppb) Ratio of the rate of wall + hν -> HONO to the NO2 photolysis rate. EPA 0.0075;
> ITC 0.0165; ETC/DTC 0.0115; XTC/CTC 0.0055" and an optional "Walls + hν → WALLVOC; rate = NO2
> photolysis rate x 0.135 ppb; WALLVOC + OH → HO2 + 0.2 HCHO" (Simonaitis et al. 1997).

Note the structural choice that matters for a chamber code: **every chamber-wall source is expressed as
a ratio to j(NO2)**, so that it scales automatically with lamp intensity and with lamp ageing.

### 3.4 MCM evaluation against SAPRC chambers — Pinho, Pio, Jenkin 2005 [S]

Search-engine snippets of the paper (full text not retrievable; consistent with Carter's method above):

> "rate constants for the photolysis reactions are calculated using the NO2 photolysis rates
> (experimentally determined, 0.21 min-1), the relative spectral distributions for the light sources,
> the absorption cross sections and quantum yields" [S]

> "the rates applied to photolysis processes were based on spectral distributions and absolute light
> intensities which are fully documented for all the SAPRC chambers" [S]

Corroboration from full texts: Metzger et al. 2008 [V]: "For the simulations in this study we used the
auxiliary mechanism of each of the chambers as described in Carter et al. (1995a) and updated by Pinho et
al. (2005)." and "The rates of a number of photolysis reactions were updated according to Pinho et al.
(2005)." Pinho et al., EGU 2006 abstract [V]: "the MCM v3 schemes were evaluated using environmental
chamber data by adaptation of MCM v3 with inclusion of a chamber-dependent auxiliary mechanism … in
conjunction with an initial evaluation of the chamber-dependent auxiliary mechanisms for the series of
relevant chambers." The MCM photolysis page states that the MCM's MACR/MVK quantum yields come from
"the optimised QY used in the evaluation of the MCM and SAPRC-99 using chamber data" [S] — i.e. the
MCM's own φ values for these carbonyls are blacklight-chamber fits made with class-A photolysis rates.

Consequence: **the MCM's isoprene/alkene/pinene/butane/carbonyl chemistry was validated against
blacklight and Xe-arc chambers with class-A photolysis rates computed by Carter's procedure, not with
the MCM solar parameterisation.**

### 3.5 MCM evaluation against EUPHORE — Wagner et al. 2003 [V]

> "For the measurement of the photolysis rates of O3 NO2, HCHO and HONO a spectroradiometer (Bentham
> DM300) was employed."
> "For O3, NO2, HCHO and HONO measured photolysis rates were taken for the model simulation. For all other
> photolysis processes, the corresponding boundary layer photolysis rates were calculated for the chamber
> location and experiment date, as described elsewhere (Saunders et al., 2002; Jenkin et al., 2002). These
> calculations do not consider the actual cloud coverage, the backscattering of the chamber floor and the
> transmission of the chamber FEP film. Therefore we have adjusted the calculated photolysis rates with a
> factor that is based on the average deviation between measured and calculated photolysis rates for O3,
> NO2, HCHO and HONO."
> "… for the surface process (NO2 → HONO) an average rate coefficient kHONO,wall of (1.15±0.45)×10-5 s-1
> has been established for the EUPHORE chamber …"

### 3.6 MCM evaluation against EUPHORE — Bloss et al. 2005a [V]

> "The measured photolysis rate for NO2 is available and for all other photolysis processes calculated
> rates were used with scaling factors applied to account for the transmission of the chamber walls,
> backscattering from the chamber floor and cloud cover. The variations from day to day and during the
> experiments resulting from cloud cover are accounted for by considering the difference between
> measured and calculated clear sky J(NO2) at any given time during the experiments. This variable
> scaling factor, f1, is applied to all photolysis rates in addition to a factor, f2, designed to account
> for wall and cloud transmission effects, which will be different for different species dependent on
> their absorption spectrum. In a previous chamber dataset (22 October 1997) the photolysis rates for O3,
> HCHO, HONO were available in addition to J(NO2), and for these species the applied scaling factors, f2,
> are based on the deviation between their measured and calculated photolysis rates, normalised by the
> deviation seen for NO2. For all other photolysis rates the average value of these factors is used."

> "An auxiliary mechanism that considers these chamber related effects has been constructed and the
> relevant reactions and associated rate coefficients are listed in Table 2. Where data from EUPHORE
> characterisation experiments are not available, literature data adapted to EUPHORE conditions are
> taken, and upper and lower limits to these rates have been assigned …"
> Table 2 (nominal / lower / upper): NO2 → 0.5 HONO + 0.5 wHNO3, 1.15e-5 / 0.7e-5 / 1.6e-5 s⁻¹ ("EUPHORE
> characterisation"); wHNO3 + hν → OH + NO2, JHNO3 / 0 / 2·JHNO3; light-induced NO2 off-gassing
> 0.1·JNO2 ppb s⁻¹ / 0.025·JNO2 / 0.2·JNO2 ("values for UNC chamber (Carter and Lurmann 1991)"); N2O5 →
> 2 wHNO3, 3e-5 / 1.3e-5 / 4.7e-5 s⁻¹ (Tuazon et al. 1983, scaled for S/V); N2O5 + H2O, 6.5e-20 /
> 1.1e-23 / 1.3e-19 cm³ s⁻¹; OH → HO2 (background reactivity) 8.4 / 4.2 / 12.5 s⁻¹ ("upper limit as UNC
> chamber, lower limit as SAPRC chamber (Carter and Lurmann 1991)"); O3 → wO3 3e-6 s⁻¹; HNO3 → wHNO3
> 8.2e-5 s⁻¹; initial HONO 0.5 / 0 / 1.0 ppb.

Bloss et al. 2005b [V] repeats the same sentence and adds measured J for the α-dicarbonyls: "the
photolysis rates for these compounds were changed to 0.0057×JNO2 in accordance with the measurements".

### 3.7 MCM evaluation against the PSI Xe-arc chamber — Metzger et al. 2008 [V]

> "3.1 Photolysis processes. The rates of photolysis processes depend on the intensity and spectral
> distribution of the light source. Spectral distributions are based on spectrometer measurements
> (Bentham Spectrograph DMc 150 FC), whereas absolute light intensities were determined by NO2
> actinometry experiments (Paulsen et al., 2005). The time trend of the light intensity was tracked using
> J(NO2) filter radiometer (FR) data and regular NO2 actinometry experiments. The light intensity
> measured by NO2 actinometry and J(NO2)-FR stayed constant (±4%) over the time period considered. An
> uncertainty of ±12% was attributed to the J(NO2) value measured with chemical actinometry."
> "The photolysis rates were determined from the measured light spectrum and using the absorption cross
> sections and quantum yield data summarized by Jenkin et al. (1997) for the MCM. The rates of a number of
> photolysis reactions were updated according to Pinho et al. (2005)."

Auxiliary mechanism (Table 2, PSI characterisation unless noted): hν + wall → HONO 9.1e6 molecule cm⁻³
s⁻¹ (5e6–13e6); [HONO]0 = 0.0074·RH + 0.006·[NO2]0; NO2 → 0.5 HONO + 0.5 wHNO3 0.53e-6 s⁻¹; hν + wall →
HCHO 5e6 cm⁻³ s⁻¹; N2O5 → 2 wHNO3 1e-5 s⁻¹ "Adopted from (Carter and Lurmann, 1991) and (Bloss et al.,
2005a)"; OH + X → HO2, [X]0 = 300 ppb "Adopted from (Bloss et al., 2005a)"; O3 wall 4e-6 s⁻¹; HNO3 wall
1e-4 s⁻¹; wHNO3 + hν → NO2 + OH at JHNO3; dilution; J(NO2) ±12 %. HONO wall production derived with
Eq. (1) "(Zador et al., 2006)".

### 3.8 MCM against SAPHIR — Fuchs et al. 2014 (same procedure as Fuchs et al. 2013) [V]

> "The model was constrained by measurements of temperature, pressure, calculated dilution rates,
> measured water vapour mixing ratios, and photolysis frequencies for NO2, HCHO, O3, MACR, and HONO.
> Photolysis frequencies and relative humidity constrained the calculated chamber sources for HONO, HCHO,
> and CH3CHO. Photolysis frequencies that were not measured were first calculated by a MCM (version 3.1)
> function for clear-sky conditions and then scaled by the ratio of the measured and calculated NO2
> photolysis frequencies to account for reductions of the solar radiation by clouds and the transmission
> of the Teflon film."
> "The dependence of source strengths for HONO, HCHO, and CH3CHO on radiation, relative humidity, and
> temperature were taken from parameterizations described earlier (Rohrer et al., 2005; Karl et al.,
> 2006). The source strengths were scaled in order to match formation of NOX, HCHO, and CH3CHO during the
> zero-air part of the experiment."

### 3.9 MCM v3.3.1 — Jenkin et al. 2015 [V]

No chamber simulations are presented; the illustrative box model uses "Photolysis coefficients …
calculated assuming clear sky conditions, using a parameterization previously applied with the MCM
(Saunders et al., 2003)" for "a latitude of 0° at equinox". Chamber evidence enters only through cited
studies (Fuchs 2013/2014 for SAPHIR; Pinho 2005 for the SAPRC chambers).

### 3.10 The MCM parameterisation itself (basis)

Wolfe et al. 2016 [V]: "solar spectra underlying this parameterization were calculated from a two-stream
radiative transfer model for clear sky conditions on 1 July at a latitude of 45° N and an altitude of
0.5 km." Sommariva et al. 2020 [V]: "the empirical parameters are obtained by fitting Eq. (3) to the
output of a two-stream isotropic scattering model, which incorporates the appropriate photolysis
cross-sections and quantum yields." Neither the MCM protocol papers nor the box-model papers claim the
l,m,n set represents any artificial light source.

---

## 4. Per-study notes beyond the table (verbatim where useful)

**CESAM, Wang et al. 2011 [V]** — the one lamp chamber in the corpus that used class B for the MCM run,
and the justification given: "The auxiliary mechanism was added without any tuning and the NO2
photolysis frequency was entered as determined in the NOx-air system simulation section. Considering
the realism of the spectral distribution of our light sources, the other photolysis frequencies used in
the MCM were scaled from their atmospheric values using JNO2 as a reference value. Only [HONO]0, the
initial concentration of nitrous acid – was adjusted as HONO measurements were not available during the
time of these experiments." Class A was used for the characterisation itself: "Photolysis frequencies
(J) were calculated by using the cross sections from Voigt et al. (2001) for NO2 and from Molina and
Molina (1986) for ozone. The quantum yields for ozone and NO2 are from DeMore et al. (1994). Irradiance
spectra taken with and without Pyrex filters with the LICOR LI-1800 instrument were then converted into
actinic flux … taking into account pseudo-zenithal angles from geometrical calculations in the chamber."
j(NO2) itself: "The fitting of modeled values to the measured data led to NO2 photolysis frequency equal
to (4.2 ± 0.2) × 10-3 s-1 … in good agreement with the single point spectroradiometer determination
presented earlier (4.5 × 10-3 s-1)." Table 3 of the paper reports JNO2/JO1D for several Pyrex
thicknesses against TUV ground level — i.e. the group checks the *ratio* j(O1D)/j(NO2) before trusting
class B. Auxiliary mechanism: "Several plausible auxiliary mechanisms were tested including those
proposed in the literature (Doyle, 1970; Jeffries et al., 1976; Akimoto et al., 1979b; Carter, 2000,
2006; Bloss et al., 2005; Carter et al., 2005; Hynes et al., 2005; Metzger et al., 2008)."

**CESAM, Picquet-Varrault et al. 2020 [V]** — direct measurement with j(NO2) as the yardstick: "as the
intensity of the irradiation in the CESAM chamber is lower than that in the ambient environment,
jnitrate/jNO2 ratios were provided in order to allow for the calculation of jnitrate under various
sunlight conditions. Hence, the intensity of the actinic flux was determined by measuring the photolysis
frequency of NO2 (JNO2) during dedicated experiments." (400 ppbv NO2 in N2, 20 min dark / 20 min light /
20 min dark, fitted with the Wang 2011 model → 2.2e-3 s⁻¹.)

**PSI mobile chamber, Platt et al. 2013 [V]** — the cleanest statement of class A in the MCM literature:
"JNO2(λ,T) = ∫ σNO2(λ,T)·φNO2(λ,T)·I(λ,T) dλ … Thus, while I(λ,T) in Fig. 6a is given by the light
photometer in arbitrary units only, absolute values for the photolysis rates of a compound i (Ji) can be
inferred from JNO2: Ji = JNO2 · ∫σi φi I dλ / ∫σNO2 φNO2 I dλ (4). Inferred photolysis rates at 25 °C from
Eq. (4) are given for several compounds in Table 5, where IUPAC recommended values of σ(λ,T) and φ(λ,T)
(Atkinson et al., 2004, 2006) were used." Also: "The photolysis rates of some compounds, notably ketones
and aldehydes, may be enhanced compared to ambient conditions when using UV lights (Carter, 1995)."

**HIRAC, Glowacki et al. 2007 [V]** — "The NO2 photolysis frequency is the product of the temperature and
wavelength dependent absorption cross section, σNO2(λ,T), multiplied by its quantum yield, φNO2(λ,T),
and the intensity of radiation at a particular wavelength and position with respect to a lamp, F(λ,x).
σNO2(λ,T) and φNO2(λ,T) are well established (DeMore et al., 1997; Voigt et al., 2002), as is F(λ) for
blacklamps (Carter et al., 2005, 1995)." Integrated j(NO2) = (1.93 ± 0.10)e-3 s⁻¹ by NO2 actinometry.

**HIRAC, Winiberg et al. 2016 [V]** — "Supporting j(NO2) measurements were conducted for all eight of the
black lamps switched on during the later time period, which gave j(NO2) = (2.4 ± 0.8) × 10-3 s-1. Using
IUPAC recommended absorption cross sections for both NO2 (Atkinson et al., 2004) and Cl2 (Atkinson et
al., 2007), j(Cl2) was estimated at (3.8 ± 1.4) × 10-4 s-1, correlating well with the observed reactant
decays. The black lamp intensity profiles as a function of time were entered into Kintecus as a
constraint for the photolysis rate … Both the predicted [OH] and [HO2] were observed to better correlate
with the measured radical concentrations when using this constraint, compared to starting the model
with a constant photolysis rate." (Lamp warm-up matters at the <40 min timescale.)

**Caltech, Schwantes et al. 2017 [V]** — "Ultraviolet (UV) lights centered ~350 nm are used in the Caltech
chamber. The calculated photolysis rates for NO2, H2O2, and NO3 are 4.4 × 10-3, 3.2 × 10-6, and
2.3 × 10-3 s-1, respectively. Because the UV lights emit only a low amount of visible light, NO3
photolysis is low and uncertain in the chamber." Nitrophenol: "This photolysis rate constant was
corrected for the difference in the photolysis of NO2 between the Caltech chamber and the atmosphere
when/where Bardini made the measurement." The F0AM paper's chamber example (same community) uses
"the 'bottom-up' method with a light spectrum corresponding to UV bulbs with output centered at 350 nm
(Crounse et al., 2011)" [V].

**CSIRO, Hynes et al. 2005 (via Ginnebaugh et al. 2010 SI) [V of the SI]** — "The chamber used a UV lamp
instead of sunlight to initiate photolysis, so the light intensity in the model was kept constant
throughout the experiment. The UV lamp did not represent the full range of wavelengths from sunlight, so
the photolysis rates in MCM had to be adjusted accordingly (Hynes et al. 2005). These rates were based on
measured nitrogen dioxide (NO2) photolysis rates and were provided by the authors." Wall reactions:
"Measured and Provided by Hynes et al. (2005)" (Table S2: OH and NO2 wall sources, HNO3/NO2/O3 wall
uptake, N2O5 + H2O → wall-HNO3). Azzi et al. 2010 (SAPRC-07 vs CSIRO chamber) describe the lamps as
blacklight-blue, ~350–390 nm, giving "lower ozone photolysis and a higher NO2 photolysis rate" [S].

**ILMARI referee exchange, Carter 2014 [V]** — the clearest statement of why class B fails for lamps:
"It shows no intensity in the UV region from 300 to about 340 nm. Although this spectrum is sufficient
for photolyzing NO2, it will result in abnormally low photolysis rates relative to atmospheric for other
photolyses that are potentially important radical sources, such as photolysis of O3 forming O1D … and
photolysis of aldehydes and other photoreactive compounds whose action spectra decrease rapidly near the
300 nm cutoff for ground-level sunlight. … it is important that they have a spectrum similar to sunlight
in the important 300-400 nm region."

**OFR, Peng et al. 2015 [V]** — "The range of photon fluxes studied at 185 and 254 nm are 10^11–10^14 and
4.2 × 10^13–8.5 × 10^15 photons cm-2 s-1, respectively, based on the characterization of Li et al.
(2015)." The model does "not assess the impact of … Hg lamp emission variations in space and time,
variable UV fluxes at different points in the reactor". Line-flux × σ is class A degenerate to two
wavelengths; the absolute anchor is the O3/SO2/OH-exposure calibration of Li et al. 2015 rather than NO2.

---

## 5. What the public MCM box models offer

**AtChem2 (Sommariva et al. 2020; manual v1.2) [V]**

- Default: `j = l·(cos SZA)^m·exp(−n·sec SZA)·τ`, "where … τ is a transmission factor. … The transmission
  factor τ can be used to account for the loss of natural or artificial light in some environmental
  chambers caused, for example, by the transmittance of the chamber walls (by default, τ = 1). In AtChem2,
  the user can customise the photolysis rates parameterisation by providing an alternative file to
  replace the values of l, m, n and τ provided by the MCM."
- "A correction factor (JFAC) can be used to account for the difference between the photolysis rates,
  which are calculated by the model under clear-sky conditions, and the measured photolysis rates …
  JFAC = jmeas / jcalc (4), where jmeas and jcalc are the measured and calculated (with the MCM
  parameterisation) photolysis rates for the reference species, usually NO2."
- Manual §3.5.1: "The typical scenario for constant photolysis rates is a lamp or solar simulator in an
  environmental chamber. All the photolysis rates in the chemical mechanism must be given a fixed value
  in the photolysisConstant.config file, otherwise they are automatically set to zero … This approach
  allows the user to model individual photolytic processes and/or to account for lamps that emit only in
  certain spectral windows." Constrained J go in `photolysisConstrained.config`; "JFAC is not applied to
  constant or constrained photolysis rates"; `ROOF = CLOSED` zeroes every J (lamps off).
- So AtChem2 does **not** compute lamp J itself: the user must bring class-A values from elsewhere, or
  accept class B at a real SZA.

**F0AM (Wolfe et al. 2016) [V]**

- "F0AM provides three options for calculating J values: bottom-up, MCM, and hybrid. In the 'bottom-up'
  method, J values are calculated by integrating the product of a user-specified actinic flux spectrum
  with literature-derived cross sections and quantum yields. Cross sections and quantum yields are taken
  from the latest IUPAC … and JPL … recommendations … This option is most useful when simulating
  photochemical chamber experiments with non-solar light sources."
- "Any J value can be constrained to observations via direct input. It is also possible to specify a
  scaling factor for all parameterized J values. Typically this scaling is taken as the ratio of an
  observed photolysis frequency to its model-calculated value."
- Chamber example (§3.1): bottom-up with the Caltech 350 nm bulb spectrum.

**KinSim / OFR models (Peng & Jimenez lineage)** — explicit per-line photon fluxes × σ (class A at two
wavelengths); the absolute flux is the calibrated quantity (Peng et al. 2015 [V]).

---

## 6. Recommendation for a general-purpose chamber code (JUICE)

Given the history above, a chamber code that wants to be defensible against both SAPRC and MCM practice
should offer, in this order of preference, and record in the output which one was used:

1. **Class A as the primary lamp mode** — store per chamber profile a *relative* lamp spectrum F(λ)
   (spectroradiometer, or the canonical Carter blacklight spectrum when none was measured), compute
   J_i = j(NO2)_meas · ∫σ_i φ_i F dλ / ∫σ_NO2 φ_NO2 F dλ for every MCM J-number with IUPAC/JPL σ, φ, and
   anchor on an actinometric j(NO2). This is what SAPRC does, what the MCM's own chamber evaluations
   (Pinho/Jenkin via Carter's software; Metzger; Platt) did, and what F0AM recommends for non-solar
   light. Make the j(NO2) anchor time-dependent (lamp warm-up, ageing, number of lamps — Carter 2005,
   Winiberg 2016).
2. **Class D override** — any J measured in the chamber (spectroradiometer or actinometry, e.g. j(O1D),
   j(HONO), j(HCHO), j(X) for a studied compound) replaces the computed value; accept j(X)/j(NO2) ratios
   as input, since that is how CESAM and EUPHORE publish them.
3. **Class B only for solar-like sources, and labelled as such** — MCM l,m,n at the real SZA (outdoor) or
   at an effective SZA (Xe/Ar arc with Pyrex filter, LED solar simulators), scaled by j(NO2)_meas /
   j(NO2)_calc, optionally with a second species-dependent factor as in Bloss 2005a. The profile should
   carry a flag "spectrum solar-like: yes/no" and the code should refuse (or warn loudly) to apply class B
   to blacklights or Hg lamps, quoting Carter's 300–340 nm argument. Report j(O1D)/j(NO2) of the chosen
   source next to the clear-sky value so the user sees the distortion (CESAM Table 3 practice).
4. **Class C as a fallback with provenance** — ship canonical ratio tables per lamp type (from
   `lamp_photolysis_literature.md`) for users without a spectrum, clearly marked as imported.
5. **Class E as a diagnostic, never silent** — allow fitting a single J (or a global factor) to a tracer
   decay (Winiberg's j(Cl2), Wang's j(NO2) from NOx-air runs) but write the fitted value and the data it
   was fitted to into the run log.
6. **Chamber auxiliary mechanism in Carter form, parameterised relative to j(NO2)** — `wall + hν → HONO`
   at R_N·j(NO2) (RN in ppb), light-induced NO2 off-gassing, HCHO off-gassing, NO2 → HONO dark
   conversion, N2O5 wall/hydrolysis loss, O3/HNO3 wall loss, background OH reactivity as a CO-like
   species, initial HONO; with the SAPHIR/Rohrer and EUPHORE/Zádor RH- and T-dependent forms as
   alternatives. Every number in a chamber profile should cite its characterisation run or literature
   source, as Metzger 2008 Table 2 and Bloss 2005a Table 2 do, with lower/upper bounds for sensitivity
   runs.
7. **Never use the MCM l,m,n at a fictitious SZA to represent lamps without saying so.** No study in this
   corpus did it; AtChem2 only does it implicitly if misconfigured.

---

## 7. Not retrieved / unverified

- Pinho et al. 2005, 2006, 2007 and the 2006 Aveiro PhD thesis: Elsevier/Springer paywalled; the Viseu
  repository copy is behind a bot check (HTTP 000 / Anubis). Quotes in §3.4 are search snippets [S].
- Hynes et al. 2005: paywalled; described via Ginnebaugh et al. 2010 SI and Azzi et al. 2010 snippets.
- Zádor et al. 2006 (J. Atmos. Chem.) and Zádor et al. 2005 (AE, ethene uncertainty): paywalled;
  Semantic Scholar rate-limited; the HONO-source form is known only from secondary descriptions.
- Fuchs et al. 2013 Nat. Geosci.: paywalled; the identical procedure is quoted from Fuchs et al. 2014 ACP.
- Paulsen et al. 2005 ES&T; Rickard et al. 2010 AE; Carter & Lurmann 1991; Carter et al. 1982 (IJCK 14,
  1071 — this is the paper the task calls "Carter et al. 1984"); Carter, Winer, Darnall, Pitts 1979 ES&T:
  not retrieved; roles inferred from citing texts.
- "Hynes et al. 2005 ACP 5, 7" does not exist (ACP returns 404 for that page); the CSIRO paper is AE 39,
  7251.
- CSIRO 2018 new-chamber paper (Environ. Chem. EN18005): server timed out twice.
- The UCR Carter web server (intra.cert.ucr.edu) was unreachable; the SAPRC-99 documentation, the
  SAPRC-07 report and the UCR EPA chamber report were read from the Internet Archive copies.

## References (DOIs)

Azzi et al. 2010 AE 44, 1707 (10.1016/j.atmosenv.2010.01.049 [?]); Bloss et al. 2005a 10.5194/acp-5-623-2005;
Bloss et al. 2005b 10.5194/acp-5-641-2005; Bohn et al. 2005 10.5194/acp-5-493-2005; Boyd et al. 2015
10.5194/acp-15-7497-2015; Carter et al. 1982 10.1002/kin.550141003; Carter, Luo, Malkina, Pierce 1995 ARB
92-345; Carter 2000 SAPRC-99 documentation (CARB 92-329/95-308); Carter et al. 2005 10.1016/j.atmosenv.2005.08.040;
Carter 2010 10.1016/j.atmosenv.2010.01.026; Fuchs et al. 2013 10.1038/ngeo1964; Fuchs et al. 2014
10.5194/acp-14-7895-2014; Glowacki et al. 2007 10.5194/acp-7-5371-2007; Hynes et al. 2005
10.1016/j.atmosenv.2005.09.005; Jenkin et al. 2015 10.5194/acp-15-11433-2015; Leskinen et al. 2015
10.5194/amt-8-2267-2015; Malkin et al. 2010 10.5194/acp-10-1441-2010; Metzger et al. 2008
10.5194/acp-8-6453-2008; Nah et al. 2016 10.5194/acp-16-9361-2016; Paulsen et al. 2005 10.1021/es0489137;
Peng et al. 2015 10.5194/amt-8-4863-2015; Picquet-Varrault et al. 2020 10.5194/acp-20-487-2020; Pinho et al.
2005 10.1016/j.atmosenv.2004.11.014; Pinho et al. 2006 10.1007/s10874-006-9025-y; Pinho et al. 2007
10.1007/s10874-007-9071-0; Platt et al. 2013 10.5194/acp-13-9141-2013; Rohrer et al. 2005
10.5194/acp-5-2189-2005; Schwantes et al. 2017 10.5194/acp-17-3453-2017; Sommariva et al. 2020
10.5194/gmd-13-169-2020; Wagner et al. 2003 10.5194/acp-3-89-2003; Wang et al. 2011 10.5194/amt-4-2465-2011;
Winiberg et al. 2015 10.5194/amt-8-523-2015; Winiberg et al. 2016 10.5194/acp-16-4023-2016; Wolfe et al.
2016 10.5194/gmd-9-3309-2016; Wyche et al. 2009 10.5194/acp-9-635-2009; Zádor et al. 2006
10.1007/s10874-006-9033-y.
