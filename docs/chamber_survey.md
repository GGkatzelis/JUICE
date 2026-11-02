# Survey of atmospheric simulation chambers (Europe, USA) for generalising JUICE

Compiled 2026-11-02 from web search and page fetches. Rule followed: a value is given only when a fetched or
snippet source states it; otherwise "not found". "[unv]" = taken from a search-result snippet / secondary listing, not
verified in the primary paper. Many primary papers (Springer book chapters, some ACP PDFs, EUROCHAMP/ACTRIS
detail pages) could not be read in full, so "not found" often means "not found in what I could access", not "does not exist".
Latitude/longitude are given only where a page states them or the chamber is outdoor; "(approx.)" marks my own rough knowledge, not a fetched source.

Main hub sources used throughout:
- EUROCHAMP chamber list: https://www.eurochamp.org/simulation-chambers (17 chambers listed; per-chamber pages e.g. .../SAPHIR)
- ICARUS chamber/experiment database (UC Davis): https://icarus.ucdavis.edu (chamber records /chamber/<id>)
- ACTRIS National Facilities labelling: https://nflabelling.actris.eu/facility/<id>
- Doussin et al. (eds.), A Practical Guide to Atmospheric Simulation Chambers, Springer 2023, doi 10.1007/978-3-031-22277-1 (open access; chapter 2 "Physical and Chemical Characterization of the Chamber"; full text not retrievable here)
- Review: "Application of smog chambers in atmospheric process studies", PMC8826053 (tables S1/S2/S5 hold the chamber and wall-loss lists; supplement not accessible)
- Review: Atmosphere (Jan 2024) covering 18 US, 24 European, 8 Asian chambers (doi not retrieved; page blocked, 403)

## 1. Comparison table

jNO2 in 1/s. "nf" = not found. T control: Y = actively controlled, N = ambient.

| Chamber | Institution (country) | Volume m3 | Material | Light source | Typical jNO2 (1/s) | T control | Data public? | Main reference |
|---|---|---|---|---|---|---|---|---|
| SAPHIR | FZ Juelich (DE) | 270 (S/V ~1 m-1) | double-wall FEP film | sunlight (shutter) | ambient, from outside actinic flux; max ~ 8e-3 (clear summer, my estimate) | N (ambient, ~10-25 C) | EUROCHAMP DC [unv format] | Rohrer 2005 ACP; Bohn 2005 ACP |
| SAPHIR-PLUS | FZ Juelich (DE) | 9.32 (gas exchange vol.) | Teflon-type, plant chamber under SAPHIR | 15 LED panels | nf | nf | nf | Hohaus 2016 AMT 9, 1247 |
| AIDA | KIT Karlsruhe (DE) | 84.5 (+3.7 prep.) | aluminium (AlMg3) | ~3000 LEDs, 305-528 nm, 1400 W (solar simulator, since 2022) | 1.58e-3 +/- 0.21e-3 (213-295 K mean) | Y, -90 to +60 C; 0-100 % RH | EUROCHAMP DC [unv] | Wagner 2003 / Vogel; Ruescher/Wagner AMT 15, 1795 (2022) |
| CESAM | LISA / CNRS Creteil (FR) | 4.2 | stainless steel | 3 x 4 kW Xe arc lamps | 4.2e-3 (also jO1D 1.4e-5) | Y (design: variable T and p); tested 290-297 K | EUROCHAMP DC | Wang 2011 AMT 4, 2465 |
| HELIOS | ICARE / CNRS Orleans (FR) | 90 | FEP film, hemispherical | sunlight + movable box of Xe lamps | nf | nf | nf | McGillen/Daele (AGU 2016 abstract) |
| EUPHORE | CEAM Valencia (ES) | 2 x ~200 | FEP film 0.13 mm, half-spherical, Al floor | sunlight (>80 % transmission 280-640 nm) | nf (O3 wall loss 3.0e-6 1/s via a model paper) | floor cooling, no range found | EUROCHAMP DC | Becker 1996 (not retrieved) |
| ILMARI | Univ. Eastern Finland, Kuopio (FI) | 29 | Teflon (FEP) | up to 160 blacklights, 350 or 365 nm; 29.7 W/m2 max | 0 - 0.62 min-1 (= 0 - 1.0e-2 1/s) | Y, uniform, 25 +/- 1 C | nf | Leskinen 2015 AMT 8, 2267 |
| IASC (CRAC) | Univ. College Cork (IE) | 27 (4.5 x 3 x 2 m) | FEP foil | ~140 UV lamps (UV-A / UV-B) | nf | Y (enclosure; range nf) | nf | none found |
| CERNESIM / ESC-Q-UAIC | Univ. Al. I. Cuza, Iasi (RO) | ~0.78 | quartz tubes (3 x, 0.48 m dia., 4.2 m) | actinic + blacklight lamps | nf | nf | nf | nf |
| PACS-C3 | PSI Villigen (CH) | 27 (also 9 m3 mobile, 9 m3 cool) | FEP bag | Xe arc lamps (4, 16 kW total) + blacklights [unv] | 2.0e-3 | Y, 15-30 C (cool chamber 10-30 C) | EUROCHAMP DC [unv] | Paulsen 2005 ES&T 39, 2668 |
| MAC | Univ. Manchester (UK) | 18 | FEP bag | 2 x 6 kW Xe arc + 112 halogen bulbs | 2.25e-3 (actinometry, new bag); 1.5e-3 (radiometry) | Y, 15-35 C; RH 25-80 % | on request | Shao 2022 AMT 15, 539 |
| CHARME | LPCA, Univ. Littoral, Dunkerque (FR) | 9.2 (S/V 2.4 m-1) | electropolished 304L steel | 4 x 5 kW Xe lamps [unv] | nf | Y, -20 to +50 C; p 0.15-1030 hPa | nf | nf |
| LPCA-ONE | LPCA Dunkerque (FR) | ~8 | PMMA cube | 10 fluorescent tubes | nf | nf | nf | nf |
| ChAMBRe | INFN / Univ. Genoa (IT) | ~2.2 | stainless steel | none (aerosol/bio-aerosol; UV not a feature) | n/a | nf | nf | Massabo 2018 AMT 11, 5885 |
| LEAK | TROPOS Leipzig (DE) | 19 (cylinder) | nf | nf | nf | RH up to 80 % | nf | nf |
| QUAREC | Univ. Wuppertal (DE) | ~1.1 (0.47 m dia. x 2 x 3 m) | quartz, Al flanges | 32 Philips TL05 (300-480 nm) + 32 TUV 40 W low-pressure Hg (254 nm) | nf | nf | nf | Barnes/Becker (Wuppertal) - not retrieved |
| FORTH-ASC | FORTH / Univ. Patras (GR) | 10-15 (in 30 m3 room); 1 m3 combustion chamber | Teflon | UV fluorescent lamps on Al walls | nf | room T controlled | nf | nf |
| HIRAC | Univ. Leeds (UK) | 2 | stainless steel | UV lamps (EUROCHAMP list silent) | nf | nf | nf | Glowacki 2007 (not retrieved) |
| MICC / RvG-ASIC | Manchester / UEA (UK) | fall tube 10 m x 1 m dia. | - | - | - | down to -55 C, 50 mbar (MICC) | - | - |
| ISAC | CNRS Paris (FR) | ~2 | Teflon in frame | UV and visible sources | nf | nf | nf | nf |
| Lund, Gothenburg, Bayreuth, Mainz MPIC | - | not found | - | - | - | - | - | no public record found in this session |
| Caltech (2 chambers) | Caltech (US) | 2 x ~18-28 (records: 18, 22, 27) | FEP 2 mil | UV blacklights (354 nm) + Sunshine fluorescents (300-900 nm) | nf | Y, 12-40 C (site); ICARUS 20-45 C | ICARUS | Cocker 2001 ES&T; Keywood 2004 |
| UCR CE-CERT (SAPRC) | UC Riverside (US) | 2 x 90 (Carter); 84.2 (ICARUS "EPA chamber") | FEP 2 mil | Ar arc lamp + 80 blacklights (350BL) | nf (blacklight and arc, adjustable) | Y, 5-45 C; RH 0.01-80 % | Carter website; ICARUS | Carter 2005 Atmos. Environ. |
| CMU CAPS | Carnegie Mellon (US) | ~10-12 [unv] | FEP | UV blacklights | nf | Y (room) | ICARUS | Hennigan 2011 ACP; Presto 2005 |
| GTEC (Ng group) | Georgia Tech (US) | 2 x 12 (12.6; S ~32.5 m2) | FEP | UVA blacklights (300-400 nm, peak 354) + sunshine fluor.; or sunlight | light-flux characterisations in ICARUS | Y, 4-40 +/- 0.5 C | ICARUS | none single paper found |
| CU Boulder CAChRF | CU Boulder (US) | 20 (1-White) | Teflon | UVA blacklights + visible fluor. | nf | Y, 4-40 C | ICARUS | Jimenez group |
| CU Boulder Ziemann ("Janet", "Roger") | CU Boulder (US) | 2 x 8 | Teflon | UVA blacklights | nf | N | ICARUS | nf |
| NCAR | NCAR Boulder (US) | 10 (3rd floor Teflon) | Teflon | UVA blacklights | nf | 22-35 C | ICARUS | Orlando/Tyndall (not retrieved) |
| UC Irvine | UC Irvine (US) | 5.6 | Teflon | UVB lamps | nf | 22-29 C | ICARUS | Nizkorodov group |
| UC Davis | UC Davis (US) | 12 | Teflon | UVA black lamps | nf | Y (range nf) | ICARUS | Nguyen group |
| UT Austin (LHR) | Univ. Texas Austin (US) | 8.92 | Teflon | UV-A blacklight | nf | Y, 5-50 C | ICARUS | Hildebrandt Ruiz group |
| Harvard | Harvard (US) | 4.7 | PFA Teflon | lamps with reflective Al enclosure (type nf) | nf | Y | nf | JoVE 55685 |
| UNC outdoor | Univ. North Carolina (US) | nf (historic) | FEP | sunlight | nf | N | historic only | Jeffries et al., Environ. Chem. 10, 349 |
| EPA RTP / mobile | US EPA (US) | 14.3 (mobile chamber) [unv] | Teflon-lined | ~120 UV bulbs [unv] | nf | nf | nf | nf |
| PNNL, MIT, Purdue, Texas A&M | - | not found | - | - | - | - | - | no public record found in this session |

## 2. Per-chamber blocks

### SAPHIR and SAPHIR-PLUS (FZ Juelich, DE)
- Site: Juelich, approx. 50.9 N, 6.4 E (approx., not from a fetched page). Sunlight used.
- 270 m3 (cylinder 5 m dia. x 18 m; other papers say 20 m long, 320 m2 surface); S/V ~1 m-1; double-wall FEP film in steel frame.
- Light: natural sunlight, shutter for dark runs. jNO2 from actinic flux measured outside, corrected by chemical actinometry
  (Bohn 2005: calibration factor 1.34 +/- 0.10, jNO2 accuracy 13 %).
- T: ambient, typically 10-25 C, no active control. Overpressure ~35 Pa keeps ambient air out.
- Dilution: from measured replenishment flow; ~4 % per hour (~1.1e-5 1/s) [unv]; one paper quotes 0.034 1/h for particles.
- Wall effects: aerosol loss ~2e-5 1/s (lifetime ~14 h); O3 loss in clean air < 3e-6 1/s; photolytic HONO source
  parameterised from T, solar radiation, humidity (Rohrer 2005: depends on humidity, light, T; fit reproduces HONO to 10 %; "a few
  hundred pptv/h" HONO, HCHO, acetone formed); background OH reactivity parameterised as a CO-like reactant, typically < 1.5 1/s.
  Dilution experiments showed no measurable MACR wall loss [unv].
- Data: EUROCHAMP Data Centre (format see section 4).
- SAPHIR-PLUS: 9.32 m3 flow-through plant chamber below SAPHIR, Teflon-type surfaces, 15 LED panels (Hohaus 2016). SAPHIR-STAR: no public record found.
- Refs: Rohrer et al. 2005, ACP 5, 2189; Bohn et al. 2005, ACP 5, 493; Hohaus et al. 2016, AMT 9, 1247.
- URLs: https://www.eurochamp.org/simulation-chambers/SAPHIR ; https://acp.copernicus.org/articles/5/2189/2005/ ; https://amt.copernicus.org/articles/9/1247/2016/

### AIDA (KIT, DE)
- 84.5 m3 aluminium vessel (+3.7 m3 prep chamber), T +50 to -90 C, 0-100 % RH, p controlled. Not transparent to sun; new
  LED solar simulator (~3000 LEDs, 11 peak wavelengths 305-528 nm, 1400 W): jNO2 = (1.58 +/- 0.21)e-3 1/s, about one third of
  clear-sky max (3.4e-3 winter, 6.0e-3 summer at Karlsruhe). Wall loss quantified only for 2,3-pentanedione; too high at 213 K to model.
- Refs: AMT 15, 1795 (2022) (LED simulator). URLs: https://amt.copernicus.org/articles/15/1795/2022/

### CESAM (LISA, FR)
- 4.2 m3 stainless steel; 3 x 4 kW Xe arc lamps filtered to solar-like spectrum; jNO2 = 4.2e-3 1/s, jO1D = 1.4e-5 1/s.
  Evacuable, variable T and p, characterised 290-297 K dry. Particle lifetime 10 h to 4 days. Authors built an auxiliary NOy
  wall mechanism for MCM (rate values not retrieved); HONO wall source scales with jNO2 (statement in a reviewer comment, not the paper).
- Ref: Wang et al. 2011, AMT 4, 2465. URLs: https://amt.copernicus.org/articles/4/2465/2011/ ; https://cesam.cnrs.fr

### HELIOS (ICARE Orleans, FR)
- Outdoor hemispherical 90 m3 FEP (>90 % solar transmission) on ICARE building roof, 47.838 N, 1.944 E (ACTRIS page). Movable box with Xe lamps allows
  artificial light and dark runs. jNO2, wall losses: not found. Ref: AGU 2016 abstract A33E-0290 (McGillen). URLs: https://icare.cnrs.fr/en/en-research/en-facilities/en-ra-helios/ ; https://nflabelling.actris.eu/facility/44

### EUPHORE (CEAM Valencia, ES)
- Two ~200 m3 half-spherical FEP (0.13 mm) bags, aluminium floor with refrigeration; sunlight, >80 % transmission 280-640 nm; stable under
  100-200 Pa overpressure. Site approx. 39.5 N, 0.4 W (Paterna, approx.). jNO2, mixing, dilution: not on CEAM page.
- Wall numbers via a UNC wall-model mechanism [unv]: O3 loss 3.0e-6 1/s; HNO3 first-order 8.2e-5 1/s; wall HONO source ~ 1.8e-3 x jNO2 scaling.
  A published parameterisation has HONO production W = a x jNO2 x exp(-T0/T) for dry air, humid expression with RH exponent 0.36 (2-15 % RH) [unv].
  S/V "close to one" [unv].
- URLs: https://www.ceam.es/WWWEUPHORE/paginas/what_is_euphore.htm

### ILMARI (UEF Kuopio, FI)
- 29 m3 FEP, up to 160 blacklights (350 or 365 nm), max 29.7 W/m2, jNO2 0-0.62 min-1 adjustable, T 25 +/- 1 C, kept in overpressure with moving top frame.
  Particle wall loss 0.088 1/h (polydisperse), 0.080 (50 nm), 0.045 (100 nm), 0.040 (200 nm). O3/NO2 wall loss: not in abstract.
- Ref: Leskinen et al. 2015, AMT 8, 2267. URL: https://amt.copernicus.org/articles/8/2267/2015/

### IASC / CRAC (UCC Cork, IE)
- 27 m3 (4.5 x 3 x 2 m) FEP in air-conditioned housing, ~140 UV lamps (UV-A and UV-B) for day/night cycles. Other UCC page says ~25 m3. 6-h experiments offered via ACTRIS. No wall numbers found.
- URLs: https://ucc.ie/en/crac/facilities/iasc ; https://meetingorganizer.copernicus.org/EGU24/EGU24-19774.html

### CERNESIM / ESC-Q-UAIC (Iasi, RO)
- ~780 L, three quartz tubes in series (0.48 m dia., 4.2 m total), actinic + blacklight lamps. All further numbers: not found. URL: https://www.eurochamp.org/simulation-chambers

### PACS-C3 (PSI, CH)
- 27 m3 FEP bag in T-controlled enclosure (15-30 C), 4 Xe arc lamps (16 kW) [unv] plus blacklights; jNO2 2.0e-3 1/s (AMT 15, 1795 table).
  Particle wall loss 0.209 +/- 0.018 1/h (number), 0.139 +/- 0.070 1/h (volume) (Paulsen 2005). Gas wall loss: not found.
- Ref: Paulsen et al. 2005, ES&T 39, 2668.

### MAC (Manchester, UK)
- 18 m3 FEP bag, 2 x 6 kW Xe arc (Osram XBO 6000 W/HSLA OFR) + 112 Solux 50 W halogens; T 15-35 C, RH 25-80 % (characterised at 25 C, 50 %).
- jNO2 2.25e-3 1/s (NO2 actinometry, new bag), 1.5e-3 1/s (spectroradiometry); clear-sky Manchester summer ~7e-3 1/s. Chamber ageing may lower jNO2 (not tested).
- Dark wall loss: NO2 9.4 +/- 7.4 e-7 1/s; O3 2.09 +/- 0.97 e-6 1/s. Particle losses size-resolved (40-600 nm), depend on RH, mixing, bag history. HONO not characterised.
- Mixing time: a few minutes (gases), ~2.5 min (seed particles). Data: "available on request".
- Ref: Shao et al. 2022, AMT 15, 539. URL: https://amt.copernicus.org/articles/15/539/2022/

### CHARME and LPCA-ONE (Dunkerque, FR)
- CHARME: cylindrical electropolished 304L steel, 9.2 m3, S/V 2.4 m-1, -20 to +50 C, 0.15-1030 hPa, 4 x 5 kW Xe lamps (single source: LPCA page). LPCA-ONE: ~8 m3 PMMA cube, 10 fluorescent tubes.
  Coupled to MULTICHARME multipass IR/THz cell. URL: https://lpca.univ-littoral.fr/?p=8993

### ChAMBRe (Genoa, IT)
- ~2.2 m3 stainless steel vertical cylinder, evacuable to ~2e-2 mbar in ~10 min; bio-aerosol focus; particle lifetime 2-10 h; reproducibility ~10 %. No photochemistry lamps. Ref: Massabo et al. 2018, AMT 11, 5885.

### LEAK-LACIS (TROPOS, DE)
- LEAK 19 m3 cylindrical, RH up to 80 %. Nothing else found. Not a photochemistry chamber (no lamps listed).

### QUAREC (Wuppertal, DE) - 254 nm and 300-480 nm
- Two quartz cylinders (47 cm dia., 3 m each) joined by a flange, ~1.1 m3 (computed). 64 lamp slots, switched in pairs: 32 Philips TL05 superactinic (300-480 nm) + 32 Philips TUV 40 W low-pressure Hg (peak 254 nm). jNO2: not found. T, wall loss: not found.
- URLs: https://quarec.uni-wuppertal.de/en/quarec/infrastructure/quarec-asc/

### FORTH-ASC (Patras, GR)
- 10 m3 Teflon reactor (EUROCHAMP: 10-15 m3 in 30 m3 room), surrounded by aluminium walls with UV fluorescent lamps, T-controlled room; plus 1 m3 combustion chamber. ACTRIS status "planned". Numbers: not found.

### Other EUROCHAMP/ACTRIS chambers (brief)
- HIRAC (Leeds): 2 m3 stainless steel, cylindrical. ISAC (CNRS): ~2 m3 Teflon. MICC (Manchester): 10 m fall tube, to -55 C, 50 mbar. RvG-ASIC (UEA): ocean-atmosphere chamber. CLOUD (CERN): 26 m3 stainless steel, Hg-Xe UV lamps, 207-310 K. Lund, Gothenburg, Bayreuth, Mainz MPIC: no record found.

### Caltech (US)
- Two flexible FEP (2 mil) chambers, ~18 m3 each on the group page (ICARUS records 27, 22 m3 for earlier bags; Cocker 2001 described 28 m3). Enclosure 12-40 C. Blacklights peak 354 nm + Sunshine fluorescent lamps (300-900 nm).
  Particle wall loss ~2-4 per day (Cocker 2001, via 2023 citation). Gas wall-loss numbers (O3, NO2, HNO3): not found; large literature on vapor wall loss (Zhang 2014).
- URLs: https://nggroup.caltech.edu/research/chamber-facility ; https://icarus.ucdavis.edu/chamber/548

### UC Riverside CE-CERT / SAPRC (US)
- Two collapsible 90 m3 FEP reactors in a 450 m3 insulated, reflective enclosure; Ar arc lamp (200 kW class) and 80 blacklights (115 W Sylvania 350BL); 5-45 C; RH 0.01-80 %. ICARUS "EPA chamber": 84.2 m3, 126 m2.
- Wall effects: early characterisation: radical source ~40-80 ppt/h NOx offgassing equivalent; O3 wall loss 0.8 %/h (~2.2e-6 1/s); chamber radical source represented as HONO offgassing proportional to jNO2 with ratio ~5-12.5 ppt HONO per unit (Carter 2004 slides; unit of the ratio not clear in snippet). Carter notes the original blacklights lack 300-340 nm output; newer lamps cover 310-400 nm.
- Refs: Carter et al. 2005, Atmos. Environ. 39 (not directly retrieved); reports at https://intra.cert.ucr.edu/~carter/pubs/chamrpt.pdf

### Carnegie Mellon (US)
- Teflon bag(s) in T-controlled room, UV blacklights; volume 7 m3 (Hennigan 2011) / ~10-12 m3 (other sources) - uncertain. Particle wall-loss rates listed in review supplement (not retrieved). ICARUS contributor (CAPS).

### Georgia Tech GTEC (US)
- Two 12 m3 FEP bags (12.6 m3, 32.5 m2 surface => S/V ~2.6 m-1), enclosure 4-40 +/- 0.5 C, blacklights (300-400 nm, peak 354) + sunshine fluorescent lamps; sunlight possible per ICARUS. ICARUS holds light-flux characterisation runs (Aug 2019, Jul 2020), particle and vapor loss runs. 2013 news item said 10 m3 (earlier configuration).
- URL: https://ng.chbe.gatech.edu/chamber-facility ; https://icarus.ucdavis.edu/chamber/542

### Colorado (CU Boulder), NCAR, UC Irvine, UC Davis, UT Austin (US)
- From ICARUS records: CU 1-White 20 m3 Teflon, UVA + visible fluorescent, 4-40 C (Jimenez); Ziemann "Janet"/"Roger" 8 m3, no T control; NCAR 10 m3 Teflon, UVA, 22-35 C; UCI 5.6 m3 Teflon, UVB lamps, 22-29 C; UC Davis 12 m3 Teflon, UVA; UT Austin 8.92 m3 Teflon, UV-A, 5-50 C. No jNO2 or gas wall-loss numbers found.

### Harvard (US)
- 4.7 m3 PFA Teflon bag in T-controlled enclosure with reflective Al sheets (JoVE 55685). Lamp type: not found.

### UNC outdoor chambers (US)
- Established 1972 (Jeffries, Kamens, Sexton); outdoor Teflon/FEP bags with sunlight; early 2.8 m3 BOM chamber (Hanst). Earlier Jeffries chambers: 2 x 156 m3 FEP (from the PMC review). Wall model mechanism ("UNC outdoor") used in the SAPRC mechanism with jNO2-scaled HONO source. Chamber is not described as currently active.

### EPA RTP, PNNL, MIT, Purdue, Texas A&M
- EPA mobile chamber ~14.3 m3 Teflon-lined [unv]. Others: no public chamber description found.

## 3. Light-source tally (see summary in the reply)
- Sunlight (outdoor or sun-capable): SAPHIR, HELIOS, EUPHORE, UNC outdoor, GTEC (optional).
- Xenon arc / solar simulator: CESAM, CHARME, MAC, PACS-C3 (with blacklights), HELIOS (box), UCR (Ar arc), AIDA (LED simulator instead).
- Blacklights 350-365 nm: ILMARI, Caltech, UCR, CMU, GTEC, CU, NCAR, UC Davis, UT Austin.
- 254 nm (low-pressure Hg, UVC): QUAREC only confirmed. CLOUD uses Hg-Xe 250-750 nm lamps. No chamber with a 185 nm lamp found.
- Fluorescent UV-A/UV-B: IASC (140 lamps), FORTH-ASC, LPCA-ONE, Nizkorodov UVB.

## 4. Data availability
- EUROCHAMP Data Centre (data.eurochamp.org, to merge into ACTRIS DC): free, registration required, DOIs available; legacy "EUROCHAMP Data Format" (EDF, human-readable ASCII; from an EGU 2008 abstract, may have changed).
- ICARUS (icarus.ucdavis.edu, mirror at NCAR GDEX): open, curated, searchable; download requires login in the web UI; file format not stated on fetched pages.
- Others: MAC data "on request"; UCR reports on the Carter website.

## 5. What a chamber profile must capture

Fields that vary across chambers (belong in a per-chamber config):

| Field | Unit | Typical range / examples |
|---|---|---|
| name, institution, lat/lon (sunlight only), altitude | - / deg | SAPHIR 50.9 N; HELIOS 47.84 N 1.94 E |
| volume V | m3 | 0.8 (CERNESIM) - 270 (SAPHIR); lamps chambers 2-30 |
| surface area A / S/V | m2, m-1 | S/V ~1 (SAPHIR) - 2.4 (CHARME) - 2.6 (GTEC); small steel/quartz 5-10 |
| wall material (FEP, steel, quartz, PMMA, Al) | enum | determines wall-loss set |
| volume mode | enum | fixed (steel/quartz), collapsible bag (V shrinks with sampling) |
| light type | enum | sunlight, Xe arc, blacklight 350/365, UVC 254, LED, fluorescent |
| jNO2 (constant or time series) | 1/s | 1e-3 to 1e-2 (lamps: 1.5e-3 MAC, 2e-3 PSI, 4.2e-3 CESAM, up to 1.0e-2 ILMARI; sun up to ~8e-3) |
| photolysis scaling for other species: ratios j(X)/jNO2 or spectrum (lamp spectral flux) | - | blacklights give a very different j(HCHO), j(O1D), j(HONO) ratio than sunlight; 254 nm lamps photolyse O3, H2O2 strongly |
| lamp on/off steps, dimmer fraction | - | 0-100 % (UCR, ILMARI 0-0.62 min-1) |
| temperature T and its control | K | ambient 283-298 (SAPHIR); 288-308 (PSI); 258-323 (CHARME); 183-333 (AIDA) |
| pressure | hPa | ~ambient (overpressure 35-200 Pa) vs evacuable (CESAM, CHARME 0.15 hPa) |
| RH / H2O | % | 0-100; important for wall HONO |
| dilution / replenishment rate kdil | 1/s | ~1e-5 (SAPHIR, ~4 %/h); bags with sampling losses 1e-5 to 1e-4; flow-through reactors much higher |
| mixing time | min | 2-10 (MAC few min); fans in steel chambers faster |
| gas wall loss k_O3, k_NO2, k_HNO3, k_organics | 1/s | O3 ~1e-6 to 3e-6 (Teflon); NO2 ~1e-6 (MAC) to 1e-4 1/min equivalents; steel: far larger (NO2 uptake up to 1e-2 1/s illuminated with ~10 % HONO yield in a 4.2 m3 steel chamber) |
| particle wall loss k_p(d) | 1/s or 1/h | 0.04-0.09 1/h (ILMARI, Teflon), 0.14-0.21 1/h (PSI), 2e-5 1/s (SAPHIR) |
| HONO wall source | pptv/h or ppt per jNO2 | SAPHIR: a few hundred pptv/h (sun, humid); EUPHORE ~1.8e-3 x jNO2 (UNC model); form: S = a * jNO2 * f(RH,T) |
| HCHO, acetone wall sources | pptv/h | SAPHIR parameterised; others nf |
| background OH reactivity | 1/s | SAPHIR < 1.5 (CO-like) |
| background / initial contaminant levels (NOx, O3, VOC) | ppb | lamp chambers 0.1-several ppb NOx |
| dark vs illuminated wall NO2 conversion | - | steel chambers need this; Teflon negligible |
| cleaning/flushing state, wall aging | - | MAC: aging may lower jNO2 |

Fields that are the same for everyone (need not be configurable): rate coefficients and temperature/pressure dependences from the chosen mechanism (MCM), absorption cross sections and quantum yields of species, gas constant and molar volume, ideal gas conversion ppb<->molecule cm-3, solar geometry (computed from lat/lon/time if sunlight), the first-order structure of dilution and wall loss (dC/dt = -k C), the HONO photolysis and OH-radical chemistry, and the input-file/time-grid conventions of the model.

## 6. Open gaps to close if time allows
- Springer book chapters (Doussin 2023, ch. 2 "Physical and Chemical Characterization", ch. on SAPHIR, CESAM, EUPHORE) hold the best per-chamber wall numbers; they could not be retrieved.
- Rohrer 2005 PDF is saved locally but no text extractor available (Python/pdftoppm missing).
- Verify CMU, Caltech, PNNL, MIT, Purdue, Texas A&M volumes; verify Lund, Gothenburg, Bayreuth, Mainz status.
