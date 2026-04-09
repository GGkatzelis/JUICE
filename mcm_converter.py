"""
mcm_converter.py — Python port of mcm331cv_FR.pro

Converts a raw MCM Master Chemical Mechanism file (downloaded from the York
website, e.g. mcm_3-3-1_facsimile_complete.fac) into the EASY-format include
files that JUICE consumes:

    <stem>.mech       — bare reaction lines (VOC chemistry)
    <stem>.rate       — k[...]=CONST(...) rate definitions
    <stem>.ro2.incl   — RO2 = CONST(0.0 + sp1 + sp2 + ...) peroxy radical pool
    <stem>.incl       — KMT/KFPAN/KBPAN/J<n> definitions for constants USED
                        by the mechanism (filtered list)

In SAPHIR mode (saphir=True), additionally produces:

    <stem>.dil.mech              — `species + DIL --> ` for every species
    <stem>.dil.rate              — `k[species + DIL -->] = CONST(DILUTE)`
    <stem>.saphir.easy_temp      — master EASY input file with all #includes
                                    and SAPHIR-specific boilerplate (HONO/HCHO
                                    chamber sources, walldep, bgreac, dilution)
    anorg.mcm331.dark.mech       — hardcoded inorganic dark reactions (37)
    anorg.mcm331.dark.rate       — corresponding rate constants
    anorg.mcm331.phot.mech       — hardcoded inorganic photolysis (8)
    anorg.mcm331.phot.rate       — corresponding photolysis frequencies

Original IDL routine: Matthias Karl 2002, updated by Cornelia Richter (2004),
Keding Lu (2011), Xin Li (2015-16), Franz Rohrer (2023). This Python port
follows the same algorithm; output is semantically equivalent (same
reactions, same rates, same RO2 pool) though formatting may differ slightly.

Usage:
    from mcm_converter import convert_mcm_file
    convert_mcm_file('ethanol_original.fac', output_dir='out/', saphir=True)
"""

import re
from pathlib import Path

# ─────────────────────────────────────────────────────────────────────────────
# Hardcoded constants from the IDL fktext array (verbatim from mcm331cv_FR.pro)
# These are the COMPRC/CONST/PHTFKT definitions emitted to <stem>.incl when
# the corresponding constant name is referenced anywhere in the mechanism.
# ─────────────────────────────────────────────────────────────────────────────

KNOWN_CONSTANTS = {
    'KFPAN':     'KFPAN   = COMPRC( 3.28D-28*M*(T/300)@(-6.87), 1.125D-11*(T/300)@(-1.105), 0.3)',
    'KBPAN':     'KBPAN   = COMPRC( 1.10D-05*M*EXP(-10100/T), 1.90D+17*EXP(-14100/T), 0.3 )',
    'KBPPN':     'KBPPN   = COMPRC( 1.7D-03*M*EXP(-11280/T), 8.3D+16*EXP(-13940/T), 0.36 )',
    'KMT01':     'KMT01   = COMPRC( 1.00D-31*M*(T/300)@(-1.6), 5.00D-11*(T/300)@(-0.3), 0.85 )',
    'KMT02':     'KMT02   = COMPRC( 1.30D-31*M*(T/300)@(-1.5), 2.30D-11*(T/300)@0.24, 0.6 )',
    'KMT03':     'KMT03   = COMPRC( 3.60D-30*M*(T/300)@(-4.1), 1.90D-12*(T/300)@0.2, 0.35 )',
    'KMT04':     'KMT04   = COMPRC( 1.30D-03*M*(T/300)@(-3.5)*EXP(-11000/T), 9.70D+14*(T/300)@0.1*EXP(-11080/T), 0.35)',
    'KMT05':     'KMT05   = CONST( 1.44D-13*(1+(M/4.2D+19)) )',
    'KMT06':     'KMT06   = CONST(1.+1.40D-21*H2O*EXP(2200/T))',
    'KMT07':     'KMT07   = COMPRC( 7.40D-31*M*(T/300)@(-2.4), 3.30D-11*(T/300)@(-0.3), 0.81 )',
    'KMT08':     'KMT08   = COMPRC( 3.20D-30*M*(T/300)@(-4.5), 3.0D-11, 0.41 )',
    'KMT09':     'KMT09   = COMPRC( 1.40D-31*M*(T/300)@(-3.1), 4.00D-12, 0.4 )',
    'KMT10':     'KMT10   = COMPRC( 4.10D-05*M*EXP(-10650/T), 6.00D+15*EXP(-11170/T), 0.4 )',
    'KMT11':     'KMT11   = CONST( 2.40D-14*EXP(460/T) + (6.50D-34*EXP(1335/T)*M)/(1+(6.50D-34*EXP(1335/T)*M/(2.70D-17*EXP(2199/T)))))',
    'KMT12':     'KMT12   = COMPRC( 2.5D-31*M*(T/300)@(-2.6), 2.0D-12, 0.53 )',
    'KMT13':     'KMT13   = COMPRC( 2.50D-30*M*(T/300)@(-5.5), 1.80D-11, 0.36 )',
    'KMT14':     'KMT14   = COMPRC( 9.00D-05*M*EXP(-9690/T), 1.10D+16*EXP(-10560/T), 0.36)',
    'KMT15':     'KMT15   = COMPRC( 8.60D-29*M*(T/300)@(-3.1), 9.00D-12*(T/300)@(-0.85), 0.48 )',
    'KMT16':     'KMT16   = COMPRC( 8.00D-27*M*(T/300)@(-3.5), 3.00D-11*(T/300)@(-1), 0.5 )',
    'KMT17':     'KMT17   =  COMPRC( 5.00D-30*M*(T/300)@(-1.5), 1.0D-12,0.17*EXP(-51/T)+EXP(-T/204))',
    'KMT18':     'KMT18  =  CONST(9.5D-39*O2*EXP(5270/T)/(1+7.5D-29*O2*EXP(5610/T)))',
    'KRO2NO3':   'KRO2NO3 = CONST( 2.30D-12 )',
    'KRO2NO':    'KRO2NO  = CONST(2.70D-12*EXP(360/T) )',
    'KRO2HO2':   'KRO2HO2 = CONST(2.91D-13*EXP(1300/T) )',
    'KAPHO2':    'KAPHO2  = CONST(5.20D-13*EXP(980/T))',
    'KAPNO':     'KAPNO   = CONST(7.50D-12*EXP(290/T))',
    'KNO3AL':    'KNO3AL  = CONST(1.40D-12*EXP(-1860/T))',
    'KDEC':      'KDEC    = CONST(1.00D+06) ',
    'KROPRIM':   'KROPRIM = CONST(2.50D-14*EXP(-300/T) )',
    'KROSEC':    'KROSEC  = CONST(2.50D-14*EXP(-300/T) )',
    'KCH3O2':    'KCH3O2 = CONST(1.03D-13*EXP(365/T))',
    'K298CH3O2': 'K298CH3O2 = CONST(3.5D-13)',
    'K14ISOM1':  'K14ISOM1 = CONST(3.00D7*EXP(-5300/T))',
    'J1':        'J1 = PHTFKT(6.073D-05,1.743,-0.474)',
    'J2':        'J2 = PHTFKT(4.775D-04,0.298,-0.080)',
    'J3':        'J3 = PHTFKT(1.041D-05,0.723,-0.279)',
    'J4':        'J4 = PHTFKT(1.165D-02,0.244,-0.267)',
    'J5':        'J5 = PHTFKT(2.485D-02,0.168,-0.108)',
    'J6':        'J6 = PHTFKT(1.747D-01,0.155,-0.125)',
    'J7':        'J7 = PHTFKT(2.644D-03,0.261,-0.288)',
    'J8':        'J8 = PHTFKT(9.312D-07,1.230,-0.307)',
    'J11':       'J11 = PHTFKT(4.642D-05,0.762,-0.353)',
    'J12':       'J12 = PHTFKT(6.853D-05,0.477,-0.323)',
    'J13':       'J13 = PHTFKT(7.344D-06,1.202,-0.417)',
    'J14':       'J14 = PHTFKT(2.879D-05,1.067,-0.358)',
    'J15':       'J15 = PHTFKT(2.792D-05,0.805,-0.338)',
    'J16':       'J16 = PHTFKT(1.675D-05,0.805,-0.338)',
    'J17':       'J17 = PHTFKT(7.914D-05,0.764,-0.364)',
    'J18':       'J18 = PHTFKT(1.482D-06,0.396,-0.298)',
    'J19':       'J19 = PHTFKT(1.482D-06,0.396,-0.298)',
    'J20':       'J20 = PHTFKT(7.600D-04,0.396,-0.298)',
    'J21':       'J21 = PHTFKT(7.992D-07,1.578,-0.271)',
    'J22':       'J22 = PHTFKT(5.804D-06,1.092,-0.377)',
    'J23':       'J23 = PHTFKT(2.4246D-06,0.395,-0.296)',
    'J24':       'J24 = PHTFKT(2.424D-06,0.395,-0.296)',
    'J31':       'J31 = PHTFKT(6.845D-05,0.130,-0.201)',
    'J32':       'J32 = PHTFKT(1.032D-05,0.130,-0.201)',
    'J33':       'J33 = PHTFKT(3.802D-05,0.644,-0.312)',
    'J34':       'J34 = PHTFKT(1.537D-04,0.170,-0.208)',
    'J35':       'J35 = PHTFKT(3.326D-04,0.148,-0.215)',
    'J41':       'J41 = PHTFKT(7.649D-06,0.682,-0.279)',
    'J51':       'J51 = PHTFKT(1.588D-06,1.154,-0.318)',
    'J52':       'J52 = PHTFKT(1.907D-06,1.244,-0.335)',
    'J53':       'J53 = PHTFKT(2.485D-06,1.196,-0.328)',
    'J54':       'J54 = PHTFKT(4.095D-06,1.111,-0.316)',
    'J55':       'J55 = PHTFKT(1.135D-05,0.974,-0.309)',
    'J56':       'J56 = PHTFKT(4.365D-05,1.089,-0.323)',
    'J61':       'J61 = PHTFKT(7.537D-04,0.499,-0.266)',
}

# Order matters: longer J names must be checked before shorter ones to avoid
# J1 matching inside J11. This list is also the canonical iteration order
# used when emitting the .incl file.
KNOWN_NAMES_ORDER = [
    'KFPAN', 'KBPAN', 'KBPPN',
    'KMT01', 'KMT02', 'KMT03', 'KMT04', 'KMT05', 'KMT06', 'KMT07',
    'KMT08', 'KMT09', 'KMT10', 'KMT11', 'KMT12', 'KMT13', 'KMT14',
    'KMT15', 'KMT16', 'KMT17', 'KMT18',
    'KRO2NO3', 'KRO2NO', 'KRO2HO2', 'KAPHO2', 'KAPNO', 'KNO3AL',
    'KDEC', 'KROPRIM', 'KROSEC',
    'KCH3O2', 'K298CH3O2', 'K14ISOM1',
    'J1', 'J2', 'J3', 'J4', 'J5', 'J6', 'J7', 'J8',
    'J11', 'J12', 'J13', 'J14', 'J15', 'J16', 'J17', 'J18', 'J19',
    'J20', 'J21', 'J22', 'J23', 'J24',
    'J31', 'J32', 'J33', 'J34', 'J35',
    'J41',
    'J51', 'J52', 'J53', 'J54', 'J55', 'J56',
    'J61',
]


# ─────────────────────────────────────────────────────────────────────────────
# Hardcoded inorganic chemistry (verbatim from IDL reactions[]/reacconst[] arrays)
# ─────────────────────────────────────────────────────────────────────────────

# Each entry is (mech_line, rate_def). Order matches IDL.
INORGANIC_DARK = [
    (' O  -->  O3 ',                  'k[O-->O3]=CONST(5.60D-34*O2*N2*(((T)/300)@-2.6) + 6.00D-34*O2*O2*(((T)/300)@-2.6))'),
    (' O + O3  -->   ',                'k[O+O3-->]=CONST(8.00D-12*EXP(-2060/(T)))'),
    (' O + NO  -->  NO2 ',             'k[O+NO-->NO2]=CONST( KMT01 )'),
    (' O + NO2  -->  NO ',             'k[O+NO2-->NO]=CONST(5.50D-12*EXP(188/(T)))'),
    (' O + NO2  -->  NO3 ',            'k[O+NO2-->NO3]=CONST( KMT02 )'),
    (' O1D  -->  O ',                  'k[O1D-->O]{SUM}=CONST(2.00D-11*N2*EXP(130/(T))+3.20D-11*O2*EXP(67/(T)) )'),
    (' NO + O3  -->  NO2 ',            'k[NO+O3-->NO2]=CONST(1.40D-12*EXP(-1310/(T)))'),
    (' NO2 + O3  -->  NO3 ',           'k[NO2+O3-->NO3]=CONST(1.40D-13*EXP(-2470/(T)))'),
    (' NO + NO  -->  NO2 + NO2 ',      'k[NO+NO-->NO2+NO2]=CONST(3.30D-39*EXP(530/(T))*O2)'),
    (' NO + NO3  -->  NO2 + NO2 ',     'k[NO+NO3-->NO2+NO2]=CONST(1.80D-11*EXP(110/(T)))'),
    (' NO2 + NO3  -->  NO + NO2 ',     'k[NO2+NO3-->NO+NO2]=CONST(4.50D-14*EXP(-1260/(T)))'),
    (' NO2 + NO3  -->  N2O5 ',         'k[NO2+NO3-->N2O5]=CONST( KMT03 )'),
    (' N2O5  -->  NO2 + NO3 ',         'k[N2O5-->NO2+NO3]=CONST( KMT04 )'),
    (' O1D  -->  OH + OH ',            'k[O1D-->OH+OH]=CONST(2.14D-10*H2O)'),
    (' OH + O3  -->  HO2 ',            'k[OH+O3-->HO2]=CONST(1.70D-12*EXP(-940/(T)))'),
    (' OH + H2  -->  HO2 ',            'k[OH+H2-->HO2]=CONST(7.70D-12*EXP(-2100/(T)))'),
    (' OH + CO  -->  HO2 + CO2 ',      'k[OH+CO-->HO2+CO2]{APPEND_CO2}=CONST(KMT05)'),
    (' OH + H2O2  -->  HO2 ',          'k[OH+H2O2-->HO2]=CONST(2.90D-12*EXP(-160/(T)))'),
    (' HO2 + O3  -->  OH ',            'k[HO2+O3-->OH]=CONST(2.03D-16*(((T)/300)@4.57*EXP(693/(T))))'),
    (' OH + HO2  -->   ',              'k[OH+HO2-->]=CONST(4.80D-11*EXP(250/(T)))'),
    (' HO2 + HO2  -->  H2O2 ',         'k[HO2+HO2-->H2O2]=CONST(2.20D-13*KMT06*EXP(600/(T)) + 1.90D-33*M*KMT06*EXP(980/(T)))'),
    (' OH + NO  -->  HONO ',           'k[OH+NO-->HONO]=CONST(KMT07)'),
    (' OH + NO2  -->  HNO3 ',          'k[OH+NO2-->HNO3]=CONST(KMT08)'),
    (' OH + NO3  -->  HO2 + NO2 ',     'k[OH+NO3-->HO2+NO2]=CONST(2.00D-11)'),
    (' HO2 + NO  -->  OH + NO2 ',      'k[HO2+NO-->OH+NO2]=CONST(3.45D-12*EXP(270/(T)))'),
    (' HO2 + NO2  -->  HO2NO2 ',       'k[HO2+NO2-->HO2NO2]=CONST(KMT09)'),
    (' HO2NO2  -->  HO2 + NO2 ',       'k[HO2NO2-->HO2+NO2]=CONST(KMT10)'),
    (' OH + HO2NO2  -->  NO2 ',        'k[OH+HO2NO2-->NO2]=CONST(3.20D-13*EXP(690/(T)))'),
    (' HO2 + NO3  -->  OH + NO2 ',     'k[HO2+NO3-->OH+NO2]=CONST(4.00D-12)'),
    (' OH + HONO  -->  NO2 ',          'k[OH+HONO-->NO2]=CONST(2.50D-12*EXP(260/(T)))'),
    (' OH + HNO3  -->  NO3 ',          'k[OH+HNO3-->NO3]=CONST(KMT11)'),
    (' O + SO2  -->  SO3 ',            'k[O+SO2-->SO3]=CONST(4.00D-32*EXP(-1000/(T))*M)'),
    (' OH + SO2  -->  HSO3 ',          'k[OH+SO2-->HSO3]=CONST(KMT12)'),
    (' HSO3  -->  HO2 + SO3 ',         'k[HSO3-->HO2+SO3]=CONST(1.30D-12*EXP(-330/(T))*O2)'),
    # Newly added since MCMv3.2 (sink reactions); not in the IDL constants
    # table but always present in the canonical anorg.mcm331.dark.mech file
    # shipped with EASY/JUICE.
    (' HNO3  -->  NA ',                'k[HNO3-->NA]=CONST(0.0)'),
    (' N2O5  -->  NA + NA ',           'k[N2O5-->NA+NA]=CONST(0.0)'),
    (' SO3  -->  SA ',                 'k[SO3-->SA]=CONST(0.0)'),
]

# OH+HCHO / NO3+HCHO are part of the MCM VOC mechanism (when HCHO appears),
# NOT part of the inorganic core. The IDL converter only used these when in
# /INITIAL graph-walk mode for filtering. We don't emit them in anorg files.
INORGANIC_SOURCE_DARK = []

# Inorganic photolysis (only emitted when not /DARK)
INORGANIC_PHOT = [
    (' O3 + hv  -->  O1D ',            'k[O3+hv-->O1D]{J<1>}=PHTFKT(6.073D-05,1.743,-0.474)'),
    (' O3 + hv  -->  O ',              'k[O3+hv-->O]{J<2>}=PHTFKT(4.775D-04,0.298,-0.080)'),
    (' H2O2 + hv  -->  OH + OH ',      'k[H2O2+hv-->OH+OH]{J<3>}=PHTFKT(1.041D-05,0.723,-0.279)'),
    (' NO2 + hv  -->  NO + O ',        'k[NO2+hv-->NO+O]{J<4>}=PHTFKT(1.165D-02,0.244,-0.267)'),
    (' NO3 + hv  -->  NO ',            'k[NO3+hv-->NO]{J<5>}=PHTFKT(2.485D-02,0.168,-0.108)'),
    (' NO3 + hv  -->  NO2 + O ',       'k[NO3+hv-->NO2+O]{J<6>}=PHTFKT(1.747D-01,0.155,-0.125)'),
    (' HONO + hv  -->  OH + NO ',      'k[HONO+hv-->OH+NO]{J<7>}=PHTFKT(2.644D-03,0.261,-0.288)'),
    (' HNO3 + hv  -->  OH + NO2 ',     'k[HNO3+hv-->OH+NO2]{J<8>}=PHTFKT(9.312D-07,1.230,-0.307)'),
]

# HCHO photolysis is part of the VOC mechanism (the MCM .fac file emits J<11>
# / J<12> rates for HCHO directly), not the anorg core. We don't emit it.
INORGANIC_SOURCE_PHOT = []

# Inorganic species that are always added to the dilution list when /DIL,
# even if they don't appear in the VOC mechanism's VARIABLE block.
SAPHIR_INORGANIC_DIL_SPECIES = [
    'H2O', 'CO', 'CO2',
    'H2O2', 'HO2', 'HONO', 'HNO3', 'HO2NO2', 'HSO3',
    'N2O5', 'NO', 'NO2', 'NO3',
    'O', 'O1D', 'O3', 'OH',
    'SO2', 'SO3',
]


# ─────────────────────────────────────────────────────────────────────────────
# Parsing
# ─────────────────────────────────────────────────────────────────────────────

def _read_lines(path: Path) -> list[str]:
    """Read MCM file using latin-1 (legacy umlaut comments)."""
    with open(path, 'r', encoding='latin-1') as f:
        return f.read().split('\n')


def _combine_continuations(lines: list[str]) -> list[str]:
    """Combine reaction continuation lines.

    A `% rate : reaction` line that lacks a closing `;` continues onto the
    next line. The IDL converter merges them into one logical line.
    """
    out = []
    i = 0
    while i < len(lines):
        line = lines[i]
        # Match IDL: if `%` is present but no `;` on this line, append next
        if '%' in line and ';' not in line and i + 1 < len(lines):
            line = line + lines[i + 1]
            i += 2
        else:
            i += 1
        out.append(line)
    return out


def _extract_variable_block(lines: list[str]) -> list[str]:
    """Pull species names from the VARIABLE ... ; block.

    The block begins with a line starting with 'VARIABLE' and continues until
    a line ending with ';'.
    """
    species = []
    in_block = False
    for line in lines:
        s = line.strip()
        if not in_block:
            if s.startswith('VARIABLE'):
                in_block = True
            continue
        # In block — accumulate tokens until ';'
        cleaned = s.rstrip(';').strip()
        toks = cleaned.split()
        species.extend(toks)
        if s.endswith(';'):
            break
    return [t for t in species if t]


def _extract_ro2_block(lines: list[str]) -> list[str]:
    """Extract the RO2 species sum from the `RO2 = sp1 + sp2 + ... ;` block.

    Multi-line; terminated by `;`. Returns ordered list of species names.
    """
    text = ''
    in_block = False
    for line in lines:
        s = line.strip()
        if not in_block:
            if s.startswith('RO2 ') and '=' in s:
                in_block = True
                # Take everything after the '='
                text = s.split('=', 1)[1]
            continue
        text += ' ' + s
        if ';' in text:
            text = text.split(';', 1)[0]
            break
    parts = [p.strip() for p in text.split('+') if p.strip()]
    return parts


def _extract_reactions(lines: list[str]) -> list[tuple[str, str]]:
    """Pull (rate_expr, reaction_str) pairs from `% rate : reaction ;` lines.

    Skips `% (...)`, `% RINJECT`, `% LT...` administrative lines per IDL.
    """
    rxns = []
    for raw in lines:
        line = raw.strip()
        if not line.startswith('%'):
            continue
        # IDL guards: skip "% (", "% RINJECT", "% LT"
        if line.startswith('% (') or line.startswith('% RINJECT') or line.startswith('% LT'):
            continue
        # Strip leading '%' and trailing ';'
        body = line[1:].rstrip(';').strip()
        if ':' not in body:
            continue
        rate_part, mech_part = body.split(':', 1)
        rxns.append((rate_part.strip(), mech_part.strip()))
    return rxns


# ─────────────────────────────────────────────────────────────────────────────
# Reaction conversion
# ─────────────────────────────────────────────────────────────────────────────

def _normalise_rate_expr(expr: str) -> str:
    """Apply the IDL rate-expression rewrites:

    1. Replace `J<n>` placeholder syntax with `Jn` (e.g. `J<13>` → `J13`)
    2. Replace `TEMP` with `(T)` (the EASY parameter name)
    3. Fix the FACSIMILE @-exponent parenthesis bug for T@2..T@5
       (`/TEMP)@2)` should be `/TEMP@2)` because @ binds tighter than `)`)
    4. Strip whitespace
    """
    e = expr
    # Strip J<n>
    e = re.sub(r'J<(\d+)>', r'J\1', e)
    # First handle the @-exponent parenthesis fix BEFORE substituting TEMP
    # IDL does this in the opposite order — but its substitution targets
    # `/TEMP)@N)` literally, so we replicate the same.
    e = e.replace('/TEMP)@2)', '/TEMP@2)')
    e = e.replace('/TEMP)@3)', '/TEMP@3)')
    e = e.replace('/TEMP)@4)', '/TEMP@4)')
    e = e.replace('/TEMP)@5)', '/TEMP@5)')
    # Now substitute TEMP -> (T)
    e = e.replace('TEMP', '(T)')
    # Remove all whitespace (matches IDL STRCOMPRESS /REMOVE_ALL)
    e = re.sub(r'\s+', '', e)
    return e


def _convert_reaction(rate_raw: str, mech_raw: str) -> list[tuple[str, str]]:
    """Convert a single MCM reaction to one or more (mech, rate_def) pairs.

    Returns a list because equilibrium reactions of the form `K1 % K2`
    expand into two reactions (forward and reverse).
    """
    # Mechanism: replace `=` with `-->`, normalise spacing
    if '=' not in mech_raw:
        return []
    lhs, rhs = mech_raw.split('=', 1)
    lhs = lhs.strip()
    rhs = rhs.strip()

    rate_expr = _normalise_rate_expr(rate_raw)

    # Handle equilibrium-rate syntax: `K1 % K2` (forward % reverse)
    if '%' in rate_expr:
        k1, k2 = rate_expr.split('%', 1)
        # Forward reaction with k1
        # If k1 contains 'J' (uppercase J followed by digit), it's photolysis
        # — but equilibrium reactions are dark/thermal, so don't add hv
        out = []
        fwd_mech = f'{lhs}-->{rhs}'
        out.append((fwd_mech, f'k[{fwd_mech}]=CONST({k1})'))
        rev_mech = f'{rhs}-->{lhs}'
        out.append((rev_mech, f'k[{rev_mech}]=CONST({k2})'))
        return [_finalise(m, r) for m, r in out]

    # Photolysis detection: rate references J1..J61 (after J<n> normalisation)
    # IDL trigger: STRPOS(expr, 'CONST(J') NE -1
    is_photolysis = bool(re.search(r'\bJ\d+\b', rate_expr))

    if is_photolysis:
        # Insert `+ hv` into the LHS so the reaction reads `A + hv -->`
        lhs_with_hv = f'{lhs}+hv'
    else:
        lhs_with_hv = lhs

    mech = f'{lhs_with_hv}-->{rhs}'
    rate_def = f'k[{mech}]=CONST({rate_expr})'
    return [_finalise(mech, rate_def)]


def _finalise(mech: str, rate_def: str) -> tuple[str, str]:
    """Final whitespace normalisation: collapse all whitespace inside the
    rate-def k[...] block but produce a readable mech line."""
    # Mech line: ' A + B --> C + D ' (spaces around tokens, matches IDL output)
    # Build by splitting on `-->`
    if '-->' in mech:
        m_lhs, m_rhs = mech.split('-->', 1)
    else:
        m_lhs, m_rhs = mech, ''
    m_lhs = ' + '.join(t.strip() for t in m_lhs.split('+'))
    m_rhs = ' + '.join(t.strip() for t in m_rhs.split('+'))
    mech_pretty = f' {m_lhs}  -->  {m_rhs} '
    # Rate def: collapse whitespace
    rate_def = re.sub(r'\s+', '', rate_def)
    return mech_pretty, rate_def


# ─────────────────────────────────────────────────────────────────────────────
# Used-constant detection
# ─────────────────────────────────────────────────────────────────────────────

def _find_used_constants(rate_defs: list[str]) -> set[str]:
    """Scan rate definitions for references to KMT/KFPAN/Jn etc.

    Uses word-boundary regex to avoid matching J1 inside J11.
    """
    used = set()
    joined = '\n'.join(rate_defs)
    for name in KNOWN_NAMES_ORDER:
        # \b...\b for word boundaries (J1 won't match inside J11)
        if re.search(r'\b' + re.escape(name) + r'\b', joined):
            used.add(name)
    return used


# ─────────────────────────────────────────────────────────────────────────────
# Output writers
# ─────────────────────────────────────────────────────────────────────────────

def _write_mech_file(path: Path, reactions: list[tuple[str, str]]) -> None:
    with open(path, 'w', encoding='latin-1', newline='\n') as f:
        for mech, _ in reactions:
            f.write(mech + '\n')


def _write_rate_file(path: Path, reactions: list[tuple[str, str]]) -> None:
    with open(path, 'w', encoding='latin-1', newline='\n') as f:
        for _, rate in reactions:
            f.write(rate + '\n')


def _write_ro2_file(path: Path, ro2_species: list[str]) -> None:
    """Emit `RO2 = CONST( 0.0 +  sp1 + sp2 + ...)` on a single line."""
    body = ' + '.join(ro2_species) if ro2_species else ''
    line = f'RO2 = CONST( 0.0 +  {body})' if body else 'RO2 = CONST( 0.0)'
    with open(path, 'w', encoding='latin-1', newline='\n') as f:
        f.write(line + '\n')


def _write_incl_file(path: Path, used: set[str], input_filename: str) -> None:
    """Emit the .incl file: header, base functions, aliases, photolysis,
    then the subset of KMT/J definitions actually referenced."""
    from datetime import datetime
    lines = [
        '  ',
        '; This is a include file for a EASY-Input File written by mcm_converter.py',
        f'; with input file {input_filename} at {datetime.now().strftime("%a %b %d %H:%M:%S %Y")}',
        '; Comments and references can be found in mcmex331.incl',
        '',
        '; Funktionen:',
        'CONST=FUNCTION[%1]',
        'DENSITY=FUNCTION[(%1)/((%2)*1.379E-19)]',
        'T =CONST(273.2)',
        'P =CONST(1013.)',
        'M =DENSITY(P,T)',
        'O2=CONST(0.20946*M)',
        'N2=CONST(0.78084*M)',
        '',
        '; Aliasse:',
        'ALIAS[HV]     = hv',
        'ALIAS[OH]     = HO',
        'ALIAS[HCHO]   = H2CO CH2O',
        'ALIAS[O]      = O3P',
        'ALIAS[CH3OOH] = CH3O2H MOOH',
        'ALIAS[HONO]   = HNO2',
        'ALIAS[HNO3]   = HONO2',
        'ALIAS[HNO4]   = HO2NO2',
        '',
        ';Photolysis:',
        'HV=CONST(1.)',
        'SZA=CONST(45.) ; Sonnen Zenit Winkel',
        'PHTFKT  = FUNCTION[( (%1)*((COS(AMIN(SZA,89.9)*0.0174533))@(%2))*EXP((%3)*1./COS(AMIN(SZA,89.9)*0.0174533)))]',
        '',
    ]
    # Emit constants in canonical order
    for name in KNOWN_NAMES_ORDER:
        if name in used:
            lines.append(KNOWN_CONSTANTS[name])
    with open(path, 'w', encoding='latin-1', newline='\n') as f:
        f.write('\n'.join(lines) + '\n')


# ─────────────────────────────────────────────────────────────────────────────
# Precursor-subset graph walk (IDL /INITIAL= equivalent)
# ─────────────────────────────────────────────────────────────────────────────

def _strip_coeff(token: str) -> str:
    """Strip a leading stoichiometric coefficient: '2*HO2' -> 'HO2'."""
    t = token.strip()
    if '*' in t:
        return t.split('*')[-1].strip()
    return t


def _parse_mech_sides(mech_line: str) -> tuple[list[str], list[str]]:
    """Split a mech line ' A + B  -->  C + D ' into ([A, B], [C, D]).

    Strips stoichiometric prefixes from each species name.
    """
    if '-->' not in mech_line:
        return [], []
    lhs_str, rhs_str = mech_line.split('-->', 1)
    lhs = [_strip_coeff(t) for t in lhs_str.split('+') if t.strip()]
    rhs = [_strip_coeff(t) for t in rhs_str.split('+') if t.strip()]
    return [t for t in lhs if t], [t for t in rhs if t]


def _graph_walk_filter(
    reactions: list[tuple[str, str]],
    precursors: list[str],
    *,
    source: bool = True,
    dark: bool = False,
    dil: bool = True,
) -> tuple[list[tuple[str, str]], set[str]]:
    """Filter reactions to only those reachable from the precursor set.

    Mirrors the IDL `/INITIAL=` graph-walk:
      1. Build initial valid_species = precursors + always-valid extras
         (H2O if /DIL, HONO+HCHO if /SOURCE, hv if not /DARK)
      2. Iterate: any reaction whose reactants are ALL valid becomes valid;
         all its products are added to valid_species. Repeat to fixed point.
      3. Return the subset of reactions whose all-reactants-valid condition
         was satisfied, plus the final valid_species set.

    Args:
        reactions: list of (mech_line, rate_def) tuples — pass the full
                   merged set of VOC reactions + inorganic core reactions
        precursors: species names the user wants to study (e.g. ['C2H5OH'])

    Returns:
        (filtered_reactions, valid_species_set)
    """
    valid = {s.upper() for s in precursors}
    # IDL /INITIAL implicit additions:
    if dil:
        valid.add('H2O')
    if source:
        valid.add('HONO')
        valid.add('HCHO')
    if not dark:
        valid.add('HV')
        valid.add('hv')  # IDL preserves the lowercase 'hv' literal in the mech
    # M / O2 / N2 are bulk gases — always valid (they appear in anorg reactants)
    valid.update({'M', 'O2', 'N2'})

    # Pre-parse all reactions once
    parsed = []  # (mech, rate, reactants_set, products_set)
    for mech, rate in reactions:
        lhs, rhs = _parse_mech_sides(mech)
        parsed.append((mech, rate, lhs, rhs))

    valid_idx = set()
    while True:
        added = False
        for i, (_, _, lhs, _) in enumerate(parsed):
            if i in valid_idx:
                continue
            # Empty LHS → can't decide; treat as invalid (zeroth-order sources
            # don't appear in MCM standard reactions)
            if not lhs:
                continue
            if all(r.upper() in {v.upper() for v in valid} for r in lhs):
                valid_idx.add(i)
                added = True
                # Add products
                for p in parsed[i][3]:
                    if p.upper() not in {v.upper() for v in valid}:
                        valid.add(p)
                        added = True
        if not added:
            break

    filtered = [(parsed[i][0], parsed[i][1]) for i in sorted(valid_idx)]
    # Drop the special hv tokens from valid set before returning (IDL does this)
    valid_clean = {s for s in valid if s.upper() != 'HV'}
    return filtered, valid_clean


# ─────────────────────────────────────────────────────────────────────────────
# SAPHIR-mode output writers
# ─────────────────────────────────────────────────────────────────────────────

def _write_anorg_files(output_dir: Path, dark: bool = False, source: bool = True) -> None:
    """Emit the hardcoded anorg.mcm331.{dark,phot}.{mech,rate} files.

    These are stable inorganic chemistry constants that don't depend on the
    input MCM file — exactly the same in every conversion.
    """
    dark_rxns = list(INORGANIC_DARK)
    if source:
        dark_rxns.extend(INORGANIC_SOURCE_DARK)

    with open(output_dir / 'anorg.mcm331.dark.mech', 'w', encoding='latin-1', newline='\n') as f:
        for mech, _ in dark_rxns:
            f.write(mech + '\n')
    with open(output_dir / 'anorg.mcm331.dark.rate', 'w', encoding='latin-1', newline='\n') as f:
        for _, rate in dark_rxns:
            f.write(rate + '\n')

    if not dark:
        phot_rxns = list(INORGANIC_PHOT)
        if source:
            phot_rxns.extend(INORGANIC_SOURCE_PHOT)
        with open(output_dir / 'anorg.mcm331.phot.mech', 'w', encoding='latin-1', newline='\n') as f:
            for mech, _ in phot_rxns:
                f.write(mech + '\n')
        with open(output_dir / 'anorg.mcm331.phot.rate', 'w', encoding='latin-1', newline='\n') as f:
            for _, rate in phot_rxns:
                f.write(rate + '\n')


def _write_dil_files(output_dir: Path, stem: str, variables: list[str], source: bool) -> int:
    """Emit `<stem>.dil.mech` and `<stem>.dil.rate` for every species.

    Merges the VOC VARIABLE block with the SAPHIR_INORGANIC_DIL_SPECIES list,
    dedupes, and sorts alphabetically (matches IDL output order).
    """
    all_species = set(SAPHIR_INORGANIC_DIL_SPECIES) | set(variables)
    if source:
        all_species.add('HCHO')
    sorted_sp = sorted(s for s in all_species if s)

    dil_mech = output_dir / f'{stem}.dil.mech'
    dil_rate = output_dir / f'{stem}.dil.rate'
    with open(dil_mech, 'w', encoding='latin-1', newline='\n') as fm, \
         open(dil_rate, 'w', encoding='latin-1', newline='\n') as fr:
        for sp in sorted_sp:
            fm.write(f'{sp} + DIL --> \n')
            fr.write(f'k[{sp} + DIL --> ] = CONST(DILUTE)\n')
    return len(sorted_sp)


def _write_easy_temp(output_dir: Path, stem: str, input_filename: str,
                     dark: bool, source: bool, walldep: bool, bgreac: bool,
                     dil: bool, measured: bool, include_anorg: bool = True) -> Path:
    """Emit the master `<stem>.saphir.easy_temp` EASY input file.

    This is the entry point that JUICE / FACSIMILE consumes — it #includes
    the .mech, .rate, .ro2.incl, .incl, .dil.mech, .dil.rate files plus
    boilerplate for the SAPHIR-specific chemistry blocks.

    When include_anorg=False (precursor/init mode), the standalone anorg
    files are not #include'd because the inorganic core has been inlined
    into the per-mech .init.mech/.init.rate files (and graph-walk filtered).
    """
    from datetime import datetime
    timestamp = datetime.now().strftime('%a %b %d %H:%M:%S %Y')

    lines = [
        '  ',
        '; This is a temporary EASY-Input File written by mcm_converter.py',
        f'; with input file {input_filename} at {timestamp}',
        ';   ',
        '; To make this file a working EASY-Input File, the following',
        '; steps have to be done at least: ',
        ';   1) Fill in the INITIAL values (and water concentration)',
        ';   2) Make sure that EASY will find the included files',
        ';   3) Save the new file with a DOS compatible name',
        '   ',
        '; INCLUDES: ',
    ]
    if include_anorg:
        lines += [
            '#include anorg.mcm331.dark.rate   ; dark anorganic reactions',
            '#include anorg.mcm331.dark.mech   ; rate constants of dark reactions',
        ]
        if not dark:
            lines += [
                '#include anorg.mcm331.phot.rate   ; anorganic photolysis reactions',
                '#include anorg.mcm331.phot.mech   ; anorganic photolysis frequencies',
            ]
    lines += [
        f'#include {stem}.incl ; MCM3-spezific constants and functions ',
        f'#include {stem}.ro2.incl ; RO2 species',
        f'#include {stem}.mech ; reaction mechanism for given VOC',
        f'#include {stem}.rate ; rate constants for given VOC reactions',
        ' ',
        '*********************************************************************',
        ' ',
    ]

    if walldep:
        lines += [
            ';/WALLDEP, Wall deposition in SAPHIR chamber',
            '#include saphir.walldep.incl',
            'WALLDEP   = CONST(1.) ; Switch fuer die Wandverluste ',
            ' ',
            '*********************************************************************',
            ' ',
        ]

    if source:
        lines += [
            ';/SOURCE, HONO and HCHO source of SAPHIR chamber',
            '',
            '; Funktion die die Wasserkonzentration in die relative',
            '; Feuchtigkeit umrechnet',
            'WVC2RH     =FUNCTION[%1*100./(exp(21.36469-5339.66/T)/P*M)]',
            '',
            '; Relative Feuchtigkeit',
            'RH = WVC2RH(H2O)',
            '',
            '; HONO Quelle aus Paper Rohrer et al.',
            'hv + SHONO --> HONO',
            'k[hv + SHONO --> HONO]  = CONST( 4.7e13 * (1. + (RH/11.6)@2) * exp(-3950./T) * jno2 )',
            'SHONO = CONST(1.)',
            '',
            '; HCHO Quelle (licht abhaengig, photochemie unabhaengig)',
            'hv + SHCHO --> HCHO',
            'k[hv + SHCHO --> HCHO] = CONST(2.8e6)',
            'SHCHO = CONST(1.)',
        ]
        if not measured:
            lines += [
                '',
                'jno2 = PHTFKT(1.165D-02,0.244,-0.267)',
            ]
        lines += [
            ' ',
            '*********************************************************************',
            ' ',
        ]

    if bgreac:
        lines += [
            ';/BGREAC, Background reactivity (see saphir.bgreac.incl) ',
            'Y + OH --> HO2',
            'k[Y + OH --> HO2]=CONST(2.4e-13)',
            'Y=CONST(200.e-9*M) ; Y in CO aequivalenten',
            ' ',
            '*********************************************************************',
            ' ',
        ]

    if dil:
        lines += [
            '; DILUTION, /DIL',
            f'#include {stem}.dil.mech ;dilution mechanism',
            f'#include {stem}.dil.rate ;dilution rate',
            'FLOW = CONST(6.)  ; Experimetal flow in m^3/h',
            'DILUTE = CONST(FLOW/(270.*3600.))',
            'DIL = CONST(1.) ; Switch for dilution',
            ' ',
            '; Water as a variable ',
            'INITIAL[H2O] = CONST(1.e-4*M)',
            ' ',
            '*********************************************************************',
            ' ',
        ]
    else:
        lines += [
            '; Water as a constant / given parameter',
            'H2O = CONST(1.e-4*M)',
            ' ',
            '*********************************************************************',
            ' ',
        ]

    if measured:
        lines += [
            ';/MEASURED, use measured photolysis frequencies (given in ENZ-file)',
            ';           and other measured parameters]',
        ]
        if not dark:
            lines += [
                ';#include saphir.jmeasured_02.incl ; definitions for further use of jno2',
                '                              ; and jo1d as scaling factors of other ',
                '                              ; photolysis frequencies',
                'jno2 = INPUT(jno2)',
                'k[NO2 + HV --> O3 + NO] = CONST( jno2 )',
                'jo1d = INPUT(jo1d)',
                'k[O3 + HV --> O1D ] = CONST( jo1d )',
                'k[HCHO + HV --> H2 + CO ] = INPUT(jhcho_m)',
                'k[HCHO + HV --> 2*HO2 + CO ] = INPUT(jhcho_r)',
                'k[H2O2 + HV --> 2*HO] = INPUT(jh2o2)',
                'k[HONO + HV --> OH + NO] = INPUT(jhono)',
                ' ',
            ]
        lines += [
            'T    = INPUT(t)   ; temperature in K',
            'P    = INPUT(p)   ; pressure in hPa',
            'SZA  = INPUT(sz)  ; solar zenith angle',
            'HV   = INPUT(hv)  ; switch for light',
            'FLOW = INPUT(dilute)  ; flow into the chamber in m^3/h',
            'H2O  = INPUT(H2O) ; mixing ratio of water',
            ' ',
            '*********************************************************************',
            ' ',
        ]

    lines += [
        '; FILES: ',
        'FILES[ENZ]   = mcm3.enz',
        ' ',
        '; FACSIMILE PARAMETERS',
        'FACS[HMAX]   = 30',
        'FACS[OUTSTEP]= 100000',
        ' ',
        '; INITIALS:',
        'INITIAL[*] = CONST(0.) ',
    ]

    path = output_dir / f'{stem}.saphir.easy_temp'
    with open(path, 'w', encoding='latin-1', newline='\n') as f:
        f.write('\n'.join(lines) + '\n')
    return path


# ─────────────────────────────────────────────────────────────────────────────
# Lightweight species peek (for GUI validation before full conversion)
# ─────────────────────────────────────────────────────────────────────────────

def list_species(input_path: str | Path) -> list[str]:
    """Quickly extract the species list from a raw MCM .fac file.

    Parses only the VARIABLE block — much faster than a full conversion.
    Used by the GUI to validate user-supplied precursor names before
    running convert_mcm_file().
    """
    raw_lines = _read_lines(Path(input_path))
    return _extract_variable_block(raw_lines)


# ─────────────────────────────────────────────────────────────────────────────
# Main entry point
# ─────────────────────────────────────────────────────────────────────────────

def convert_mcm_file(
    input_path: str | Path,
    output_dir: str | Path | None = None,
    *,
    saphir: bool = False,
    dark: bool = False,
    measured: bool = False,
    walldep: bool | None = None,
    source: bool | None = None,
    bgreac: bool | None = None,
    dil: bool | None = None,
    precursors: list[str] | None = None,
) -> dict:
    """Convert a raw MCM .fac file into EASY include files.

    Args:
        input_path: path to raw MCM .fac file (downloaded from York website)
        output_dir: directory to write output files (default: input dir)
        saphir: master switch — enables walldep, source, bgreac, dil and emits
                the .saphir.easy_temp master file plus dilution files
        dark: skip photolysis (no .phot files, no source)
        measured: use measured J-values from ENZ file (jno2/jo1d/jhcho/...)
        walldep, source, bgreac, dil: individual SAPHIR sub-switches; default
                None means "follow saphir flag". Set to True/False to override.
        precursors: list of precursor species names (e.g. ['C2H5OH']). If
                given, applies the IDL /INITIAL= graph walk: only reactions
                reachable from the precursor set are kept, the inorganic core
                is inlined into the per-mech file (no separate anorg files),
                output filenames get a '.init.' insert, and the RO2 / dilution
                lists are filtered to only the reachable species.

    Returns a dict with the output file paths and counts.
    """
    input_path = Path(input_path)
    if output_dir is None:
        output_dir = input_path.parent
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    stem = input_path.stem
    init_mode = bool(precursors)
    init_tag = '.init' if init_mode else ''

    # Resolve SAPHIR sub-switches
    if saphir:
        if walldep is None: walldep = True
        if source  is None: source  = not dark
        if bgreac  is None: bgreac  = True
        if dil     is None: dil     = True
    else:
        if walldep is None: walldep = False
        if source  is None: source  = False
        if bgreac  is None: bgreac  = False
        if dil     is None: dil     = False

    # Parse the MCM input
    raw_lines = _read_lines(input_path)
    lines = _combine_continuations(raw_lines)
    variables = _extract_variable_block(lines)
    ro2_species = _extract_ro2_block(lines)
    raw_rxns = _extract_reactions(lines)

    # Convert VOC reactions
    voc_reactions = []
    for rate_raw, mech_raw in raw_rxns:
        voc_reactions.extend(_convert_reaction(rate_raw, mech_raw))

    # ── Precursor mode: graph-walk filter ───────────────────────────────────
    if init_mode:
        # Merge VOC reactions with the inorganic core so the walk can
        # reach inorganic chemistry through degradation products.
        inorganic = list(INORGANIC_DARK)
        if not dark:
            inorganic += list(INORGANIC_PHOT)
        merged = voc_reactions + inorganic

        filtered, valid_species = _graph_walk_filter(
            merged, precursors, source=source, dark=dark, dil=dil)
        if not filtered:
            raise RuntimeError(
                f"Graph walk found NO reactions reachable from precursors "
                f"{precursors}. Check spelling and that the species exists "
                f"in {input_path.name}.")
        reactions = filtered
        # Restrict variables list and RO2 pool to the reachable subset
        variables = sorted(s for s in valid_species if s)
        ro2_species = [r for r in ro2_species if r.upper() in valid_species]
    else:
        reactions = voc_reactions

    # For .incl constant detection
    rate_defs_for_scan = [r for _, r in reactions]
    if saphir and not init_mode:
        # In non-init mode the inorganic core lives in separate anorg files
        # but we still scan its rates so e.g. KMT01 referenced by the anorg
        # block gets emitted in our .incl too.
        rate_defs_for_scan += [r for _, r in INORGANIC_DARK]
        if not dark:
            rate_defs_for_scan += [r for _, r in INORGANIC_PHOT]
    used = _find_used_constants(rate_defs_for_scan)

    # ── Write VOC chemistry files ────────────────────────────────────────────
    mech_path = output_dir / f'{stem}{init_tag}.mech'
    rate_path = output_dir / f'{stem}{init_tag}.rate'
    ro2_path  = output_dir / f'{stem}{init_tag}.ro2.incl'
    incl_path = output_dir / f'{stem}{init_tag}.incl'

    _write_mech_file(mech_path, reactions)
    _write_rate_file(rate_path, reactions)
    _write_ro2_file(ro2_path, ro2_species)
    _write_incl_file(incl_path, used, input_path.name)

    result = {
        'mech': mech_path,
        'rate': rate_path,
        'ro2':  ro2_path,
        'incl': incl_path,
        'n_species':   len(variables),
        'n_reactions': len(reactions),
        'n_ro2':       len(ro2_species),
        'n_constants': len(used),
    }
    if init_mode:
        result['precursors'] = list(precursors)

    # ── SAPHIR-mode extras ───────────────────────────────────────────────────
    if saphir:
        # In init mode we do NOT emit the standalone anorg files — the
        # inorganic core is already inlined into the .init.mech/.rate file
        # and the easy_temp does not #include them.
        if not init_mode:
            _write_anorg_files(output_dir, dark=dark, source=source)
            result['anorg_dark_mech'] = output_dir / 'anorg.mcm331.dark.mech'
            result['anorg_dark_rate'] = output_dir / 'anorg.mcm331.dark.rate'
            if not dark:
                result['anorg_phot_mech'] = output_dir / 'anorg.mcm331.phot.mech'
                result['anorg_phot_rate'] = output_dir / 'anorg.mcm331.phot.rate'

        if dil:
            dil_stem = f'{stem}{init_tag}'
            n_dil = _write_dil_files(output_dir, dil_stem, variables, source=source)
            result['dil_mech'] = output_dir / f'{dil_stem}.dil.mech'
            result['dil_rate'] = output_dir / f'{dil_stem}.dil.rate'
            result['n_dil_species'] = n_dil

        easy_temp = _write_easy_temp(
            output_dir, stem + init_tag, input_path.name,
            dark=dark, source=source, walldep=walldep,
            bgreac=bgreac, dil=dil, measured=measured,
            include_anorg=not init_mode)
        result['easy_temp'] = easy_temp

    return result


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description="Convert raw MCM .fac file to EASY include files")
    ap.add_argument('input', help='Input MCM .fac file')
    ap.add_argument('output_dir', nargs='?', default=None, help='Output directory')
    ap.add_argument('--saphir', action='store_true', help='Enable SAPHIR mode (walldep, source, bgreac, dil, easy_temp)')
    ap.add_argument('--dark', action='store_true', help='Disable photolysis')
    ap.add_argument('--measured', action='store_true', help='Use measured J values from ENZ file')
    ap.add_argument('--precursors', default=None,
                    help='Comma-separated list of precursor species (e.g. C2H5OH or "C2H5OH,CH4"). '
                         'Triggers IDL /INITIAL= graph walk: only reactions reachable from the '
                         'precursor set are kept; output filenames get a .init insert; the '
                         'inorganic core is inlined into the per-mech file.')
    args = ap.parse_args()

    precursors = None
    if args.precursors:
        precursors = [s.strip() for s in args.precursors.split(',') if s.strip()]

    result = convert_mcm_file(
        args.input, args.output_dir,
        saphir=args.saphir, dark=args.dark, measured=args.measured,
        precursors=precursors)

    print(f"[mcm_converter] Wrote:")
    for k in ('mech', 'rate', 'ro2', 'incl'):
        print(f"  {result[k]}")
    if 'easy_temp' in result:
        print(f"  {result['easy_temp']}")
        if 'dil_mech' in result:
            print(f"  {result['dil_mech']}  ({result['n_dil_species']} species)")
            print(f"  {result['dil_rate']}")
        for k in ('anorg_dark_mech', 'anorg_dark_rate', 'anorg_phot_mech', 'anorg_phot_rate'):
            if k in result:
                print(f"  {result[k]}")
    print(f"  species: {result['n_species']}  reactions: {result['n_reactions']}"
          f"  RO2 pool: {result['n_ro2']}  constants: {result['n_constants']}")
