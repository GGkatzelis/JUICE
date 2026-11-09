#!/usr/bin/env python
"""Generate SOA_FZJ_2026.include -- reversible gas<->particle<->wall partitioning for every
condensable SIM-HOM/LIMO/BPIN species (SVOC+LVOC+ELVOC, closed-shell, C*<300 ug/m3).

Kinetics = Krechmer, Day & Jimenez 2020 (ES&T 54:12890). ALL species get the SAME 6-reaction
reversible set; the per-species C* alone sets reversibility (low C* -> ~irreversible automatically),
so no SVOC/ELVOC boundary has to be chosen. The VBS group label only routes each species into the
right diagnostic element sum.

Per species X (6 reactions):
  1  X       --> SOA_X                              CONST(KCOND)            condensation onto particles
  2  SOA_X   --> X                                  CONST(KCOND*C*/COA)     particle evaporation
  3  X       --> WALL_X                             CONST(KFWALL)           gas->wall deposition
  4  WALL_X  --> X                                  CONST(KFWALL*C*/Cw_X)   wall evaporation (Cw_X = C*-dep.)
  5  SOA_X   --> WALL_PL_X                          CONST(KPWALL)           particle wall loss (irrev., PER-SPECIES)
  6  SOA_X + DIL -->                                CONST(DILUTE)           dilution of suspended SOA
Plus:  COA (live organic-aerosol mass, ug/m3) and 9 SOA_ + 9 WALL_ + 9 WALL_PL_ element sums (C/O/N by group).

Importable: read_condensable(param) -> species list; emit_soa(sp, out_path) -> writes include.
Used by make_subset.py to regenerate a reduced include into a subset folder (never the canonical path).
"""
import re
ROOT=r"E:/Python Code/JUICE archive/EASY/include files/SIM_HOM"
PARAM=ROOT+"/SIMHOM_volatility.param"
OUT=ROOT+"/SOA_FZJ_2026.include"
# C*-dependent wall capacity Cw -- Krechmer et al. 2016 (ES&T 50,5757, Table 1):
#   Cw = 16*(C*)^0.6 ug/m3 for C*>1 ; Cw = 16 for C*<1 (plateau). Cw is per-chamber-VOLUME so it
#   scales with A/V: Krechmer chamber = 8 m3 (A/V~3/m); SAPHIR = 270 m3 (A/V~1/m) -> factor ~1/3.
CW0=16.0          # Krechmer prefactor (ug/m3)
CW_N=0.6          # Krechmer volatility exponent
SV_RATIO=1.0/3.0  # A/V_SAPHIR / A/V_Krechmer (~1/3); tune with the real surface areas
CW_SCALE=1.0      # <<< SINGLE TUNING KNOB for wall capacity Cw.  Cw_X = CW_SCALE*CW0*SV_RATIO*max(1,C*^CW_N).
                  #     RAISE it to push SVOC+LVOC to the walls (less semi-volatile SOA on the aerosol;
                  #     ELVOC is C*-trapped and barely responds). e.g. 1 -> 3 -> 10 to match AMS organics.
# --- rate-constant knobs (s-1); each is written into the include as CONST(...) ---
KCOND =5.5E-3     # condensation sink CS, tau_CS ~ 3 min (SAPHIR, held constant)
KFWALL=1.0E-3     # gas->wall deposition -- vapor-wall-loss SPEED knob
KPWALL=5.0E-5     # <<< KNOB: particle wall loss (irrev.). SAPHIR seed sulfate (2026-07-02) -> tau~20h => set ~1.4E-5.
# --- inert injected SEED aerosol: non-volatile organic ELVOC, injected in the .easy via q_seed. Here it only gets
#     its FATE (particle->wall + dilution), a COA mass term, and inclusion in the ELVOC element sums. It carries NO
#     gas-phase / condensation / evaporation reaction (never in the gas). Harmless when not injected (SOA_SEED=0).
SEED=('SEED',280.20,10,9)   # (name, MW g/mol, nC, nO) = C10H16O9 limonene-type ELVOC; set to () to disable entirely

def nN_of(fo):
    d={}
    for el,ct in re.findall(r'(Cl|[CHONS])(\d*)',fo): d[el]=d.get(el,0)+(int(ct) if ct else 1)
    return d.get('N',0)

def read_condensable(param=PARAM):
    """Return list of (name, MW, nC, nO, nN, Cstar_str, group) for condensable closed-shell species."""
    sp=[]
    for l in open(param,encoding='utf-8'):
        if l.startswith(';') or l.startswith('name\t'): continue
        p=l.rstrip('\n').split('\t')
        if len(p)<11: continue
        name,fo,MW,nC,nO,cstar,grp,shell=p[0],p[1],p[2],p[4],p[5],p[6],p[7],p[8]
        if shell!='closed' or grp not in ('SVOC','LVOC','ELVOC') or not MW or not nC: continue
        sp.append((name,float(MW),int(nC),int(nO),nN_of(fo),cstar,grp))
    return sp

def _joined(terms,per=8,ind="        "):
    if not terms: return "0.0"
    segs=[]; row=[]
    for t in terms:
        row.append(t)
        if len(row)==per: segs.append(" + ".join(row)); row=[]
    if row: segs.append(" + ".join(row))
    return (" $\n"+ind+"+ ").join(segs)

SUM_CHUNK=500   # terms per CONST(...) sum: Python's compiler hits "maximum recursion depth" at ~3000 chained '+' terms
                # (found 2026-11-08 with 3190 GECKO-A octane SOA species) -> long sums are split into partial sums

def _sum_lines(name, terms):
    """`name = CONST( t1 + t2 + ... )`, split into name_1, name_2, ... partial sums when longer than SUM_CHUNK."""
    if len(terms)<=SUM_CHUNK:
        return [f"{name} = CONST( {_joined(terms)} )"]
    L=[]; subs=[]
    for k in range(0,len(terms),SUM_CHUNK):
        sub=f"{name}_{k//SUM_CHUNK+1}"; L.append(f"{sub} = CONST( {_joined(terms[k:k+SUM_CHUNK])} )"); subs.append(sub)
    L.append(f"{name} = CONST( {' + '.join(subs)} )")
    return L

def _coa_lines(coa_floor, terms):
    """COA = floor + 1.6606E-12*(sum MW_X*SOA_X), the sum chunked like _sum_lines (COA_1, COA_2, ...)."""
    if len(terms)<=SUM_CHUNK:
        return [f"COA = CONST( {coa_floor:g} + 1.6606E-12*( {_joined(terms)} ) )"]
    L=[]; subs=[]
    for k in range(0,len(terms),SUM_CHUNK):
        sub=f"COA_{k//SUM_CHUNK+1}"; L.append(f"{sub} = CONST( {_joined(terms[k:k+SUM_CHUNK])} )"); subs.append(sub)
    L.append(f"COA = CONST( {coa_floor:g} + 1.6606E-12*( {' + '.join(subs)} ) )")
    return L

# GECKO-A's own gas-wall module (boxmodel4gecko, reactions.dum WIN/WOU lines, read 2026-11-08): k_in = 1.111E-3 s-1 for
# every species and ONE wall capacity for every species, Cw/(Mw*gamma) = 1.2E-4 mol/m3 (Matsunaga & Ziemann 2010 bag,
# A/V ~3.4 m-1) -> Cw_X = 1.2E-4*MW*1e6 ug/m3 = 5-18 mg/m3, independent of C*.  cw_mode='gecko' reproduces that.
GECKO_CW_MOL=1.2E-4      # mol/m3
GECKO_KFWALL=1.111E-3    # s-1 (1/900 s)

def emit_soa(sp, out_path, cw_scale=None, kcond=None, kfwall=None, kpwall=None, seed=None, coa_floor=0.5, cw_mode='krechmer'):
    """Write the SOA condensation include. Tuning knobs default to the module constants:
      cw_scale -> Cw wall capacity (SVOC+LVOC particle<->wall);  kpwall -> particle wall loss (s-1);
      kfwall -> gas->wall speed (s-1);  kcond -> condensation sink (s-1).
      coa_floor -> non-volatile COA baseline (ug/m3) added to the live organic mass (default 0.5).
      seed -> (name,MW,nC,nO) inert injected seed (defaults to module SEED); pass () to disable."""
    if cw_scale is None: cw_scale=CW_SCALE
    if kcond  is None: kcond =KCOND
    if kfwall is None: kfwall=KFWALL
    if kpwall is None: kpwall=KPWALL
    _seed = SEED if seed is None else seed
    if not _seed: _seed=None
    L=["; ==============================================================================",
       "; SOA_FZJ_2026.include  -- GENERATED by make_soa_condensation.py (FZ Juelich, 2026).",
       "; Reversible gas <-> particle <-> wall partitioning (Krechmer, Day & Jimenez 2020, ES&T 54:12890).",
       f"; {len(sp)} condensable species (SVOC+LVOC+ELVOC, closed-shell, C*<300). 6 reactions each.",
       ";   SOA_X = suspended aerosol ; WALL_X = gas-deposited wall organic (both reversible).",
       ";   Particle wall loss is IRREVERSIBLE -> per-species sink WALL_PL_X (+ WALL_PL_<group>_{C,O,N} diagnostic sums).",
       "; C* alone sets reversibility (tiny C* -> evaporation ~0 -> effectively irreversible).",
       "; INCLUDE AFTER the gas-phase mechanism and SIM_HOM.saphir.incl (needs DILUTE + the gas species).",
       "; ==============================================================================","",
       f"KCOND  = CONST({kcond:.2E})   ; condensation sink CS (s-1), tau_CS ~ 3 min (SAPHIR, held constant)",
       f"KFWALL = CONST({kfwall:.2E})   ; gas-wall deposition kf (s-1) -- vapor-wall-loss speed, TUNE to obs",
       f"KPWALL = CONST({kpwall:.2E})   ; particle wall loss (s-1) -- KNOB; SAPHIR seed => tau~20h ~ 1.4E-5",
       f"CW_SCALE = CONST({cw_scale:g})   ; wall-capacity scale -- RUNTIME KNOB (override in .easy). Cw_X = CW_SCALE*{CW0*SV_RATIO:.4g}*max(1,C*^{CW_N}) ug/m3",
       "; Cw (wall absorbing capacity, ug/m3) is C*-DEPENDENT -> embedded PER-SPECIES in each wall-evap",
       ";   rate below as a literal, not a single constant. Krechmer et al. 2016 (ES&T 50,5757, Table 1):",
       f";     Cw = {CW0:.0f}*(C*)^{CW_N} for C*>1 ; Cw = {CW0:.0f} for C*<1 (plateau). Cw is per-VOLUME so",
       f";     it scales with A/V: Krechmer chamber 8 m3 (A/V~3/m) -> SAPHIR 270 m3 (A/V~1/m), x{SV_RATIO:.3f}.",
       f";     => Cw_X = {cw_scale*CW0*SV_RATIO:.3g}*max(1, C*^{CW_N}) ug/m3 (plateau {cw_scale*CW0*SV_RATIO:.3g} below C*=1).  CW_SCALE={cw_scale:g}",
       ";        (SINGLE tuning knob in make_soa_condensation.py; raise to push SVOC+LVOC to walls, ELVOC ~unchanged).",
       ";   ATTRIBUTION: Krechmer 2016 gives Cw(C*) for THEIR 8 m3 chamber and explicitly scales only the",
       ";     TIME SCALE tau_GWE with A/V (weakly, ~20%, eddy-diffusion-compensated). The A/V scaling of Cw",
       ";     here is OUR physical inference (Cw = per-volume absorbing capacity ~ wall-Teflon/volume = A/V;",
       ";     thermodynamic, no ke compensation -> strong/linear). Cite Krechmer for Cw(C*), justify A/V separately.",
       f";   *** This C*^{CW_N} rise IS ACTIVE: it is Krechmer's SVOC wall-condensation amplification -- higher C*",
       ";       sees a larger wall reservoir -> MORE SVOC partitions to the walls -> higher gas-phase loss.",
       ";       ELVOC/LVOC (C*<1) sit on the plateau; only C*>1 (SVOC) are amplified. NOT simplified away. ***","",
       f"; --- COA: organic-aerosol mass (ug/m3) = cOA_start({coa_floor:g}) + 1.6606e-12 * sum(MW_X * SOA_X) ---",
       *_coa_lines(coa_floor, ([f'{_seed[1]:.3f}*SOA_{_seed[0]}'] if _seed else [])+[f'{MW:.3f}*SOA_{n}' for (n,MW,*_ ) in sp]),"",
       "; ============================ per-species partitioning (6 reactions each) ============================"]
    for (n,MW,nC,nO,nN,cstar,grp) in sp:
        SOA,WALL=f"SOA_{n}",f"WALL_{n}"
        if cw_mode=='gecko':
            cw_base=GECKO_CW_MOL*MW*1e6                          # GECKO-A/M&Z: constant per species, no C* dependence (A/V 3.4 bag)
        else:
            cw_base=CW0*SV_RATIO*max(1.0, float(cstar)**CW_N)   # CW_SCALE=1 wall capacity; runtime-scaled by the CW_SCALE variable below
        WPL=f"WALL_PL_{n}"   # per-species deposited (irreversible). Was element pools nC*WALL_PL_grp_C + ...
        L+=[f"{n} --> {SOA}",                    f"k[{n}-->{SOA}]=CONST(KCOND)",
            f"{SOA} --> {n}",                    f"k[{SOA}-->{n}]=CONST(KCOND*{cstar}/COA)",
            f"{n} --> {WALL}",                   f"k[{n}-->{WALL}]=CONST(KFWALL)",
            f"{WALL} --> {n}",                   f"k[{WALL}-->{n}]=CONST(KFWALL*{cstar}/({cw_base:.4g}*CW_SCALE))",
            f"{SOA} --> {WPL}",                  f"k[{SOA}-->{WPL}]=CONST(KPWALL)",
            f"{SOA} + DIL --> ",                 f"k[{SOA}+DIL-->]=CONST(DILUTE)"]
    if _seed:
        sn=_seed[0]
        L+=["","; --- inert injected SEED (SOA_"+sn+"): non-volatile organic ELVOC; particle->wall + dilution ONLY.",
            ";     Inject it in the .easy:  q_seed=input(q_seed) ; q_seed-->SOA_"+sn+f" ; k=CONST(1.0/(1.6606E-12*{_seed[1]:.2f}*300.))  (q_seed in ug/m3).",
            f"SOA_{sn} --> WALL_PL_{sn}",   f"k[SOA_{sn}-->WALL_PL_{sn}]=CONST(KPWALL)",
            f"SOA_{sn} + DIL --> ",         f"k[SOA_{sn}+DIL-->]=CONST(DILUTE)"]
    L+=["","; ====================== diagnostic element sums (9 SOA + 9 WALL + 9 WALL_PL) ======================"]
    for pref in ("SOA","WALL","WALL_PL"):
        for grp in ("SVOC","LVOC","ELVOC"):
            for el in ("C","O","N"):
                terms=[]
                for (n,MW,nC,nO,nN,cstar,g) in sp:
                    if g!=grp: continue
                    cnt={'C':nC,'O':nO,'N':nN}[el]
                    if cnt>0: terms.append(f"{cnt}*{pref}_{n}")
                if _seed and grp=='ELVOC' and pref in ('SOA','WALL_PL'):    # inert seed counts as ELVOC (C/O only)
                    _sc={'C':_seed[2],'O':_seed[3],'N':0}[el]
                    if _sc>0: terms=[f"{_sc}*{pref}_{_seed[0]}"]+terms
                L+=_sum_lines(f"{pref}_{grp}_{el}", terms)
    open(out_path,'w',encoding='utf-8').write("\n".join(L)+"\n")
    return {'species':len(sp),'phase_species':2*len(sp),'reactions':6*len(sp)}

if __name__=='__main__':
    from collections import Counter
    sp=read_condensable()
    st=emit_soa(sp,OUT)
    print(f"species: {st['species']} condensable -> {st['phase_species']} phase species + per-species WALL_PL_X")
    print(f"reactions: {st['reactions']} ; s-entries: COA + 27 element sums (9 SOA + 9 WALL + 9 WALL_PL)")
    print("by group:",dict(Counter(g for *_,g in sp)))
    print("wrote",OUT)
