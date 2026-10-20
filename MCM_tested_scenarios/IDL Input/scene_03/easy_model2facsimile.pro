; Copyright (c) 1998, Forschungszentrum Juelich GmbH ICG-3
; All rights reserved.
; Unauthorized reproduction prohibited.
;
;+
; NAME:
;   easy_model2facsimile
;
; PURPOSE:
;   This function returns a string array with html code created from
;   the entries of the model structure
;
; CATEGORY:
;   PROG_TOOLS/STRINGS
;
; CALLING SEQUENCE:
;   Result = easy_model2facsimile(Str)
;
; INPUTS:
;   Str:  structure (see easy_userfile2model)
;
; KEYWORD PARAMETERS:
;   ERROR: set to zero if successfully done
;
; OUTPUTS:
;   This function returns the string array with html code
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;   Written by: Theo Brauers, Franz Rohrer
;   1999 Jul 20 IMPORTANT MODIFICATIONS:
;               1.  call of parset 3 time in order to provide the correct settings if
;                   settings are functions of other settings
;               2.  change in time management previously an error occured with timesteps
;                   different from 10 (?). This error was partly caused by the poor ducumentation
;   1999 Jul 23 IMPORTANT MODIFICATIONS:
;               1.  remove CALL init2 from COMPILE initial
;               2.  move CALL table up to 1st line in COMPILE  print
;   1999 AUG 6 IMPORTANT MODIFICATIONS:
;               1.  new version for steady state.
;   1999 NOV 11 IMPORTANT MODIFICATIONS:
;               1. in compile in: out=[out, '  time=tnull' ]
;   1999 NOV 25 IMPORTANT MODIFICATIONS:
;               output flag in time dep calculations
;   2003 NOV 06 Zeile 329 um 3 'COMPILE parset' erweitert, siehe n�chste drei Zeilen
;# compile equations was already started as buffer for commands
;# at the compile parset location
;#out=[out, 'COMPILE equations', out_eq,'CALL parset','CALL parset','CALL parset']
;   2004 JUN 25 FR von out=[out, 'COMPILE equations', out_eq,'CALL parset','CALL parset','CALL parset']
;            zur�ck zu out=[out, 'COMPILE equations', out_eq]
;  FR 2026-05-06 Eintrag: out=[out, ' CALL parset']
;
;-
FUNCTION easy_model2facsimile, m, ERROR=err

err = 1
; --- outstep_follow flag ---------------------------------------------------
; 1 (default): concentration-dependent k-values (KMT06, k_O1D+H2O, ...) are
;              re-evaluated from the current ODE state at every outstep.
;              Physically correct when ODE species (e.g. H2O via dilution)
;              change between ENZ rows.
; 0 (old):     settings only updated when an ENZ row is read (CALL in).
;              k-values stay frozen between ENZ boundaries.
outstep_follow = 1
; ---------------------------------------------------------------------------

; write comment header
COMMON easy_message_var, easy_mv
out=['* Facsimile-Script / '+easy_mv.version_ident, '* '+SYSTIME(),'']

; get the number of variables
a = WHERE(m.substances.type EQ 'V',nvar)
IF nvar LE 0 THEN RETURN, [out, 'nvar=0']
out=[out,'INTEGER #nvar '+easy_n2s(nvar)]

; number of boxes
nbox=m.geometry.nbox
out=[out,'INTEGER #nbox ' +easy_n2s(nbox)]

; number of reactions
nreac=N_ELEMENTS(m.reactions)
IF nreac LE 0 THEN RETURN, [out, 'nreac=0']
out=[out,'INTEGER #nreac '+easy_n2s(nreac)]

; number of settings
IF (SIZE(m.settings))[0] NE 0 THEN BEGIN
    nset=N_ELEMENTS(m.settings)
    out=[out,'INTEGER #nset ' +easy_n2s(nset)]
ENDIF ELSE nset=0

; string for reading the first line of dat files
out=[out,'TEXT #title "TITLE"','']

; facsimile parameter for the integration order
out=[out, 'KMAX 3']
; facsimile parameter for # of retries after numerical problems
out=[out, 'SETRS 25', '']


out=[out,'VARIABLE <#nvar,#nbox> c '+easy_n2s(/FAC, nvar*nbox)+'1.']
out=[out,'PARAMETER <'+easy_n2s(LONG(nvar*nbox*6.+10))+'> worksp']
out=[out,'PARAMETER <#nvar,#nbox> sc '+easy_n2s(/FAC, nvar*nbox)+'1.']
out=[out,'PARAMETER <#nvar,#nbox> p '+easy_n2s(/FAC, nvar*nbox)+'0.']
out=[out,'PARAMETER <#nvar,#nbox> d '+easy_n2s(/FAC, nvar*nbox)+'0.']

; definition transportparameter
IF m.geometry.nbox GT 1 THEN BEGIN
    out=[out,easy_model2transport_code(m, CODE=0, ERR=err2)]
    IF err2 THEN RETURN, [out, 'transport code failure']
ENDIF
;
out=[out,'PARAMETER <#nreac,#nbox> v '+easy_n2s(/FAC, nreac*nbox)+'0.']
out=[out,'PARAMETER <#nreac,#nbox> k '+easy_n2s(/FAC, nreac*nbox)+'1.e-20']

; setting are required
IF nset GT 0 THEN $
    out=[out,'PARAMETER <#nset,#nbox> s '+easy_n2s(/FAC, nset*nbox)+'1.']

; h is the time increment
out=[out,'PARAMETER h']

; tread is the time which was read from file and corrected with the offset
out=[out,'PARAMETER tread -1E25']

;
out=[out,'PARAMETER zeitflag ']
out=[out,'PARAMETER zeitschwelle']

; this parameter holds the sign for the whenever procedure: required since
; each call needs adifferent gradient, work-around for horror facsimile
; new set -1 after change in code 20 jul 1999
out=[out,'PARAMETER zsgn -1']
IF DOUBLE(m.facs.STEADY_STATE_TIME) LE 0 THEN $
    out=[out,'PARAMETER zeitintervall '+easy_n2s(m.facs.outstep)] ELSE $
out=[out,'PARAMETER zeitintervall '+easy_n2s(m.facs.STEADY_STATE_TIME)+' tfile 1E25']

; temporary variable to hold corrected time for output
out=[out,'PARAMETER tout']

; time offset which is set to the first available input time: required since
; negative values for time do not reliably work as input to facsimile
IF DOUBLE(m.facs.STEADY_STATE_TIME) GT 0 THEN $
    out=[out,'PARAMETER tnull -1E25'] ELSE $
    out=[out,'PARAMETER tnull 1E25']

; define the number of input columns to be read from file
out=[out,'PARAMETER <'+easy_n2s(m.facs.datacol)+'>input 1E25']

; stops execution after last input line
out=[out,'PARAMETER fertig 0.0']
out=[out,'']


;compile initial
out=[out,'COMPILE initial']
out=[out,' HMAX '+m.facs.hmax]
IF m.facs.monitor EQ '1' THEN BEGIN
    out=[out,' OPEN 17 "monitor.txt"']
ENDIF
out=[out,' FOR #99='+easy_n2s(m.facs.usubr_msg)]

; enz file read
out=[out, ' OPEN 9 "'+m.files.enz+'"', ' DO 21 FOR #88=9', '  ENTER', ' LABEL 21' ]

; write a status line to the facsimile out file
out=[out,' WRITE 1=6, "start reading from input file"']
out=[out,' CALL in']
out=[out,' CALL init2']

; time is a pre defined variable for the integration
IF DOUBLE(m.facs.STEADY_STATE_TIME) LE 0 THEN BEGIN
    out=[out,' time = tread']
ENDIF


out=[out,' CALL zeitsprung']

; ignore the the secon call to in if steady state is requested
IF DOUBLE(m.facs.STEADY_STATE_TIME) LE 0 THEN BEGIN
    out=[out,' CALL in']
ENDIF

; open the header file and assign unit number
out=[out,' OPEN 14 "'+m.files.hd+'"']
out=[out,' DO 16 FOR #88=14']

; open the cs file and assign unit number and copy from hd to cs
out=[out,'  OPEN 10 "'+m.files.cs+'"']
out=[out,'  DO 12 FOR #87=10', '    ENTER', '  LABEL 12']

; open the kv file and assign unit number and copy from hd to kv
out=[out,'  OPEN 11 "'+m.files.kv+'"']
out=[out,'  DO 13 FOR #87=11', '    ENTER', '  LABEL 13']

; open the pd file and assign unit number and copy from hd to pd
out=[out,'  OPEN 12 "'+m.files.pd+'"']
out=[out,'  DO 14 FOR #87=12', '    ENTER', '  LABEL 14']

; open a transport parameter out files
IF m.geometry.nbox GT 1 THEN BEGIN
    out=[out,'  OPEN 13 "'+m.files.tr+'"']
    out=[out,'  DO 15 FOR #87=13', '    ENTER', '  LABEL 15']
ENDIF

; close the header file
out=[out,' LABEL 16']
out=[out,' CLOSE 14']
out=[out,'**', '']

; compile parset
; must appear before compile instant
out=[out,'COMPILE parset']


; save all command which contain variables in array for compile equation
out_eq=['* from compile parset']
IF (SIZE(m.settings))[0] NE 0 THEN BEGIN
    FOR i=0L, N_ELEMENTS(m.settings)-1 DO BEGIN
        IF m.settings[i].name EQ 'M' THEN BEGIN
            out=[out,'* convert to mixing ratios']
            out=[out,' DO 50 FOR #1=0 (1) #nbox-1']
            out=[out,' ARRAY <#nvar> worksp']
            out=[out,' c<*,#1>=c<*,#1>/'+replace_string(m.settings[i].fc_name,'*','#1')]
            out=[out,' ARRAY end']
            out=[out,' LABEL 50']
        ENDIF
        FOR j=0,m.geometry.nbox-1 DO BEGIN
            IF m.settings[i].scale NE '' THEN $
                out2 = replace_string(m.settings[i].fc_call[j], '=','=(')+')*'+m.settings[i].scale ELSE $
                out2 = m.settings[i].fc_call[j]
            out = [out, out2]
            IF STRPOS(str_sep2(m.settings[i].fc_call[j], '=', ONLY=1), 'c<') GE 0 THEN out_eq=[out_eq, out2]
        ENDFOR
        IF m.settings[i].name EQ 'M' THEN BEGIN
            out=[out,'* convert back to conc']
            out=[out,' DO 51 FOR #1=0 (1) #nbox-1']
            out=[out,' ARRAY <#nvar> worksp']
            out=[out,' c<*,#1>=c<*,#1>*'+replace_string(m.settings[i].fc_name,'*','#1')]
            out=[out,' ARRAY end']
            out=[out,' LABEL 51']
        ENDIF
    ENDFOR
ENDIF


;Belegung der Transportparameter
IF m.geometry.nbox GT 1 THEN BEGIN
    out=[out,easy_model2transport_code(m, ERR=err2, CODE=3)]
    IF err2 THEN RETURN, [out, 'transport code failure']
    out=[out,easy_model2transport_code(m, Err=err2, CODE=1)]
    IF err2 THEN RETURN, [out, 'transport code failure']
ENDIF

out=[out,'ARRAY <#nbox> worksp']
FOR i=0,N_ELEMENTS(m.reactions)-1 DO IF m.reactions[i].scale NE '' THEN $
    out=[out, replace_string(m.reactions[i].fc_call, '=','=(')+')*'+m.reactions[i].scale] ELSE $
    out=[out, m.reactions[i].fc_call ]
out=[out,'ARRAY end']
out=[out,'**', '']

;init 2 was compile instant
out=[out,'COMPILE init2']
out=[out,' CALL parset']
out=[out,' CALL parset']
out=[out,' CALL parset']
FOR i=0L, N_ELEMENTS(m.substances)-1 DO $
    IF m.substances[i].type EQ 'V' AND m.substances[i].scale NE '' THEN $
        FOR j=0L,m.geometry.nbox-1 DO out=[out, ' s'+str_sep2(m.substances[i].fc_init[j],'=', ONLY=0)+'='+m.substances[i].scale ]
FOR i=0L, N_ELEMENTS(m.substances)-1 DO $
    IF m.substances[i].type EQ 'V' THEN $
        FOR j=0L,m.geometry.nbox-1 DO out=[out, ' '+m.substances[i].fc_init[j]]
out=[out,' DO 47 FOR #1=0 (1) #nbox-1' $
        ,'  ARRAY <#nvar> worksp' $
        ,'  c<*,#1>=c<*,#1>/sc<*,#1>' $
        ,'  ARRAY end'  $
        ,' LABEL 47'   ]
out=[out,'**', '']

; compile print
; this is completely new code after 1999 Jul 20
; FR 2026-05-06: idx needed here for ENTER logic (time-dep. mode)
;
; Two checks suppress unwanted output rows without a flag variable:
;   IF(fertig)     : skip write when WHENEVER fertig=1 fires (last duplicate row)
;   IF(tout-tnull) : skip write when time=0 (first row = initial conditions)
;                    tout = time+tnull is computed first so the check is current.
;
; FACSIMILE IF syntax: IF (expr) label_neg, label_zero, label_pos
idx=WHERE(STRUPCASE(m.input.name) EQ 'TIME', cdx)
IF cdx EQ 0 THEN RETURN, [out, 'time not defined']
out=[out,'COMPILE print']
; Skip write if fertig=1 (WHENEVER call -> no duplicate last row)
out=[out, ' IF (fertig) *,60,64']
out=[out, ' LABEL 60']
; Update tout so next check reflects current integration time
out=[out, ' tout = time + tnull']
; Skip write if time=0 (first call before any integration -> no initial-conditions row)
out=[out, ' IF (tout-tnull) *,64,63']
out=[out, ' LABEL 63']
out=[out, ' CALL table']
out=[out, ' LABEL 64']
IF DOUBLE(m.facs.STEADY_STATE_TIME) LE 0 THEN $
    out=[out, ' IF (tread-time) *,*,65']
out=[out, ' CALL in']
out=[out, ' LABEL 65']
IF outstep_follow EQ 1 THEN BEGIN
    ; FR 2026-05-06: Update only concentration-dependent settings and k-values.
    ; Uses s<> (set at last ENZ boundary) and c<> (current ODE state).
    ; NOT input<> -> ENZ-driven quantities (HV, T, DIL) are not touched.
    FOR i=1, N_ELEMENTS(out_eq)-1 DO out=[out, ' '+out_eq[i]]
    out=[out, ' ARRAY <#nbox> worksp']
    FOR i=0,N_ELEMENTS(m.reactions)-1 DO $
        IF m.reactions[i].scale NE '' THEN $
            out=[out, ' '+replace_string(m.reactions[i].fc_call,'=','=(')+')*'+m.reactions[i].scale] $
        ELSE $
            out=[out, ' '+m.reactions[i].fc_call]
    out=[out, ' ARRAY end']
ENDIF
out=[out, ' CALL zeitsprung']
out=[out, ' zsgn = -zsgn']
out=[out, ' RESTART']
out=[out,'**', '']

;compile table
out=[out,'COMPILE table']
IF DOUBLE(m.facs.STEADY_STATE_TIME) GT 0 THEN $
    out=[out,' tout = tnull + tfile'] ELSE $
    out=[out,' tout = time + tnull']

; convert w/ scale factor
out=[out,' DO 48 FOR #1=0 (1) #nbox-1' $
        ,'  ARRAY <#nvar> worksp' $
        ,'  c<*,#1>=c<*,#1>*sc<*,#1>' $
        ,'  ARRAY end'  $
        ,' LABEL 48'   ]

out=[out,' DO 50 FOR #50=#nvar,#51=#nbox,#52=#nset,#53=#nbox,#54=10']
out=[out,'  ENTER tout c s ']
out=[out,' LABEL 50']

; convert w/ scale factor
out=[out,' DO 49 FOR #1=0 (1) #nbox-1' $
        ,'  ARRAY <#nvar> worksp' $
        ,'  c<*,#1>=c<*,#1>/sc<*,#1>' $
        ,'  ARRAY end'  $
        ,' LABEL 49'   ]

out=[out,' DO 51 FOR #50=#nreac,#51=#nbox,#52=#nreac,#53=#nbox,#54=11']
out=[out,'  ENTER tout k v ']
out=[out,' LABEL 51']
out=[out,' DO 52 FOR #50=#nvar,#51=#nbox,#52=#nvar,#53=#nbox,#54=12']
out=[out,'  ENTER tout p d ']
out=[out,' LABEL 52']
IF m.geometry.nbox GT 1 THEN BEGIN
    out=[out,' DO 53 FOR #50=#nvar,#51=#nbox,#52=#nvar,#53=#nbox,#54=13']
    out=[out,'  ENTER tout tr mxr ']
    out=[out,' LABEL 53']
ENDIF
out=[out,' IF <#98> 54,54,*']
out=[out,'  fertig=1']
out=[out,' LABEL 54']
out=[out,'**', '']

; compile general
out=[out,'COMPILE general']
out=[out,' zeitflag=(zeitschwelle-time)*zsgn']
out=[out,'**', '']

; compile equations was already started as buffer for commands
; at the compile parset location
out=[out, 'COMPILE equations', out_eq]

;Transportmatrix berechnen
IF m.geometry.nbox GT 1 THEN BEGIN
    out=[out,easy_model2transport_code(m, ERR=err2, CODE=2)]
    IF err2 THEN RETURN, [out, 'transport code failure']
ENDIF

out=[out,' ARRAY <#nbox> worksp']
out=[out, easy_model2facs_pd_code(m)]
out=[out,' ARRAY <#nvar,#nbox> worksp']

; include the transport matrix
IF m.geometry.nbox GT 1 THEN str='+tr' ELSE str=''

; the equation for the time derivative
out=[out," 'c=p-d"+str]

out=[out,' ARRAY end']
out=[out,'**']
out=[out,'']

;compile zeitsprung
out=[out,'COMPILE zeitsprung']
IF DOUBLE(m.facs.STEADY_STATE_TIME) LE 0 THEN BEGIN
    out=[out,' h=zeitschwelle+zeitintervall']
    out=[out,' zeitschwelle=AMIN(h,tread)']
ENDIF ELSE $
    out=[out,' zeitschwelle=tnull+zeitintervall']
out=[out,'**']
out=[out,'']

;compile in
; this is the input reader from unit # 9 which is the
; .dat or .enz file opened, at fertig the input unit is closed
out=[out, 'COMPILE in']
; idx already defined before COMPILE print

IF DOUBLE(m.facs.STEADY_STATE_TIME) LE 0 THEN BEGIN
    out=[out, '  CALL parset']
    out=[out, '  CALL parset']
    out=[out, '  CALL parset']
    out=[out, '  DO 97 FOR #49='+easy_n2s(m.facs.datacol)+',#88=9']
    out=[out, '    ENTER input', '  LABEL 97']
    out=[out, '  IF (tnull-1.e25) 98, *, *' $
            , '  tnull'+m.input[idx[0]].fc_call $
            , '  LABEL 98']
    out=[out, '  WRITE 1=6, "t(input)=" ,'+str_sep2(m.input[idx[0]].fc_call,'=', ONLY=1)+' %']
    out=[out, '  IF ('+str_sep2(m.input[idx[0]].fc_call,'=', ONLY=1)+'-tnull) *,99,99']
    out=[out, '  fertig=1']
    out=[out, '  CLOSE 9', '  CALL table', '  STOP']
    out=[out, '  LABEL 99']
    out=[out, '  tread'+m.input[idx[0]].fc_call+'-tnull']
ENDIF ELSE BEGIN
    out=[out, '  DO 97 FOR #49='+easy_n2s(m.facs.datacol)+',#88=9']
    out=[out, '    ENTER input', '  LABEL 97']
    out=[out, '  IF (tfile-1.E25) 98,*,*']
    out=[out, '  tfile'+m.input[idx[0]].fc_call ]
    out=[out, '  WRITE 1=6, "tfile=", tfile %']
    out=[out, '  LABEL 98']
    out=[out, '  IF ('+str_sep2(m.input[idx[0]].fc_call,'=', ONLY=1)+'-tread) *,99,99']
    out=[out, '  CLOSE 9', '  STOP']
    out=[out, '  LABEL 99']
    out=[out, '  WRITE 1=6, "t(input)=" ,'+str_sep2(m.input[idx[0]].fc_call,'=', ONLY=1)+' %']
    out=[out, '  tread'+m.input[idx[0]].fc_call ]
    out=[out, '  tnull'+m.input[idx[0]].fc_call+'-tfile' ]
    ; nov 11 T.B.
    out=[out, '  time=tnull' ]
    out=[out, '  tout=tnull+zeitinterval']
   ;out=[out, '  WRITE 1=6, "tnull=", tnull, "tout=", tout %']
    out=[out, '  CALL parset']
    out=[out, '  CALL parset']
    out=[out, '  CALL parset']
ENDELSE

out=[out,'**', '']

; define the end of facsimile execution:
; CALL print is called with zeitflag=0. Inside COMPILE print,
; IF(fertig)>0 jumps to LABEL 64 and skips CALL table -> no duplicate last row.
out=[out,'WHENEVER fertig=1 %']
out=[out,' zeitflag=0 CALL print']
out=[out,'**', '']
IF m.facs.monitor EQ '1' THEN BEGIN
    out=[out,' MONiTOR 17']
ENDIF
out=[out,'BEGIN','STOP', '']

; add a ';' at the end of each program line
a=WHERE(out NE '', c)
IF c GT 0 THEN out[a]=out[a]+';'

; break lines into pieces of maximum 70 characters
out=break_lines(STRCOMPRESS(out),70, LIST=[',',' ','/','(',')','|'])

; reset error flag
err = 0
RETURN, out
END
