%% run_scene_T1.m - F. Rohrer's scene_T1 in F0AM (see ../scene_T1.easy for the scenario and its solutions)
% X + X --> XX, k = c*X, c = 19375/(421*14400*1e18); X(0) = 1e9 molec/cm3; P = 1013 hPa;
% T = 280 K for the first hour, 300 K for the second.
%
% F0AM cannot follow X inside k (Sect. "F0AM" of ../README.md), so the run is a loop of F0AM calls of
% length dt: each call holds k = c*X(start of the call), with X from F0AM's previous call. The
% concentrations are passed from call to call as mixing ratios, so when T steps from 280 to 300 K F0AM
% converts them with its own density M = NumberDensity(P,T): the density step at 3600 s is F0AM's own.
% dt = 3600 s is "k held per row"; as dt shrinks the result should approach "k follows X".
% Writes scene_T1_F0AM.csv: dt, X and XX (molec/cm3) at 3600 s (before the step) and at 7200 s.

clear
global KT1_
here = fileparts(mfilename('fullpath'));
cd(here)

c  = 19375/(421*14400*1e18);
P  = 1013;
DTS = [3600 600 60 10];

ModelOptions.Verbose       = 0;
ModelOptions.EndPointsOnly = 1;
ModelOptions.LinkSteps     = 0;
ModelOptions.Repeat        = 1;
ModelOptions.SavePath      = 'DoNotSave';             % one call per step: do not fill F0AM's Runs folder
ChemFiles = {''; 'scene_T1_J(Met)'; 'scene_T1_mech'};
BkgdConc  = {'DEFAULT' 0};

res = zeros(numel(DTS), 5);
timer_ = tic;
for d = 1:numel(DTS)
    dt = DTS(d);
    ModelOptions.IntTime = dt;
    px = 1e9/NumberDensity(P, 280)*1e9; pxx = 0;     % ppb, carried from call to call (F0AM's own units)
    x3600 = NaN; xx3600 = NaN;
    for n = 0:round(7200/dt)-1
        t = n*dt;
        T = 280 + 20*(t >= 3600);                      % the ENZ row the call starts in
        M = NumberDensity(P, T);                       % F0AM's own density for this call
        KT1_ = c*px*M/1e9;                             % k = c*X, X in molec/cm3 at the start of the call
        Met = {'P' P; 'T' T; 'H2O' 0; 'kdil' 0};
        InitConc = {'X' px 0; 'XX' pxx 0};
        S = F0AM_ModelCore(Met, InitConc, ChemFiles, BkgdConc, ModelOptions);
        px = S.Conc.X(end); pxx = S.Conc.XX(end);
        if abs(t + dt - 3600) < 1e-9                   % end of the first hour, still at 280 K
            x3600 = px*M/1e9; xx3600 = pxx*M/1e9;
        end
    end
    x = px*NumberDensity(P, 300)/1e9; xx = pxx*NumberDensity(P, 300)/1e9;
    res(d,:) = [dt x3600 xx3600 x xx];
    fprintf('dt %5d s: X(3600) %.6e  XX(3600) %.6e  X(7200) %.6e  XX(7200) %.6e\n', res(d,:));
end
fprintf('F0AM run time %.1f s\n', toc(timer_));
fid = fopen('scene_T1_F0AM.csv', 'w');
fprintf(fid, 'dt_s,X_3600,XX_3600,X_7200,XX_7200\n');
fprintf(fid, '%g,%.7e,%.7e,%.7e,%.7e\n', res');
fclose(fid);
