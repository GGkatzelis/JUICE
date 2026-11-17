% scene_T1_mech.m - F0AM mechanism for F. Rohrer's scene_T1: X + X --> XX with k = c*X.
% F0AM computes rate constants before it integrates and allows two reactants, so a rate constant that
% follows a species cannot be written here. run_scene_T1.m sets k for each short step from F0AM's own
% X at the start of that step (global KT1_), so k is held per step and the steps are made short.

global KT1_

SpeciesToAdd = {'X'; 'XX'};
RO2ToAdd = {};
AddSpecies

i=i+1;
Rnames{i} = 'X+X=XX';
k(:,i) = KT1_;
Gstr{i,1} = 'X';
Gstr{i,2} = 'X';
fX(i)=fX(i)-2; fXX(i)=fXX(i)+1;
