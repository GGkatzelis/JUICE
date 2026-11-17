function J = scene_T1_J(Met)
% placeholder J function: scene_T1 has no photolysis, but F0AM's InitializeChemistry needs one.
J.J0 = zeros(size(Met.T(:)));
