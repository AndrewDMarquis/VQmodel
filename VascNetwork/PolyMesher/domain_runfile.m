clear; close all; clc;

[Node,Element,Supp,Load,P] = PolyMesher(@quarter_circle_domain,10000,100);

% save('SEED_reuleaux_triangle.mat')
save('SEED_quarter_circ.mat')
% save('SEED_half_circ.mat')
