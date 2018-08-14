clearvars; close all; clc

USP_Path = '..\..\UnifiedSagittalPlane\';
addpath(genpath(USP_Path));

pnl = [USP_Path, 'ExampleData\'];
fn = dir([pnl, '*.mat']);

for f=1:length(fn)
load([pnl, fn(f).name]); 

[USPTFM, PFEA, CEA] = USP(Vertices, Faces, Side, InitialRot, 'Subject', Subject, 'Visualization', false, 'Verbose', true);

save(fn(f).name, 'Vertices', 'Faces', 'Side', 'USPTFM', 'PFEA'); 

clearvars -except fn pnl

end


