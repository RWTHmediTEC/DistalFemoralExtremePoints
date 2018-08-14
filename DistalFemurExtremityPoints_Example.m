clearvars; close all; clc; opengl hardware;
% [List.f, List.p] = matlab.codetools.requiredFilesAndProducts([mfilename '.m']); 
% List.f = List.f'; List.p = List.p';

addpath(genpath([fileparts([mfilename('fullpath'), '.m']) '\src']));
addpath('data');

FileName = dir('data\*.mat');

Sigma.Medial = 4;
Sigma.Intercondylar = 6;
Sigma.Lateral = 4;

Vis = 1; % Visualization

for f=1 %:length(FileName)
load(FileName(f).name); 

Vertices = transformPoint3d(Vertices, USPTFM);

ExPoints = DistalFemurExtremityPoints(Vertices, Faces, Side, PFEA, Sigma, Vis);

end