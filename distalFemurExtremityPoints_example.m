clearvars; close all; clc; opengl hardware;

addpath(genpath([fileparts([mfilename('fullpath'), '.m']) '\src']));
addpath('data');

FileName = dir('data\*.mat');

Idx = 1;

load(FileName(Idx).name); 
distalFemur.vertices = Vertices;
distalFemur.faces = Faces;
distalFemurUSP = transformPoint3d(distalFemur, USPTFM);

ExPoints = distalFemurExtremityPoints(distalFemurUSP, Side, PFEA, 'visu', 1, 'debug', 0);

% [List.f, List.p] = matlab.codetools.requiredFilesAndProducts([mfilename '.m']); 
% List.f = List.f'; List.p = List.p';