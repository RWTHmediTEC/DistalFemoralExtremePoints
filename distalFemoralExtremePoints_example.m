clearvars; close all; opengl hardware

% Add src path
addpath(genpath([fileparts([mfilename('fullpath'), '.m']) '\src']));

%% Clone example data
if ~exist('VSD', 'dir')
    try
        !git clone https://github.com/RWTHmediTEC/VSDFullBodyBoneModels VSD
        rmdir('VSD/.git', 's')
    catch
        warning([newline 'Clone (or copy) the example data from: ' ...
            'https://github.com/RWTHmediTEC/VSDFullBodyBoneModels' newline 'to: ' ...
            fileparts([mfilename('fullpath'), '.m']) '\VSD' ...
            ' and try again!' newline])
        return
    end
end

%% Load subject names
Subjects = dir('data\*.mat');
Subjects = strrep({Subjects.name}','.mat','');
Subjects(1:2:20,2) = {'L'}; Subjects(2:2:20,2) = {'R'};

for s=1%:size(Subjects, 1)
    name = Subjects{s,1};
    side = Subjects{s,2};
    
    % Prepare distal femur
    load(['VSD\Bones\' name '.mat'], 'B');
    load(['data\' name '.mat'],...
        'inertiaTFM','uspInitialRot','distalCutPlaneInertia','USPTFM','PFEA');
    femurInertia = transformPoint3d(B(ismember({B.name}, ['Femur_' side])).mesh, inertiaTFM);
    distalFemurInertia = cutMeshByPlane(femurInertia, distalCutPlaneInertia);
    distalFemurUSP = transformPoint3d(distalFemurInertia, USPTFM);
    
    ExPoints = distalFemoralExtremePoints(distalFemurUSP, side, PFEA, 'visu',1, 'debug',0);
end

% [List.f, List.p] = matlab.codetools.requiredFilesAndProducts([mfilename '.m']); 
% List.f = List.f'; List.p = List.p';