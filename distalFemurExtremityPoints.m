function EP = distalFemurExtremityPoints(distalFemurUSP, Side, PFEA, varargin)
% TODO
% 
% INPUT:
%   - REQUIRED:
%     TODO
% 
%   - OPTIONAL:
%     TODO
% 
% OUTPUT:
%     TODO
% 
% EXAMPLE:
%     Run the file 'DistalFemurExtremityPoints_Example.m'
% 
% TODO:
%   - Parsing
%   - Revise documentation
%   - Use axes handles during visualization
% 
% AUTHOR: Maximilian C. M. Fischer
% 	mediTEC - Chair of Medical Engineering, RWTH Aachen University
% VERSION: 1.0.1
% DATE: 2018-08-14
% LICENSE: Modified BSD License (BSD license with non-military-use clause)

addpath(genpath([fileparts([mfilename('fullpath'), '.m']) '\src']));

%% Parse inputs
p = inputParser;
logParValidFunc=@(x) (islogical(x) || isequal(x,1) || isequal(x,0));
addRequired(p,'distalFemurUSP',@(x) isstruct(x) && isfield(x, 'vertices') && isfield(x,'faces'))
addParameter(p,'visualization',false,logParValidFunc);
addParameter(p,'debugVisualization',false,logParValidFunc);
parse(p,distalFemurUSP,varargin{:});
visu = logical(p.Results.visualization);
debugVisu = logical(p.Results.debugVisualization);

%% Settings
% Sigma: For explanation see the function: BOMultiScaleCurvature2D_adapted.m
sigmaStart = 1;
sigmaDelta = 1;
sigmaMedial = 4;
sigmaIntercondylar = 6;
sigmaLateral = 4;

%% Find boundary points of the condyles
% Intersection points between the mesh an the posterior foci elliptical axis 
PFEA_IS_Pts = intersectLineMesh3d(PFEA, distalFemurUSP.vertices, distalFemurUSP.faces);

% Should be always 4 points
if size(PFEA_IS_Pts,1) ~= 4
    error('PFEA intersection points with the distal femur mesh ~= 4')
end

% The intersection points devide the distal femur into 3 parts in Z direction:
% Negative Z = NZ
% Middle     = MZ
% Positive Z = PZ
switch Side
    case 'Right'
        % RIGHT knee: medial - neg. Z values (NZ), lateral - pos. Z values (PZ)
        SC(1).Zone = 'NZ';
        SC(2).Zone = 'MZ';
        SC(3).Zone = 'PZ';
    case 'Left'
        %  LEFT knee: medial - pos. Z values (PZ), lateral - neg. Z values (NZ)
        SC(1).Zone = 'PZ';
        SC(2).Zone = 'MZ';
        SC(3).Zone = 'NZ';
end

% Assign them to the zones
PISP{1,2} = ismember(PFEA_IS_Pts(:,3), min(PFEA_IS_Pts((PFEA_IS_Pts(:,3) < 0),3)));
PISP{2,2} = ismember(PFEA_IS_Pts(:,3), max(PFEA_IS_Pts((PFEA_IS_Pts(:,3) < 0),3)));
PISP{3,2} = ismember(PFEA_IS_Pts(:,3), min(PFEA_IS_Pts((PFEA_IS_Pts(:,3) > 0),3)));
PISP{4,2} = ismember(PFEA_IS_Pts(:,3), max(PFEA_IS_Pts((PFEA_IS_Pts(:,3) > 0),3)));

for i=1:size(PFEA_IS_Pts,1)
    PISP{i,3} = PFEA_IS_Pts(PISP{i,2},:);
end

% Start and end point of the zones, rounded to an integer value in Z direction
NZS = PISP{1,3}; NZS(3) = ceil(NZS(3));
MZS = PISP{2,3}; MZS(3) = ceil(MZS(3));
PZS = PISP{3,3}; PZS(3) = ceil(PZS(3));
NZE = PISP{2,3}; NZE(3) = floor(NZE(3));
MZE = PISP{3,3}; MZE(3) = floor(MZE(3));
PZE = PISP{4,3}; PZE(3) = floor(PZE(3));

% The number of sagittal cuts per zone
SC(1).NoC = abs(NZS(3)-NZE(3))+1;
SC(2).NoC = abs(MZS(3)-MZE(3))+1;
SC(3).NoC = abs(PZS(3)-PZE(3))+1;

SC(1).Origin = NZS;
SC(2).Origin = MZS;
SC(3).Origin = PZS;

LS = size(SC,2);
ZoneColors = parula(LS*2);

%% Sagittal Cuts (SC)
ZVector = [0, 0, 1];

Contours = cell(1,LS);
for s=1:LS
    SC(s).Color = ZoneColors(s,:);
    % Create cutting plane origins
    SC(s).PlaneOrigins = repmat(SC(s).Origin, SC(s).NoC, 1) + ...
        [zeros(SC(s).NoC,2), (0:SC(s).NoC-1)'];
    % Create SC(s).NoC Saggital Contour Profiles (SC(s).P)
    Contours{s} = IntersectMeshPlaneParfor(distalFemurUSP, SC(s).PlaneOrigins, ZVector);
    for c=1:SC(s).NoC
        % If there is more than one closed contour after the cut, use the longest one
        [~, IobC] = max(cellfun(@length, Contours{s}{c}));
        SC(s).P(c).xyz = Contours{s}{c}{IobC}';
        % Close contour: Copy start value to the end
        if ~isequal(SC(s).P(c).xyz(1,:),SC(s).P(c).xyz(end,:))
            SC(s).P(c).xyz(end+1,:) = SC(s).P(c).xyz(1,:);
        end
        % If the contour is sorted clockwise
        if sign(varea(SC(s).P(c).xyz')) == -1
            % Sort the contour counter-clockwise
            SC(s).P(c).xyz = flipud(SC(s).P(c).xyz);
        end
        % Set Start of the contour to the maximum Y value
        [~, IYMax] = max(SC(s).P(c).xyz(:,2));
        if IYMax ~= 1
            SC(s).P(c).xyz = SC(s).P(c).xyz([IYMax:length(SC(s).P(c).xyz),1:IYMax-1],:);
        end
    end
end

%% SagExPts = Sagittal Extremity Points
% A pattern-recognition algorithm for identifying the articulating surface
for s=1:LS
    for c=1:SC(s).NoC
        % Only the X & Y values are needed because the countours are parallel to the X-Y plane.
        Contour = SC(s).P(c).xyz(:,1:2);
        % Get the anterior and posterior extremity points of the articulating surface
        switch SC(s).Zone
            case 'NZ'
                [SC(s).P(c).A, SC(s).P(c).B, SC(s).P(c).H] = ...
                    SagExPts_MedCond(Contour, sigmaStart, sigmaDelta, sigmaMedial, debugVisu);
            case 'MZ'
                [SC(s).P(c).A, SC(s).P(c).B, SC(s).P(c).H] = ...
                    SagExPts_IntCond(Contour, sigmaStart, sigmaDelta, sigmaIntercondylar, debugVisu);
            case 'PZ'
                [SC(s).P(c).A, SC(s).P(c).B, SC(s).P(c).H] = ...
                    SagExPts_LatCond(Contour, sigmaStart, sigmaDelta, sigmaLateral, debugVisu);
        end
    end
end


%% Find global extremity points of the sagittal contours (SCs)
NZ = strcmp({SC(:).Zone}, 'NZ');
MZ = strcmp({SC(:).Zone}, 'MZ');
PZ = strcmp({SC(:).Zone}, 'PZ');

% Exclude border of the zones
BORD_COL = 6; BORD_INT = 4;
SC(NZ).ExRange = BORD_COL:SC(NZ).NoC-BORD_INT;
SC(MZ).ExRange = BORD_INT:SC(MZ).NoC-BORD_INT;
SC(PZ).ExRange = BORD_COL:SC(PZ).NoC-BORD_INT;

% All start points (A) of zone MZ
Distal_MZ = nan(length(SC(MZ).ExRange),3);
for c = SC(MZ).ExRange
    Distal_MZ(c-BORD_INT+1,:) = SC(MZ).P(c).xyz(SC(MZ).P(c).A,:);
end
% Most anterior point of all starting points (A) of zone MZ
[~, I_Distal_MZ_Xmax] = max(Distal_MZ(:,1));
EP.Intercondylar = Distal_MZ(I_Distal_MZ_Xmax,:);

% All start points (A) of zone NZ
ProximalPosterior_NZ = nan(length(SC(NZ).ExRange),3);
for c = SC(NZ).ExRange
    ProximalPosterior_NZ(c-BORD_COL+1,:) = SC(NZ).P(c).xyz(SC(NZ).P(c).A,:);
end
% Start point RC-NZ: Most proximal point of all start points (A) of zone NZ
[~, I_ProximalPosterior_NZ_Ymax] = max(ProximalPosterior_NZ(:,2));
EP.Medial = ProximalPosterior_NZ(I_ProximalPosterior_NZ_Ymax,:);

% All start points (A) of zone PZ
ProximalPosterior_PZ = nan(length(SC(PZ).ExRange),3);
for c = SC(PZ).ExRange
    ProximalPosterior_PZ(c-BORD_COL+1,:) = SC(PZ).P(c).xyz(SC(PZ).P(c).A,:);
end
% Start point RC-PZ: Most proximal point of all start points (A) of zone PZ
[~, I_ProximalPosterior_PZ_Ymax] = max(ProximalPosterior_PZ(:,2));
EP.Lateral = ProximalPosterior_PZ(I_ProximalPosterior_PZ_Ymax,:);

% All end points (B) of the zones NZ, MZ, PZ
ProximalAnterior = nan(length([SC.ExRange]),3);
countIdx = 1;
for s=1:LS
    for c = SC(s).ExRange
        ProximalAnterior(countIdx,:) = SC(s).P(c).xyz(SC(s).P(c).B,:);
        countIdx = countIdx+1;
    end; clear c
end
% End point RC-MZ: Most proximal point of all end points (B) of the zones NZ, MZ, PZ
[~, I_ProximalAnterior_Ymax] = max(ProximalAnterior(:,2));
EP.Anterior = ProximalAnterior(I_ProximalAnterior_Ymax,:);

%% Visualization
if visu == 1
    FigColor = [1 1 1];
    MonitorsPos = get(0,'MonitorPositions');
    Fig = figure('Units','pixels',...
        'Color',FigColor,'ToolBar','figure',...
        'WindowScrollWheelFcn',@M_CB_Zoom,'WindowButtonDownFcn',@M_CB_RotateWithLeftMouse,...
        'renderer','opengl');
    if     size(MonitorsPos,1) == 1
        set(Fig,'OuterPosition',MonitorsPos(1,:));
    elseif size(MonitorsPos,1) == 2
        set(Fig,'OuterPosition',MonitorsPos(2,:));
    end
    
    H_Axes = axes;
    axis on; xlabel('X [mm]'); ylabel('Y [mm]'); zlabel('Z [mm]');
    set(H_Axes,'Color',FigColor);
    daspect([1 1 1])
    cameratoolbar('SetCoordSys','none')
    hold on
    
    BoneProps.EdgeColor = 'none';
    BoneProps.FaceColor = [0.882, 0.831, 0.753];
    BoneProps.FaceAlpha = 0.7;
    patch('Faces',distalFemurUSP.faces, 'Vertices',distalFemurUSP.vertices, BoneProps);
    
    drawLine3d(PFEA, 'b')
    
    for i=1:size(PFEA_IS_Pts,1)
        scatter3(PISP{i,3}(1),PISP{i,3}(2),PISP{i,3}(3),'m','filled');
    end; clear i
    
    for s=1:LS
        for c=1:SC(s).NoC
            % Plot contour-part in 3D
            CP3D = SC(s).P(c).xyz;
            EPA = SC(s).P(c).A; EPB = SC(s).P(c).B;
            plot3(CP3D(EPA:EPB,1),CP3D(EPA:EPB,2),CP3D(EPA:EPB,3),'Color', SC(s).Color,'Linewidth',1);
        end; clear c
    end; clear s
    
    T = EP.Medial;  scatter3(T(1),T(2),T(3),'r','filled');
    text(T(1),T(2),T(3), 'Medial')
    T = EP.Intercondylar; scatter3(T(1),T(2),T(3),'r','filled');
    text(T(1),T(2),T(3), 'Intercondylar')
    T = EP.Lateral;  scatter3(T(1),T(2),T(3),'r','filled');
    text(T(1),T(2),T(3), 'Lateral')
    T = EP.Anterior;  scatter3(T(1),T(2),T(3),'r','filled');
    text(T(1),T(2),T(3), 'Anterior')
    
    hold off
    
    set(H_Axes,'CameraTarget',[0, 0, 0])
    CamPos = [-0.6499, 0.4339, 0.6240] * norm(get(H_Axes,'CameraPosition'));
    set(H_Axes,'CameraPosition',CamPos)
    set(H_Axes,'CameraUpVector',[0, 1, 0])
    
    light1 = light; light('Position', -1*(get(light1,'Position')));
    lighting phong
end
end


%% Subfunction: Mouse stuff
function M_CB_RotateWithLeftMouse(src,~)
if strcmp(get(src,'SelectionType'),'normal')
    cameratoolbar('SetMode','orbit')
end
end

function M_CB_Zoom(~,evnt)
if evnt.VerticalScrollCount > 0
    CVA_old = get(gca,'CameraViewAngle');
    CVA_new = CVA_old + 1;
    draw
elseif evnt.VerticalScrollCount < 0
    CVA_old = get(gca,'CameraViewAngle');
    CVA_new = CVA_old - 1;
    draw
end
    function draw
        set(gca,'CameraViewAngle',CVA_new)
        drawnow
    end
end