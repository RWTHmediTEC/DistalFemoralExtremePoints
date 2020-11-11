function EP = distalFemoralExtremePoints(distalFemurUSP, side, PFEA, varargin)
% DISTALFEMORALEXTREMEPOINTS detects landmarks of the distal femur based on
%   sagittal cutting contours.
%
%   EP = distalFemurExtremePoints(distalFemurUSP, side, PFEA) returns the
%   struct EP containing multiple landmarks of the distal femur, such as
%   the intercondylar notch (ICN), the proximoposterior point of the 
%   medial and lateral condyle (MPPC, LPPC) [beta version!]. 
% 
% INPUT
%   REQUIRED:
%     distalFemurUSP - struct: A clean mesh of the distal femur defined by  
%           the fields vertices (double [Nx3]) and faces (integer [Mx3]) 
%           transformed into the USP coordinate system (see USP.m)
%     side - char: 'Left' or 'Right' femur. Should start with 'L' or 'R'.
%     PFEA - double [1x6]: A line fitted through the posterior foci of the
%            ellipses with minimum dispersion (see USP.m).
% 
%   OPTIONAL:
%     'Visualization' - Logical: Figure output. Default is false.
% 
% OUTPUT
%     EP - struct containing fields with the name of the the landmarks with 
%           the xyz-coordinates [1x3] in the USP coordinate system:
%           ICN: intercondylar notch
%           MPPC: medial proximoposterior condyle [beta version!]
%           LPPC: lateral proximoposterior condyle [beta version!]
% 
% EXAMPLE:
%     Run the file DistalFemurExtremePoints_Example.m
%
% REFERENCE:
%     Based on: 2010 - Li - Automating Analyses of the Distal Femur 
%     Articular Geometry Basedon Three-Dimensional Surface Data
% 
% TODO:
%   - Complete documentation
%   - Maybe add input of landmarks for additonal sanity checks
% 
% AUTHOR: Maximilian C. M. Fischer
% 	mediTEC - Chair of Medical Engineering, RWTH Aachen University
% VERSION: 1.0.2
% DATE: 2020-09-06
% COPYRIGHT (C) 2020 Maximilian C. M. Fischer
% LICENSE: EUPL v1.2

addpath(genpath([fileparts([mfilename('fullpath'), '.m']) '\src']));

%% Parse inputs
p = inputParser;
logParValidFunc=@(x) (islogical(x) || isequal(x,1) || isequal(x,0));
addRequired(p,'distalFemurUSP',@(x) isstruct(x) && isfield(x, 'vertices') && isfield(x,'faces'))
addRequired(p,'side',@(x) any(validatestring(upper(x(1)),{'R','L'})));
addParameter(p,'visualization',false,logParValidFunc);
addParameter(p,'debugVisualization',false,logParValidFunc);
parse(p, distalFemurUSP, side, varargin{:});
side = upper(p.Results.side(1));
visu = logical(p.Results.visualization);
debugVisu = logical(p.Results.debugVisualization);

%% Settings
% Sigma: For explanation see the function: BOMultiScaleCurvature2D_adapted.m
sigmaStart = 1;
sigmaDelta = 1;
sigma = 4;

%% Find boundary points of the condyles
% Check that PFEA points in positive z direction
if PFEA(6)<0; PFEA(4:6) = -PFEA(4:6); end
% Intersection points between the mesh an the posterior foci elliptical axis
[PFEA_IS_Pts, PFEA_Pos] = intersectLineMesh3d(PFEA, distalFemurUSP.vertices, distalFemurUSP.faces);
% Sort inters. points along the PFEA
[PFEA_Pos, sortIdx] = sort(PFEA_Pos);
PFEA_IS_Pts = PFEA_IS_Pts(sortIdx,:);

% Should always be four inters. points
PFEA_error = 'PFEA intersection points with the distal femur mesh ~= 4.';
if size(PFEA_IS_Pts,1) > 4
    % Otherwise try to merge close points
    PFEA_THRESHOLD = 5; % Threshold in [mm] to merge points
    % Find close points
    diffIdx = find(diff(PFEA_Pos) < PFEA_THRESHOLD);
    % Delete them based on the position
    for d=1:length(diffIdx)
        if PFEA_IS_Pts(diffIdx(d),3)<0
            PFEA_IS_Pts(diffIdx(d),:) = nan;
        else
            PFEA_IS_Pts(diffIdx(d)+1,:) = nan;
        end
    end
    PFEA_IS_Pts(any(isnan(PFEA_IS_Pts),2),:) = [];
    % Check again
    if size(PFEA_IS_Pts,1) ~= 4
        error(PFEA_error)
    end
elseif size(PFEA_IS_Pts,1) < 4
    error(PFEA_error)
end

if visu
    patchProps.FaceAlpha = 1;
    [~, axH, figH] = visualizeMeshes(distalFemurUSP, patchProps); %#ok<ASGLU>
    drawPoint3d(axH, PFEA_IS_Pts, 'MarkerFaceColor','r', 'MarkerEdgeColor','r','MarkerSize',10)
end

% The intersection points devide the distal femur into 3 parts in Z direction:
% Negative Z = NZ
% Middle     = MZ
% Positive Z = PZ
switch side
    case 'R'
        % RIGHT knee: medial - neg. Z values (NZ), lateral - pos. Z values (PZ)
        SC(1).Zone = 'NZ';
        SC(2).Zone = 'MZ';
        SC(3).Zone = 'PZ';
    case 'L'
        %  LEFT knee: medial - pos. Z values (PZ), lateral - neg. Z values (NZ)
        SC(1).Zone = 'PZ';
        SC(2).Zone = 'MZ';
        SC(3).Zone = 'NZ';
end

% Start and end point of the zones, rounded to an integer value in Z direction
BORD_COL = 6; BORD_INT = 4;
NZS = PFEA_IS_Pts(1,:); NZS(3) = ceil(NZS(3))+BORD_COL;
NZE = PFEA_IS_Pts(2,:); NZE(3) = floor(NZE(3));
MZS = PFEA_IS_Pts(2,:); MZS(3) = ceil(MZS(3));
MZE = PFEA_IS_Pts(3,:); MZE(3) = floor(MZE(3));
PZS = PFEA_IS_Pts(3,:); PZS(3) = ceil(PZS(3));
PZE = PFEA_IS_Pts(4,:); PZE(3) = floor(PZE(3))-BORD_COL;

% The number of sagittal cuts per zone
SC(1).NoC = abs(NZS(3)-NZE(3))+1;
SC(2).NoC = abs(MZS(3)-MZE(3))+1;
SC(3).NoC = abs(PZS(3)-PZE(3))+1;

SC(1).Origin = NZS;
SC(2).Origin = MZS;
SC(3).Origin = PZS;

LS = size(SC,2);
ZoneColors = [.3 .3 1; 0 0 1; .6 .6 1];

%% Sagittal Cuts (SC)
ZVector = [0, 0, 1];

Contours = cell(1,LS);
for s=1:LS
    SC(s).Color = ZoneColors(s,:);
    % Create cutting plane origins
    SC(s).PlaneOrigins = repmat(SC(s).Origin, SC(s).NoC, 1) + [zeros(SC(s).NoC,2), (0:SC(s).NoC-1)'];
    % Create SC(s).NoC Saggital Contour Profiles (SC(s).P)
    Contours{s} = IntersectMeshPlaneParfor(distalFemurUSP, SC(s).PlaneOrigins, ZVector);
    for c=1:SC(s).NoC
        % If there is more than one closed contour after the cut, use the longest one
        [~, IobC] = max(cellfun(@length, Contours{s}{c}));
        SC(s).P(c).xyz = Contours{s}{c}{IobC}';
        SC(s).P(c).xyz = unique(SC(s).P(c).xyz,'rows','stable');
        % Set Start of the contour to the maximum Y value
        [~, IYMax] = max(SC(s).P(c).xyz(:,2));
        if IYMax ~= 1
            SC(s).P(c).xyz = SC(s).P(c).xyz([IYMax:length(SC(s).P(c).xyz),1:IYMax-1],:);
        end
        % Close contour: Copy start value to the end
        if ~isequal(SC(s).P(c).xyz(1,:),SC(s).P(c).xyz(end,:))
            SC(s).P(c).xyz(end+1,:) = SC(s).P(c).xyz(1,:);
        end
        % If the contour is sorted clockwise
        if sign(varea(SC(s).P(c).xyz')) == -1
            % Sort the contour counter-clockwise
            SC(s).P(c).xyz = flipud(SC(s).P(c).xyz);
        end
    end
end

%% sagittalExPts = Sagittal Extreme Points
% A pattern-recognition algorithm for identifying the articulating surface
for s=1:LS
    for c=1:SC(s).NoC
        % Only the X & Y values are needed because the countours are parallel to the X-Y plane.
        Contour = SC(s).P(c).xyz(:,1:2);
        % Get the anterior and posterior extreme points of the articulating surface
        switch SC(s).Zone
            case 'NZ'
                [SC(s).P(c).A, SC(s).P(c).B, SC(s).P(c).H] = ...
                    sagittalExPts_MedCond(Contour, sigmaStart, sigmaDelta, sigma, debugVisu);
            case 'MZ'
                [SC(s).P(c).A, SC(s).P(c).B, SC(s).P(c).H] = ...
                    sagittalExPts_IntCond(Contour, sigmaStart, sigmaDelta, sigma, debugVisu);
            case 'PZ'
                [SC(s).P(c).A, SC(s).P(c).B, SC(s).P(c).H] = ...
                    sagittalExPts_LatCond(Contour, sigmaStart, sigmaDelta, sigma, debugVisu);
        end
    end
end


%% Find global extreme points of the sagittal contours (SCs)
NZ = strcmp({SC(:).Zone}, 'NZ');
MZ = strcmp({SC(:).Zone}, 'MZ');
PZ = strcmp({SC(:).Zone}, 'PZ');

% Exclude border of the zones
if SC(MZ).NoC > 2*BORD_INT
    SC(MZ).ExRange = BORD_INT:SC(MZ).NoC-BORD_INT;
else
    SC(MZ).ExRange = 1:SC(MZ).NoC;
end
switch side
    case 'R'
        SC(NZ).ExRange = 1:SC(NZ).NoC-BORD_INT;
        SC(PZ).ExRange = BORD_INT:SC(PZ).NoC;
    case 'L'
        SC(NZ).ExRange = BORD_INT:SC(NZ).NoC;
        SC(PZ).ExRange = 1:SC(PZ).NoC-BORD_INT;
end

%% Intercondylar Notch (ICN)
% Take start points (A) of zone MZ specified by ExRange
A_MZ = nan(SC(MZ).NoC,3);
for c = SC(MZ).ExRange
    A_MZ(c,:) = SC(MZ).P(c).xyz(SC(MZ).P(c).A,:);
end
A_MZ = A_MZ(~any(isoutlier(A_MZ),2),:);
% Median of the points A
A_MZ_median = median(A_MZ,'omitnan');

% Clip region of the ICN
ICNmesh = cutMeshByPlane(distalFemurUSP, [MZS, 0 1 0, 0 0 1]);
ANTERIOR_CUT = 5; %[mm]
ICNmesh = cutMeshByPlane(ICNmesh, [A_MZ_median(1)+ANTERIOR_CUT A_MZ_median(2:3), 0 1 0, 0 0 -1]);
ICNmesh = cutMeshByPlane(ICNmesh, [MZS, 1 0 0, 0  1 0]);
ICNmesh = cutMeshByPlane(ICNmesh, [MZE, 1 0 0, 0 -1 0]);
ICNmesh = cutMeshByPlane(ICNmesh, [MZE, 1 0 0, 0 0 1]);

% Get contour of the ICN mesh in the y-z plane
ICNcontour = meshSilhouette(ICNmesh, [0 0 0 0 1 0 0 0 1],'visu', 0);

% Shrink contour on both sides and on the top
ICNcontour(isBelowPlane(ICNcontour,[MZS(1:2) MZS(3)+2, 1 0 0, 0  1 0]),:) = [];
ICNcontour(isBelowPlane(ICNcontour,[MZE(1:2) MZE(3)-2, 1 0 0, 0 -1 0]),:) = [];
ICNcontour(isBelowPlane(ICNcontour,[MZE(1) MZE(2)-3 MZE(3), 1 0 0, 0 0 1]),:) = [];
ICNcontour(any(isnan(ICNcontour),2),:) = [];

% Remove outliers of the distal contour and take the longest part
if size(ICNcontour,1) >= 3
    ICNcontour = angleSort3d(ICNcontour, [0 0 0]);
end
[~, ICNcontourSortIdx] = sort(ICNcontour(:,3));
ICNcontour = ICNcontour(ICNcontourSortIdx,:);
[~, ICNcontourUniqueIdx] = uniquetol(ICNcontour(:,3), 3e-2);
ICNcontour = ICNcontour(ICNcontourUniqueIdx,:);
ICNcontour(isoutlier(ICNcontour(:,2),'movmedian',5),:) = nan;
ICNcontour = splitPolygons(ICNcontour);
[~, ICNcontourMax] = max(cellfun(@(x) size(x,1),ICNcontour));
ICNcontour = ICNcontour{ICNcontourMax};

% Find the maximum of the contour in y direction
[~, ICNconYmaxIdx] = max(ICNcontour(:,2));

% Construct the final ICN point
ICNpoint = [A_MZ_median(1) ICNcontour(ICNconYmaxIdx,2:3)];
[~, ICN_Idx] = pdist2(A_MZ,ICNpoint,'euclidean','Smallest',1);
[~, ICN_Idx] = pdist2(distalFemurUSP.vertices, A_MZ(ICN_Idx,:),'euclidean','Smallest',1);
EP.Intercondylar = distalFemurUSP.vertices(ICN_Idx,:);

%% Take start points (A) of zone NZ specified by ExRange
ProximalPosterior_NZ = nan(SC(NZ).NoC,3);
for c = SC(NZ).ExRange
    ProximalPosterior_NZ(c,:) = SC(NZ).P(c).xyz(SC(NZ).P(c).A,:);
end
ProximalPosterior_NZ(any(isnan(ProximalPosterior_NZ),2),:)=[];
% Median of the start points
[~, PPNZ_Idx] = pdist2(distalFemurUSP.vertices,median(ProximalPosterior_NZ),'euclidean','Smallest',1);
EP.Medial = distalFemurUSP.vertices(PPNZ_Idx,:);

%% Take start points (A) of zone PZ specified by ExRange
ProximalPosterior_PZ = nan(SC(PZ).NoC,3);
for c = SC(PZ).ExRange
    ProximalPosterior_PZ(c,:) = SC(PZ).P(c).xyz(SC(PZ).P(c).A,:);
end
ProximalPosterior_PZ(any(isnan(ProximalPosterior_PZ),2),:)=[];
% Median of the start points
[~, PPPZ_Idx] = pdist2(distalFemurUSP.vertices,median(ProximalPosterior_PZ),'euclidean','Smallest',1);
EP.Lateral = distalFemurUSP.vertices(PPPZ_Idx,:);

%% All end points (B) of the zones NZ, MZ, PZ
% ProximalAnterior = nan(sum([SC.NoC]),3);
% countIdx = 1;
% for s=1:LS
%     for c = 1:SC(s).NoC
%         ProximalAnterior(countIdx,:) = SC(s).P(c).xyz(SC(s).P(c).B,:);
%         countIdx = countIdx+1;
%     end; clear c
% end
% ProximalAnterior(any(isoutlier(ProximalAnterior,'movmedian',10),2),:) = nan;
% % End point RC-MZ: Most proximal point of all end points (B) of the zones NZ, MZ, PZ
% [~, I_ProximalAnterior_Ymax] = max(ProximalAnterior(:,2));
% EP.Anterior = ProximalAnterior(I_ProximalAnterior_Ymax,:);

%% Visualization
if visu
    drawEdge3d(axH, ...
        [PFEA_IS_Pts(1,1:2) PFEA_IS_Pts(1,3)-5],...
        [PFEA_IS_Pts(4,1:2) PFEA_IS_Pts(4,3)+5],...
        'Color', 'r','Linewidth',2)
    
    for s=1:LS
        for c=1:SC(s).NoC
            % Plot contour-part in 3D
            CP3D = SC(s).P(c).xyz;
            EPA = SC(s).P(c).A; EPB = SC(s).P(c).B;
            plot3(axH, CP3D(EPA:EPB,1),CP3D(EPA:EPB,2),CP3D(EPA:EPB,3),...
                'Color', SC(s).Color,'Linewidth',1,'LineStyle','--');
        end
        for c=SC(s).ExRange
            % Plot contour-part in 3D
            CP3D = SC(s).P(c).xyz;
            EPA = SC(s).P(c).A; EPB = SC(s).P(c).B;
            plot3(axH, CP3D(EPA:EPB,1),CP3D(EPA:EPB,2),CP3D(EPA:EPB,3),...
                'Color', SC(s).Color,'Linewidth',2);
        end
    end
    
    drawPolyline3d(axH, A_MZ,'Color','b','LineWidth',1,'LineStyle','--')
    drawPolyline3d(axH, A_MZ,'Color','b','LineWidth',3)
    
    pointProps.MarkerFaceColor = 'k';
    pointProps.MarkerEdgeColor = 'k';
    pointProps.Marker = 'o';
    pointProps.MarkerSize = 12;
    
    drawPoint3d(axH, EP.Medial,pointProps)
    drawLabels3d(axH, EP.Medial, 'Medial')
    drawPoint3d(axH, EP.Intercondylar,pointProps)
    drawLabels3d(axH, EP.Intercondylar, 'ICN',...
        'FontSize',14,'VerticalAlignment','Top')
    drawPoint3d(axH, EP.Lateral,pointProps)
    drawLabels3d(axH, EP.Lateral, 'Lateral')
    % drawPoint3d(axH, EP.Anterior,pointProps)
    % drawLabels3d(axH, EP.Anterior, 'Anterior')
    
    drawPoint3d(axH, ICNpoint, pointProps, 'MarkerFaceColor','r', 'MarkerEdgeColor','r')
    drawLabels3d(axH, ICNpoint, 'ICN_{temp}')
    drawPolyline3d(axH, ICNcontour,'Color','g','LineWidth',3)
    patchProps.EdgeColor = 'none';
    patchProps.FaceColor = 'g';
    patchProps.FaceAlpha = 0.5;
    patchProps.FaceLighting = 'gouraud';
    patch(axH, ICNmesh, patchProps)
    
    anatomicalViewButtons(axH, 'ASR')
    
    % % For publication
    % axis(axH, 'off')
    % camTar = [0 0 0];
    % camNorm = [-1 0 0];
    % set(axH, 'CameraTarget',camTar);
    % set(axH, 'CameraPosition',camTar + camNorm*700);
    % set(axH, 'CameraUpVector',normalizeVector3d([0 1 0]));
    % set(axH, 'CameraViewAngle',10)
    % set(figH, 'GraphicsSmoothing','off')
    % export_fig('Figure7', '-tif', '-r300')
end

end
