function [ExA, ExB, cpfh] = sagittalExPts_LatCond(Contour, sigmastart, sigmadelta, sigma, vis)
%
%    - A Pattern-Recognition Algorithm for Identifying the Articulating Surface
%   
%   HELP:                                       
%
%   INPUT:
%       Contour    - nx2 double: X- & Y-coordinates of the contour
%                      	Requirements: - Sorting: counter-clockwise, 
%                                     - Start point: Max Y-value
%       sigmastart - Starting sigma value (see BOMultiScaleCurvature2D_adapted)
%       sigmadelta - Delta value of sigma (see BOMultiScaleCurvature2D_adapted)
%       sigma      - Sigma, that is used to detect the local maxima
%       vis        - visualization options:
%                               - 0: Plot nothing
%                               - 1: Plot contour
%
%   OUTPUT:
%       ExA      - integer: extreme point A
%       ExB      - integer: extreme point B
%       cpfh     - figure handle: empty if vis == 0 
%
%   USAGE:
% 
%   AUTHOR: MCMF
%
%   VERSION:
%       

%% Calculations
% Calculate the multi-scale curvature & the curvature scale-space image
[K,S,~,~,Xsm,Ysm,zcp] = BOMultiScaleCurvature2D_adapted(Contour,sigmastart,sigmadelta);

% Boundaries of the contour
[XMin, IXMin] = min(Contour(:,1));
[XMax, IXMax] = max(Contour(:,1));
[YMin, IYMin] = min(Contour(:,2));
[YMax, IYMax] = max(Contour(:,2));
% IYMax should always be 1, because the contour should start there
if IYMax ~=1
    warning('Contour should start at the max. Y value (YMax): Algorithm1 won''t work!')
end
if IXMax == 1
    XMax = Contour(end,1);
    IXMax = size(Contour,1);
end
    
%% Find the posterior extreme point A (ExA)
% To find the curvature zerocrossing point of ExA:
% Find the sigma that meets the following condition:
% Only 2 zero crossings lie between IYMax and IXMin of the contour
sigma_ExA = 1;
for z=1:length(zcp(1:end-1))
    if sum(zcp{z}>IYMax & zcp{z}<IXMin) > 2
        sigma_ExA = z+1;
    end
end
zero_ExA = max(zcp{sigma_ExA}(zcp{sigma_ExA}>IYMax & zcp{sigma_ExA}<IXMin));
% If no zero crossing is found set IYMax as crossing point
if isempty(zero_ExA)
    zero_ExA = IYMax;
end

% Which Gaussian of width (sigma) should be used to find the local maxima
if sigma == 0
    sigma = sigma_ExA;
elseif sigma > length(zcp)
    sigma = length(zcp);
end
[~, Local_Maxima_Indcs] = findpeaks(K{sigma});

% Find the posterior extreme point A (ExA) defined as:
% - local maximum curvature point inferior to the curvature zerocrossing point (zero_pExA)
ExA_Candidates = Local_Maxima_Indcs(Local_Maxima_Indcs>zero_ExA & Local_Maxima_Indcs<IXMin);
% If no candiates are found set zero_ExA as pExA
if ~isempty(ExA_Candidates)
    ExA  = ExA_Candidates(1);
else
    ExA = zero_ExA;
    warning('ExA = zero_ExA!')
end

%% Find the anterior extreme point B
% - If there are zero crossing points between the most anterior point and 
%   the most inferior point of the contour, use the maximum of kappa
%   between the most inferior point and the first zero crossing point after
%   the most inferior point
if sum(zcp{sigma}>IYMin & zcp{sigma}<IXMax) > 1
    IYMin_IXMax_zcp = zcp{sigma}(zcp{sigma}>IYMin & zcp{sigma}<IXMax);
%     [~, alExB] = min(K{sigma}(IYMin_IXMax_zcp(1):IYMin_IXMax_zcp(2)));
%     alExB = IYMin_IXMax_zcp(1)-1+alExB;
    [~, ExB] = max(K{sigma}(IYMin:IYMin_IXMax_zcp(1)));
    ExB = IYMin-1+ExB;
else
% - else use the local maximum curvature point closest to the most anterior point of the contour.
ExB = Local_Maxima_Indcs(knnsearch(Local_Maxima_Indcs, IXMax));
end

%% Visualization
cpfh = [];
if vis == 1 || vis == 2   
    %% Plot: Contour
    cpfh = figure('name', 'Contour');
    title('Contour');
    plot(Contour(:,1),Contour(:,2),'k-','LineWidth',2);
    hold on;
    for i=1:length(zcp)
        Xs = Xsm{i};
        Ys = Ysm{i};
        zc = zcp{i};
        sch(1) = plot(Xs,Ys,'r-');
        hold on;
        sch(2) = plot(Xs(zc),Ys(zc),'go');
        pause(0);
        delete(sch);
    end
    
    % Visualization of the Running direction: Arrow -> at YMax
    % Should be counter-clockwise
    quiver(Contour(1,1),Contour(1,2),Contour(6,1)-Contour(1,1),Contour(6,2)-Contour(1,2),...
        'g','LineWidth',3,'AutoScale','off','MaxHeadSize',30);
    scatter(Contour(zero_ExA,1),Contour(zero_ExA,2), 'filled');
    text(Contour(zero_ExA,1),Contour(zero_ExA,2), 'Zero crossing point',...
        'VerticalAlignment','bottom');
    
    %% Normals of the contour
    % Get the normals (already normed)
    [~, Normals, ~, ~, ~] = frenet(Contour(:,1),Contour(:,2));
    Normals = Normals(:,1:2);
    % Unify normals: All normals have to point outside the contour
    % If normal points inside resp. the tip of the normal is inside the
    % contour, change the direction
    Indcs_Inside_Normals = isPointInPolygon( (Contour+Normals) , Contour);
    % Change the direction of all normals pointing inside the contour
    Normals(Indcs_Inside_Normals,1:2)=Normals(Indcs_Inside_Normals,1:2) * -1;
    % A scaling factor for the length of the normals
    ScalingFactor = max(abs([XMin,XMax,YMin,YMax]));
    % Multiply the normals with kappa (K), to visualize K on the contour
    Normals = repmat(K{sigma},1,2).*Normals*ScalingFactor;
    quiver(Contour(:,1),Contour(:,2),Normals(:,1),Normals(:,2),...
        'color','k','ShowArrowHead','off','AutoScale','off','LineStyle','--')
    % Connect the tips of the normals
    NormalEnds = Normals + Contour;
    plot(NormalEnds(:,1),NormalEnds(:,2),'k-.','LineWidth',1.5)
    
    %% Plot extreme points of the articulating surface
    % Plot the posterior extreme point A (pExA)
    scatter(Contour(ExA,1),Contour(ExA,2), 'filled');
    text(Contour(ExA,1),Contour(ExA,2), ['A for \sigma = ' num2str(sigma)],...
        'VerticalAlignment','top');
    % Plot the anterior lateral extreme point B (alExB)
    scatter(Contour(ExB,1),Contour(ExB,2), 'filled');
    text(Contour(ExB,1),Contour(ExB,2), ['B for \sigma = ' num2str(sigma)],...
        'VerticalAlignment','bottom','HorizontalAlignment','right');
    
    scatter(Contour(IXMin,1),Contour(IXMin,2), 'k', 'filled');
    text(Contour(IXMin,1),Contour(IXMin,2), 'X_{Min}','HorizontalAlignment','left');
    scatter(Contour(IXMax,1),Contour(IXMax,2), 'k', 'filled');
    text(Contour(IXMax,1),Contour(IXMax,2), 'X_{Max}','HorizontalAlignment','right');
    scatter(Contour(IYMin,1),Contour(IYMin,2), 'k', 'filled');
    text(Contour(IYMin,1),Contour(IYMin,2), 'Y_{Min}','VerticalAlignment','bottom');
    scatter(Contour(IYMax,1),Contour(IYMax,2), 'k', 'filled');
    text(Contour(IYMax,1),Contour(IYMax,2), 'Y_{Max}','VerticalAlignment','top');
    
    axis equal;
    title('Lateral')
    
    if vis == 2 % <- Set this to 1 if Plots are needed
        %% Plot: Curvature Scale Space (CSS) Image
        figure('name', 'Curvature Scale Space (CSS) Image'),
        for i=1:length(zcp)
            plot(zcp{i},S{i},'k.');%,'MarkerFaceColor',C(i,:),'MarkerEdgeColor',C(i,:));
            hold on
        end
        hold off
        xlim([1 length(K{1})]);
        title('Curvature Scale Space (CSS) Image');
        
        if length(K) > 1
            %% Plot: Curvature kappa
            figure
            surf(cell2mat(K))
            xlabel('\sigma');ylabel('u');zlabel('\kappa(u,\sigma)');
            title('Curvature of \Gamma_\sigma: \kappa(u,\sigma)');
        end
    end
end

end