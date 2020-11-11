function [ExP, ExA, cpfh] = sagittalExPts_IntCond(Contour, sigmastart, sigmadelta, sigma, vis)
% A pattern-recognition algorithm for identifying the articulating surface
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
%       ExP        - integer: posterior extreme point P
%       ExA        - integer: anterior extreme point A
%       cpfh       - figure handle: empty if vis == 0
%
% AUTHOR: Maximilian C. M. Fischer
% COPYRIGHT (C) 2020 Maximilian C. M. Fischer
% LICENSE: EUPL v1.2
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

%% Find the posterior extreme point P (ExP)
% Which Gaussian of width (sigma) should be used to find the local maxima
if sigma > length(zcp)
    sigma = length(zcp);
end
[~, Local_Maxima_Indcs] = findpeaks(K{sigma});

% Find the inferior extreme point P (ExP)
ZCP_Candidates = zcp{sigma}(zcp{sigma}>IXMin & zcp{sigma}<IYMin);
zero_ExP = [];
if ~isempty(ZCP_Candidates)
    zero_ExP = ZCP_Candidates(1);
    lowerBound = zero_ExP;
else
    lowerBound = IXMin;
end
ExP_Candidates = Local_Maxima_Indcs(Local_Maxima_Indcs>lowerBound & Local_Maxima_Indcs<IYMin);
if ~isempty(ExP_Candidates)
    % Use the largest peak
    [~, ExP_CandMaxIdx] = max(K{sigma}(ExP_Candidates));
    ExP = ExP_Candidates(ExP_CandMaxIdx);
else
    % If no candiates are found
    ExP  = lowerBound + round((Local_Maxima_Indcs(knnsearch(Local_Maxima_Indcs, IYMin)) - lowerBound)/2);
end

%% Anterior extremities A
zero_ExA = [];
ExA = IXMax;
ZCP_A_Candidates = zcp{sigma}(zcp{sigma}>IYMin);
if ~isempty(ZCP_A_Candidates)
    zero_ExA = ZCP_A_Candidates(1);
    ExA_Candidates = Local_Maxima_Indcs(Local_Maxima_Indcs>IYMin & Local_Maxima_Indcs<zero_ExA);
    if ~isempty(ExA_Candidates)
        ExA  = ExA_Candidates(end);
    end
end

%% Visualization
cpfh = [];
if vis == 1 || vis == 2
    %% Plot: Contour
    cpfh = figure('Name','Intercondylar Contour', 'Color','w', 'WindowState','Maximized');
    axH = axes();
    plot(axH, Contour(:,1),Contour(:,2),'k-','LineWidth',2);
    hold(axH, 'on')
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
    
    % Visualize the counter-clockwise running direction: Arrow at YMax.
    quiver(axH, Contour(1,1),Contour(1,2),Contour(6,1)-Contour(1,1),Contour(6,2)-Contour(1,2),...
        'g','LineWidth',3,'AutoScale','off','MaxHeadSize',30);
    if ~isempty(zero_ExP)
        scatter(axH, Contour(zero_ExP,1),Contour(zero_ExP,2), 'filled');
        text(axH, Contour(zero_ExP,1),Contour(zero_ExP,2), 'Zero crossing point P',...
            'VerticalAlignment','bottom','HorizontalAlignment','right');
    end
    if ~isempty(zero_ExA)
        scatter(axH, Contour(zero_ExA,1),Contour(zero_ExA,2), 'filled');
        text(axH, Contour(zero_ExA,1),Contour(zero_ExA,2), 'Zero crossing point A',...
            'VerticalAlignment','bottom','HorizontalAlignment','left');
    end
    
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
    quiver(axH, Contour(:,1),Contour(:,2),Normals(:,1),Normals(:,2),...
        'color','k','ShowArrowHead','off','AutoScale','off','LineStyle','--')
    % Connect the tips of the normals
    NormalEnds = Normals + Contour;
    plot(axH, NormalEnds(:,1),NormalEnds(:,2),'k-.','LineWidth',1.5)
    
    %% Plot extreme points of the articulating surface
    % Plot the posterior extreme point P (ExP)
    scatter(axH, Contour(ExP,1),Contour(ExP,2), 'filled');
    text(axH, Contour(ExP,1),Contour(ExP,2), ['P for \sigma = ' num2str(sigma)],...
        'VerticalAlignment','bottom');
    
    % Plot the anterior medial extreme point A (ExA)
    scatter(axH, Contour(ExA,1),Contour(ExA,2), 'filled');
    if ExA~=IXMax
        text(axH, Contour(ExA(1),1),Contour(ExA(1),2), ['intercondylar A for \sigma = ' num2str(sigma)],...
            'HorizontalAlignment','right');
    end
    
    scatter(axH, Contour(IXMin,1),Contour(IXMin,2), 'k', 'filled');
    text(axH, Contour(IXMin,1),Contour(IXMin,2), 'X_{Min}','HorizontalAlignment','left');
    scatter(axH, Contour(IXMax,1),Contour(IXMax,2), 'k', 'filled');
    text(axH, Contour(IXMax,1),Contour(IXMax,2), 'X_{Max}','HorizontalAlignment','right');
    scatter(axH, Contour(IYMin,1),Contour(IYMin,2), 'k', 'filled');
    text(axH, Contour(IYMin,1),Contour(IYMin,2), 'Y_{Min}','VerticalAlignment','bottom');
    scatter(axH, Contour(IYMax,1),Contour(IYMax,2), 'k', 'filled');
    text(axH, Contour(IYMax,1),Contour(IYMax,2), 'Y_{Max}','VerticalAlignment','top');
    
    axis(axH,'equal');
    title(axH,'Intercondylar Contour')
    
    if vis == 2
        %% Plot: Curvature Scale Space (CSS) Image
        figure('name', 'Curvature Scale Space (CSS) Image'),
        for i=1:length(zcp)
            plot(zcp{i},S{i},'k.');
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