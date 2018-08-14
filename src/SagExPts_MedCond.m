function [ExA, ExB, cpfh] = SagExPts_MedCond(Contour, sigmastart, sigmadelta, sigma, vis)

%    - A Pattern-Recognition Algorithm for Identifying the Articulating Surface
%   
%   REFERENCE:
%       Li et al. - Automating Analyses of the Distal Femur Articular 
%       Geometry Basedon Three-Dimensional Surface Data
%       Annals of Biomedical Engineering, Vol. 38, No. 9, September 2010 
%       pp. 2928–2936                                     
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
%       ExA      - integer: posterior extremity A
%       ExB      - integer: anterior extremity B
%       cpfh     - figure handle: empty if vis == 0 
%
%   USAGE:
% 
%   AUTHOR: MCMF
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
    warning('Contour should start at the max. Y value (YMax)!')
end

%% Find the posterior extremity ExA
% To find the curvature zerocrossing point of ExA:
% Find the sigma that meets the following condition:
% Only 2 zero crossings lie between IYMax and IXMin of the contour
sigma_ExP = 1;
for z=1:length(zcp(1:end-1))
    if sum(zcp{z}>IYMax & zcp{z}<IXMin) > 2
        sigma_ExP = z+1;
    end
end
zero_ExP = max(zcp{sigma_ExP}(zcp{sigma_ExP}>IYMax & zcp{sigma_ExP}<IXMin));
% If no zero crossing is found set IYMax as crossing point
if isempty(zero_ExP)
    zero_ExP = IYMax;
end

% Which Gaussian of width (sigma) should be used to find the local maxima
if sigma == 0
    sigma = sigma_ExP;
elseif sigma > length(zcp)
    sigma = length(zcp);
end
[~, Local_Maxima_Indcs] = findpeaks(K{sigma});

% Find the posterior extremity P (ExA) defined as:
% - local maximum curvature point inferior to the curvature zerocrossing point (zero_ExP)
ExA_Candidates = Local_Maxima_Indcs(Local_Maxima_Indcs>zero_ExP & Local_Maxima_Indcs<IXMin);
% If no candiates are found set IYMax as ExA
if isempty(ExA_Candidates)
    ExA = IYMax;
    warning('ExA = IYMax!')
else
    ExA  = ExA_Candidates(1);
end

%% Find the medial anterior extremity ExB
% Find ExB defined as:
% - the local maximum curvature point with the largest curvature value in
%   the region between the most anterior point and the most inferior point of the contour;
[~, ExB] = max(K{sigma}(IYMin:IXMax));
ExB = IYMin-1+ExB;

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
    scatter(Contour(zero_ExP,1),Contour(zero_ExP,2), 'filled');
    text(Contour(zero_ExP,1),Contour(zero_ExP,2), 'Zero crossing point',...
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
    
    %% Plot extremity points of the articulating surface
    % Plot the posterior extremity P (ExP)
    scatter(Contour(ExA,1),Contour(ExA,2), 'filled');
    text(Contour(ExA,1),Contour(ExA,2), ['P for \sigma = ' num2str(sigma)],...
        'VerticalAlignment','top');
    
    % Plot the anterior medial extremity A (ExA)
    scatter(Contour(ExB,1),Contour(ExB,2), 'filled');
    text(Contour(ExB(1),1),Contour(ExB(1),2), ['medial A for \sigma = ' num2str(sigma)],...
        'HorizontalAlignment','right');
    
    scatter(Contour(IXMin,1),Contour(IXMin,2), 'k', 'filled');
    text(Contour(IXMin,1),Contour(IXMin,2), 'X_{Min}','HorizontalAlignment','left');
    scatter(Contour(IXMax,1),Contour(IXMax,2), 'k', 'filled');
    text(Contour(IXMax,1),Contour(IXMax,2), 'X_{Max}','HorizontalAlignment','right');
    scatter(Contour(IYMin,1),Contour(IYMin,2), 'k', 'filled');
    text(Contour(IYMin,1),Contour(IYMin,2), 'Y_{Min}','VerticalAlignment','bottom');
    scatter(Contour(IYMax,1),Contour(IYMax,2), 'k', 'filled');
    text(Contour(IYMax,1),Contour(IYMax,2), 'Y_{Max}','VerticalAlignment','top');
    
    axis equal;
    title('Medial')
    
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