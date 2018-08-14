function [K,S,X,Y,Xsm,Ysm,zcp] = BOMultiScaleCurvature2D_adapted(im,sigmastart,sigmadelta)
%% ADAPTED - NOT THE ORIGINAL FILE
% - Boundary section was changed, because the input is a already a contour instead of an image

%%  BOMultiScaleCurvature2D - calculates multi-scale curvature 
%                           and curvature scale-space image
%   
%   REFERENCE:
%       Farzin Mokhtarian and Alan Mackworth,
%       Scale-based description and recognition of planar curves and 
%       two-dimensional shapes,
%       IEEE Transactions on Pattern Analysis and Machine Intelligence,
%       1, 34-43, 1986.
%   
%   HELP:                                       
%
%   INPUT:
%       im          - binary image
%       sigmastart  - starting sigma value 
%       sigmadelta  - delta value of sigma
%
%   OUTPUT:
%       K           - curvatures calculated for each sigma
%       S           - sigmas
%       X           - x coor of contour
%       Y           - y coor of contour
%       Xsm         - x coor of contour convolved with Gaussian
%       Ysm         - x coor of contour convolved with Gaussian
%       zcp         - indexes of zero-crossing points
%
%   USAGE:
%       sigmastart = 1; sigmadelta = 0.1; p = 0.1;
%       [K,S,X,Y,Xsm,Ysm,zcp] = 
%                   BOMultiScaleCurvature2D(im,sigmastart,sigmadelta);
%       BOMultiScaleCurvaturePlot2D(K,S,X,Y,Xsm,Ysm,zcp,p);
% 
%   AUTHOR:
%       Boguslaw Obara, http://boguslawobara.net/
%
%   VERSION:
%       0.1 - 30/06/2009 First implementation

%% Boundary
% [row, col] = find(im,1);
% contour = bwtraceboundary(im, [row, col], 'N');
% X = contour(:,2); Y = contour(:,1);
X = im(:,1); Y = im(:,2);
%% WHILE
flag = 1; sigma = sigmastart; i = 1; zcp = [];
while flag == 1
    flag = 0;
%% Define Gauss Kernel
    % full width at half maximum (FWHM) of the peak     
    r = ceil(2*sqrt(2*log(2))*sigma); 
    r = 2*r;
    % Make the 1-D Gaussian 
    x = (-r:r)';
    G = (1/sqrt(2*pi*sigma.^2)).*exp(-x.^2./(2*sigma.^2));  
    G = G ./ sum(G);
%% Cutoff Points
    cop = r;
    Xn = []; Yn = [];
    Xn= [Xn; X(length(X)+1-cop:length(X))];
    Xn = [Xn; X];
    Xn= [Xn; X(1:cop)];
    Yn= [Yn; Y(length(Y)+1-cop:length(Y))];
    Yn = [Yn; Y];
    Yn= [Yn; Y(1:cop)];    
%% Calculate Derivatives of Gaussian
    % First Derivatives
    Gp = gradient(G);
    % Second Derivatives
    Gpp = gradient(Gp);
%% Calculate Derivatives
    % Smooth
    Xs = convn(Xn,G,'same');
    Ys = convn(Yn,G,'same');
    % First Derivatives
    Xp = convn(Xn,Gp,'same');
    Yp = convn(Yn,Gp,'same');
    % Second Derivatives
    Xpp = convn(Xn,Gpp,'same');
    Ypp = convn(Yn,Gpp,'same');
%% Calculate Curvature
    kappa = ( Xp.*Ypp - Yp.*Xpp ) ./ ( Xp.^2 + Yp.^2 ).^(3/2);
%% Crop     
    Xs = Xs(cop+1:length(X)+cop);
    Ys = Ys(cop+1:length(X)+cop);
    kappa = kappa(cop+1:length(X)+cop);
    %% Zero-Crossing points
%     zc = crossing(kappa);
    zc = crossing_adapted(kappa);
    %% Outputs
    if ~isempty(zc)
        K{i} = kappa;
        Xsm{i} = Xs;
        Ysm{i} = Ys;
        zcp{i} = zc;
        S{i} = sigma;
        flag = 1;
    end    
    %% Index ++
    sigma = sigma + sigmadelta;
    i = i + 1;
end
%% End
end