function cfh = BOMultiScaleCurvaturePlot2D_adapted(K,S,X,Y,Xsm,Ysm,zcp,p)

%% ADAPTED - NOT THE ORIGINAL FILE
% - first figure handle output added
% - third figure added

%% BOMultiScaleCurvaturePlot2D - plots curvature scale-space image
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
%       K           - curvatures calculated for each sigma
%       S           - sigmas
%       X           - x coor of contour
%       Y           - y coor of contour
%       Xsm         - x coor of contour convolved with Gaussian
%       Ysm         - x coor of contour convolved with Gaussian
%       zcp         - indexes of zero-crossing points
%       p           - pause
%
%   OUTPUT:
%       out     - image with detected points
%       zcp     - zero-crossing points
%
%   USAGE:
% 
%   AUTHOR:
%       Boguslaw Obara, http://boguslawobara.net/
%
%   VERSION:
%       0.1 - 30/06/2009 First implementation

%% Plot: Contour
cfh = figure('name', 'Contour'); 
for i=1:length(zcp)
    Xs = Xsm{i};
    Ys = Ysm{i};
    zc = zcp{i};
    plot(X,Y,'k-','LineWidth',2);
    hold on;
    plot(Xs,Ys,'r-');
    hold on;
    plot(Xs(zc),Ys(zc),'go');
    pause(p);
    hold off
end
plot(X,Y,'k-','LineWidth',2);
title('Contour');
axis equal
%% Plot: Curvature Scale Space (CSS) Image
% figure('name', 'Curvature Scale Space (CSS) Image'),
% for i=1:length(zcp)
%     plot(zcp{i},S{i},'k.');%,'MarkerFaceColor',C(i,:),'MarkerEdgeColor',C(i,:));
%     hold on
% end
% hold off    
% xlim([1 length(K{1})]);
% title('Curvature Scale Space (CSS) Image');
%% Plot: Curvature kappa
% figure
% surf(cell2mat(K))
% xlabel('\sigma');ylabel('u');zlabel('\kappa(u,\sigma)');
% title('Curvature of \Gamma_\sigma: \kappa(u,\sigma)');
%% End
end