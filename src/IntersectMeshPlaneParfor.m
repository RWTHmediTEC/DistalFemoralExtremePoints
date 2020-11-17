function Contours = IntersectMeshPlaneParfor(Mesh, PlaneOrigins, PlaneNormals)
%
% AUTHOR: Maximilian C. M. Fischer
% COPYRIGHT (C) 2020 Maximilian C. M. Fischer
% LICENSE: EUPL v1.2
%

if size(PlaneOrigins,2)~=3 
    error('Size of PlaneOrigins has to be Nx3')
end

if size(PlaneNormals,2)~=3 
    error('Size of PlaneNormals has to be Nx3')
end

if size(PlaneNormals,1)~=size(PlaneOrigins,1) && size(PlaneNormals,1)~=1
    error('Size of PlaneNormals has to the same as PlaneOrigins or 1x3')
end

if size(PlaneNormals,1)==1
    PlaneNormals = repmat(PlaneNormals, size(PlaneOrigins,1),1);
end

Contours = cell(size(PlaneOrigins,1),1);
parfor p=1:size(PlaneOrigins,1)
    % To speed things up, use only the faces in the cutting plane as input for intersectPlaneSurf
    Plane = createPlane(PlaneOrigins(p,:), PlaneNormals(p,:));
    MeshInPlane = cutMeshByPlane(Mesh, Plane, 'part','in');
    Contours{p} = intersectPlaneSurf(MeshInPlane, PlaneOrigins(p,:), PlaneNormals(p,:));
end; clear p
