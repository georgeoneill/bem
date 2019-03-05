function surf = completesurfaces(meshes)

% Function to do some bookeeping on the previously generated surfaces,
% assumes that these are correct (ie normals point outwards, meshes have a
% solid angle of 4pi sr etc.)


type = {'iskull','oskull','scalp'};


for ii = 1:length(meshes)
    
    surf(ii).name           = type(ii);
    surf(ii).rr             = meshes(ii).vertices;
    surf(ii).tris           = meshes(ii).faces;
    surf(ii).ntri           = length(meshes(ii).faces);
    surf(ii).np             = length(meshes(ii).vertices);
    
    rr      = meshes(ii).vertices;
    tris    = meshes(ii).faces;
    r1      = rr(tris(:,1),:);
    r2      = rr(tris(:,2),:);
    r3      = rr(tris(:,3),:);
    
    surf(ii).tri_cent       = (r1+r2+r3)/3;
    tmp                     = cross((r2-r1),(r3-r1));
    surf(ii).tri_area       = 0.5*(vnorm(tmp,2));
    surf(ii).tri_nrms       = tmp./(2*surf(ii).tri_area*ones(1,3));
    surf(ii).tri_neighbours = triangle_neighbours(surf(ii).tris,surf(ii).np);
    surf(ii).nrms           = normstri2vert(surf(ii).tris,surf(ii).tri_nrms,surf(ii).tri_area,surf(ii).np);
    
    
end

% We want the flip the order so the scalp is the first mesh in the 3 shell
% case, just to keep convention with MNE's way of doing things.

if length(meshes) == 3
    tmp = surf;
    surf = struct;
    surf = tmp(3);
    surf(2) = tmp(2);
    surf(3) = tmp(1);
end