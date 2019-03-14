function pass = check_meshes(meshes)

% Test 1: are the surfaces complete?
fprintf('Checking boundary mesh integrity: ')
for ii = 1:length(meshes)
    
    sr = get_solids_fast(meshes(ii).vertices,meshes(ii).faces);
    if abs(sr - (2*pi)) < 1e-5
        tmp_pass(ii) = 1;
    else
        tmp_pass(ii) = 0;
    end
    
end
pass(1) = sum(tmp_pass) == length(meshes);
if pass(1)
    fprintf('OK!\n')
else
    fprintf('FAIL\n')
    error('One or more boundaries are not completely closed surfaces');
end

% Test 2: Is surface 1 entirely within 2, which itself is in 3?
% Obvs, skip the test is there is a single layer BEM.
if length(meshes) > 1
    fprintf('Checking boundaries are entirely nested within each other: ')
    clear tmp_pass
    for ii = 1:length(meshes)-1
        fro = meshes(ii+1);
        to = meshes(ii);
        clear sr
        sr = get_solids_fast(fro.vertices,fro.faces,to.vertices);
        out = not(abs(sr - (2*pi)) < 1e5);
        tmp_pass(ii) = not(sum(out));
    end
    pass(2) = sum(tmp_pass) == length(meshes)-1;
    if pass(2)
        fprintf('OK!\n')
    else
        fprintf('FAIL\n')
        error('One or more are boundaries are not nested properly.');
    end
else
    pass(2) = 1;
end


end

% function [tot_solids, solids]= get_solids_slow(verts,faces,varargin)
%
% if nargin==1
%     fros = varargin{1};
%     npoints = length(fros);
% else
%     fros = mean(verts,1);
%     npoints = 1;
% end
% solids = zeros(1,length(faces));
%
% for ii = 1:length(faces)
%     tri_rrs = verts(faces(ii,:),:);
%     v1 = fros - tri_rrs(1,:);
%     v2 = fros - tri_rrs(2,:);
%     v3 = fros - tri_rrs(3,:);
%     triples = sum(fast_cross(v1,v2).*v3,2);
%     l1 = vnorm(v1,2);
%     l2 = vnorm(v2,2);
%     l3 = vnorm(v3,2);
%     ss = l1.*l2.*l3;
%     s1 = l3.*sum(v1.*v2,2);
%     s2 = l2.*sum(v1.*v3,2);
%     s3 = l1.*sum(v2.*v3,2);
%     ss = ss+s1+s2+s3;
%     solids(ii) = bsxfun(@atan2,triples,ss);
% end
%
% tot_solids = sum(-1*solids);
%
% end

function [solids] = get_solids_fast(verts,faces,varargin)

if nargin==3
    fros = varargin{1};
    npoints = length(fros);
else
    fros = mean(verts,1);
    npoints = 1;
end

% need to order the verices in order of all the triangles.

vtris = zeros(length(faces),3,3);
for ii = 1:length(faces)
    vtris(ii,:,:) = verts(faces(ii,:),:);
end

slices = arange(1,100,npoints);
if npoints == 1
    slices = [1 1];
end

solids = zeros(1,npoints);
for ii = 1:length(slices)-1
    tmp1 = reshape(fros(slices(ii):slices(ii+1),:),1,1,numel(slices(ii):slices(ii+1)),[]);
    tmp2 = reshape(permute(vtris,[2 1 3]),size(vtris,2),size(vtris,1),1,size(vtris,3));
    vs = bsxfun(@minus,tmp1,tmp2);
    v1 = squeeze(vs(1,:,:,:));
    v2 = squeeze(vs(2,:,:,:));
    v3 = squeeze(vs(3,:,:,:));
    triples = fast_cross_nd_sum(v1,v2,v3);
    ls = vnorm(vs,4);
    ss = squeeze(prod(ls,1));
    switch(ndims(v1))
        case 2 % just incase there is one point to compare to.
            s1 = ls(3,:,:)'.*sum(v1.*v2,ndims(v1));
            s2 = ls(2,:,:)'.*sum(v1.*v3,ndims(v1));
            s3 = ls(1,:,:)'.*sum(v2.*v3,ndims(v1));
            ss = ss'+s1+s2+s3;
            solids(slices(ii):slices(ii+1)) = -sum(bsxfun(@atan2,triples,ss));
        case 3
            s1 = squeeze(ls(3,:,:)).*sum(v1.*v2,ndims(v1));
            s2 = squeeze(ls(2,:,:)).*sum(v1.*v3,ndims(v1));
            s3 = squeeze(ls(1,:,:)).*sum(v2.*v3,ndims(v1));
            ss = ss+s1+s2+s3;
            solids(slices(ii):slices(ii+1)) = -sum(bsxfun(@atan2,triples,squeeze(ss)),1);
    end
end

end

function X = fast_cross(a,b)
i = a(:,2).*b(:,3)-a(:,3).*b(:,2);
j = a(:,3).*b(:,1)-a(:,1).*b(:,3);
k = a(:,1).*b(:,2)-a(:,2).*b(:,1);
X = [i j k];
end

function triples = fast_cross_nd_sum(a,b,c)
otherdims = repmat({':'},1,ndims(a)-1);
triples = (a(otherdims{:},2).*b(otherdims{:},3)-a(otherdims{:},3).*b(otherdims{:},2)).*c(otherdims{:},1)...
    + (a(otherdims{:},3).*b(otherdims{:},1)-a(otherdims{:},1).*b(otherdims{:},3)).*c(otherdims{:},2)...
    + (a(otherdims{:},1).*b(otherdims{:},2)-a(otherdims{:},2).*b(otherdims{:},1)).*c(otherdims{:},3);
end

function vec = arange(lo,step,hi)
vec = lo:step:hi;
if vec(end) ~= hi
    vec = [vec hi];
end
end