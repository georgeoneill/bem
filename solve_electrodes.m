function bem = solve_electrodes(surf,bem,els)



elpos = els.elecpos;
scalp = surf(1);

% The electrodes need to be projected onto the surface, just incase they
% are floating in space due to a ropey coregistration. Luckily fieldtrip
% has our back here with and does most of the legwork for us. 
[tri_weights, tri_idx] = project_onto_surface(elpos,scalp);

sol = zeros(length(elpos),size(bem.solution,1));
for ii = 1:length(elpos)
    tmp = bem.solution(scalp.tris(tri_idx(ii),:),:);
    sol(ii,:) = sum(bsxfun(@times,tmp,tri_weights(ii,:)'));
end

mults = [];
for ii = 1:length(surf)
    mults = cat(1,mults,bem.source_mult(ii)*ones(surf(ii).np,1));
end
mults = mults/(4*pi);
sol = bsxfun(@times,sol,mults');
bem.electrode_solution = sol;

end

function [weights, idx] = project_onto_surface(elrs,scalp)

geom = get_supplementary_geometry(scalp);
a = project_elec(elrs,scalp.rr,scalp.tris);
idx = a(:,1);
coords = zeros(length(idx),3);
for ii = 1:length(idx)
    coords(ii,:) = triangle_coordinates(elrs(ii,:),geom,idx(ii));
end
weights = [1-coords(:,1)-coords(:,2),coords(:,1),coords(:,2)];

end

function coo = triangle_coordinates(r,geom,id)

% unpack the variables we need 

tri_nn = geom.nn(id,:);
r1 = geom.r1(id,:);
r12 = geom.r12(id,:);
r13 = geom.r13(id,:);
a = geom.a(id,:);
b = geom.b(id,:);
c = geom.c(id,:);
rr = r-r1;
z = sum(rr.*tri_nn);
v1 = sum(rr.*r12);
v2 = sum(rr.*r13);
det = a.*b - c.^2;
x = (b*v1 - c*v2)/det;
y = (a*v2 - c*v1)/det;
coo = [x y z]; 

end

function geom = get_supplementary_geometry(surf)

r1 = surf.rr(surf.tris(:,1),:);
r12 = surf.rr(surf.tris(:,2),:) - r1;
r13 = surf.rr(surf.tris(:,3),:) - r1;
r1213(1,:,:) = r12;
r1213(2,:,:) = r13;
r1213 = permute(r1213,[2 1 3]);
a = sum(r12.*r12,2);
b = sum(r13.*r13,2);
c = sum(r12.*r13,2);

mat(1,1,:) = b;
mat(1,2,:) = -c;
mat(2,1,:) = -c;
mat(2,2,:) = a;

mat = permute(mat,[3 1 2]);
nrm = (a.*b - c.*c);
nrm(nrm==0) = 1;

mat = bsxfun(@rdivide,mat,nrm);

nn = fast_cross(r12,r13);
nn = bsxfun(@rdivide,nn,vnorm(nn,2)); 

geom.r1 = r1;
geom.r12 = r12;
geom.r13 = r13;
geom.r1213 = r1213;
geom.a = a;
geom.b = b;
geom.c = c;
geom.nn = nn;
geom.mat = mat;

end

function c = fast_cross(a,b)

assert(size(a,2)==3)
assert(size(b,2)==3)

c = [a(:,2).*b(:,3)-a(:,3).*b(:,2) ...
    a(:,3).*b(:,1)-a(:,1).*b(:,3)...
    a(:,1).*b(:,2)-a(:,2).*b(:,1)];
end