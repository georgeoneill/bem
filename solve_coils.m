function bem = solve_coils(surf,bem,coils)

offsets = [0 cumsum([surf.np])];
sol = zeros(length(coils.r),sum([surf.np]));
for ii = 1:length(surf)
    tic
    coeff = lin_field_coeff(surf(ii),coils);
    coeff = coeff*bem.field_mult(ii);
    tmp = bem.solution((offsets(ii)+1):(offsets(ii+1)),:);
    sol = sol + coeff*tmp;
    toc
end

mults = [];
for ii = 1:length(surf)
    mults = cat(1,mults,bem.source_mult(ii)*ones(surf(ii).np,1));
end
mults = mults/(4*pi);
sol = bsxfun(@times,sol,mults');
bem.coil_solution = sol;
end

function coeff = lin_field_coeff(surf,coils)

o = coils.o;
r = coils.r;


bem_rr = surf.rr;
tris = surf.tris;
tn = surf.tri_nrms;
ta = surf.tri_area;
coeff = zeros(length(r),length(bem_rr));

for ii = 1:surf.ntri
    tri = tris(ii,:);
    tri_rr = bem_rr(tri,:);
    tnn = tn(ii,:);
    zz = zeros(length(r),3);
    for jj = 1:3
        trr = tri_rr(jj,:);
        diff = bsxfun(@minus,r,trr);
        d2 = sum(diff.*diff,2);
        c = fast_cross(diff,tnn);
        x = ta(ii)*sum(c.*o,2)./(3.*d2.*sqrt(d2));
        zz(:,jj) = x;
    end
coeff(:,tri) = coeff(:,tri)+zz;

end
end

function c = fast_cross(a,b)

assert(size(a,2)==3)
assert(size(b,2)==3)

c = [a(:,2).*b(:,3)-a(:,3).*b(:,2) ...
    a(:,3).*b(:,1)-a(:,1).*b(:,3)...
    a(:,1).*b(:,2)-a(:,2).*b(:,1)];
end
