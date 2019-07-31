function coeff = fwd_bem_lin_pot_coeff(surf)

% Calculate the coefficients for linear collocation approach.



nsurfs = length(surf);
nps = [surf.np];
coeff = zeros(sum([surf.np]));
offsets = [0 cumsum(nps)];

for ii = 1:nsurfs
    
    rr_ord = 1:nps(ii);
    
    for jj = 1:nsurfs
        
        tmp = zeros(surf(ii).np,surf(jj).np);
        fprintf('lin pot coeffs: %s -> %s\n',cell2mat(surf(ii).name),cell2mat(surf(jj).name));
        tri_nrms    = surf(jj).tri_nrms;
        tri_area    = surf(jj).tri_area;
        
        for kk = 1:surf(jj).ntri
            tri         = surf(jj).tris(kk,:);
            tri_rr      = surf(jj).rr(tri,:);
            if ii == jj
                skip_idx = rr_ord == tri(1) | rr_ord == tri(2) | rr_ord == tri(3);
            else 
                skip_idx = [];
            end
            coeffs = lin_pot_coeff(surf(ii).rr, tri_rr, tri_nrms(kk,:), tri_area(kk));
            coeffs(skip_idx,:) = 0;
            tmp(:,tri) = tmp(:,tri) - coeffs;
        end
        
        if ii == jj
            tmp = correct_auto_elements(surf(ii),tmp);
        end
        
        % check for NaNs or Infs.
        chk = isnan(tmp) | isinf(tmp);
        if sum(chk(:)) > 0
            error('NaNs or Infs detected in potential coefficients, check meshes!')
        end
        
        coeff((offsets(ii)+1):(offsets(ii+1)),(offsets(jj)+1):(offsets(jj+1))) = tmp;
    end

end
end

function omega = lin_pot_coeff(fros, tri_rr, tri_nrms, tri_area)

omega = zeros(size(fros));

v1 = ones(length(fros),1)*tri_rr(1,:) - fros;
v2 = ones(length(fros),1)*tri_rr(2,:) - fros;
v3 = ones(length(fros),1)*tri_rr(3,:) - fros;
triples = fast_cross_nd_sum(v1,v2,v3);
l1 = vnorm(v1,2);
l2 = vnorm(v2,2);
l3 = vnorm(v3,2);
ss = l1.*l2.*l3;
s1 = l3.*sum(v1.*v2,2);
s2 = l2.*sum(v1.*v3,2);
s3 = l1.*sum(v2.*v3,2);
ss = ss+s1+s2+s3;
% ss = ss + l3.*sum(v1.*v2,2);
% ss = ss + l2.*sum(v1.*v3,2);
% ss = ss + l1.*sum(v2.*v3,2);
solids = zeros(size(triples));
solids = bsxfun(@atan2,triples,ss);
bad_mask = (abs(solids) < (pi/1e6));
l1(bad_mask) = 1;
l2(bad_mask) = 1;
l3(bad_mask) = 1;

beta = zeros(length(fros),3);
beta(:,1) = calc_beta(v1,l1,v2,l2);
beta(:,2) = calc_beta(v2,l2,v3,l3);
beta(:,3) = calc_beta(v3,l3,v1,l1);


vec_omega = ((beta(:,3) - beta(:,1))*ones(1,3)).*v1;
vec_omega = vec_omega + ((beta(:,1) - beta(:,2))*ones(1,3)).*v2;
vec_omega = vec_omega + ((beta(:,2) - beta(:,3))*ones(1,3)).*v3;

area2 = 2*tri_area;
n2 = 1/(area2.*area2);
yys = cat(3,v1,v2,v3);
idx = [3 1 2; 2 3 1];

for ii = 1:3
    diff = yys(:,:,idx(1,ii)) - yys(:,:,idx(2,ii));
    zdots = fast_cross_nd_sum(yys(:,:,idx(2,ii)),yys(:,:,idx(1,ii)),tri_nrms);
    omega(:,ii) = -n2 * ( 2 * area2 * zdots .* solids - triples .* sum(diff.*vec_omega,2));
end
omega(bad_mask,:) = 0;

end

function mat = correct_auto_elements(surf,mat)
pi2 = 2*pi;
misses = pi2 - sum(mat,2);
for ii = 1:length(misses);
    miss = misses(ii);
    n_memb = length(surf.tri_neighbours{ii});
    mat(ii,ii) = miss/2;
    miss = miss/(2*n_memb);
    neighbours = unique(surf.tris(surf.tri_neighbours{ii},:));
    neighbours(neighbours==ii) = [];
    for jj = 1:length(neighbours)   
        mat(ii,neighbours(jj)) = mat(ii,neighbours(jj))+miss;
    end
end
end

function triples = fast_cross_nd_sum(a,b,c)
    triples = (a(:,2).*b(:,3)-a(:,3).*b(:,2)).*c(:,1)...
        + (a(:,3).*b(:,1)-a(:,1).*b(:,3)).*c(:,2)...
        + (a(:,1).*b(:,2)-a(:,2).*b(:,1)).*c(:,3);
end

function beta = calc_beta(rk,rk_norm,rk1,rk1_norm)
rkk1    = rk1(1,:) - rk(1,:);
sz      = vnorm(rkk1,2);
rkk1    = rkk1./sz;
tmp     = sum(bsxfun(@times,[rk;rk1],rkk1),2);
num     = rk_norm + tmp(1:length(rk));
den     = rk1_norm + tmp(length(rk)+1:end);
tmp1    = num./den;
tmp2    = log(tmp1);
beta    = tmp2./sz;
end