function solids_inv = fwd_bem_multi_solution(solids, gamma, nps)

nsurf = length(nps);
pi2 = 1/(2*pi);
n_tot = sum(nps);
defl = 1/n_tot;
offsets = [0 cumsum(nps)];

tmp = zeros(size(solids));

for ii = 1:nsurf
    for jj = 1:nsurf
        
        if ~isempty(gamma);
            mult = pi2*gamma(ii,jj);
        else
            mult = pi2;
        end
        tmp((offsets(ii)+1):(offsets(ii+1)),(offsets(jj)+1):(offsets(jj+1))) =...
            defl - mult*solids((offsets(ii)+1):(offsets(ii+1)),(offsets(jj)+1):(offsets(jj+1)));
    end
end
   
tmp = tmp+eye(size(tmp));
solids_inv = inv(tmp);

