function fixed = fwd_bem_ip_modify_solution(broken,ip,ip_mult,nps)


n_last = nps(end);
offsets = [0 cumsum(nps)];
mult = (1+ip_mult)/ip_mult;
disp('Combining')
for ii = 1:length(nps)
    % The BEM solution can be thought of as 3x3 block matrix, the right
    % most column of the matrix needs modifying block by block
    
    tmp = broken((offsets(ii)+1):(offsets(ii+1)),(offsets(length(nps))+1):end);
    tmp = tmp - 2*tmp*ip;
    broken((offsets(ii)+1):(offsets(ii+1)),(offsets(length(nps))+1):end) = tmp;
        
end

% special treatment for the bottom right block
tmp = broken((offsets(ii)+1):(offsets(ii+1)),(offsets(length(nps))+1):end);
tmp = tmp + mult*ip;
broken((offsets(ii)+1):(offsets(ii+1)),(offsets(length(nps))+1):end) = tmp;

fixed = ip_mult*broken;
