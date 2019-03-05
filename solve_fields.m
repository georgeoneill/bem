function L = solve_fields(locs,bem)

% Right this is a multimodal solver EEG/MEG so we need to do a bit of
% housekeeping to ensure that we dont mix up what *type* of field we are
% solving for...

type = isfield(bem,{'coil_solution','electrode_solution'});
ismeg = type(1);
iseeg = type(2);

if xor(ismeg,iseeg); % Checks if only one we can pass it through now.
    if ismeg
        L.meg = compute_meeg_fields(locs,bem);
    else
        L.eeg = compute_meeg_fields(locs,bem);
    end
else
    
    % Quick sanity check to see if both have been solved for, instead of
    % neither
    if ~sum(type)
        error('You appear to not have solved BEM for either electrods or coils!')
    end
    
    % calculate indvidually, starting with MEG
    bem_tmp = struct;
    bem_tmp.surfs = bem.surfs;
    bem_tmp.coils = bem.coils;
    bem_tmp.coil_solution = bem.coil_solution;
    L.meg = compute_meeg_fields(locs,bem_tmp);
    
    % now EEG.
    bem_tmp = struct;
    bem_tmp.surfs = bem.surfs;
    bem_tmp.els = bem.els;
    bem_tmp.electrode_solution = bem.electrode_solution;
    L.eeg = compute_meeg_fields(locs,bem_tmp);
    
end

end

function L = compute_meeg_fields(locs,bem)

bem_rr = [];
for ii = 1:length(bem.surfs)
    bem_rr = cat(1,bem_rr,bem.surfs(ii).rr);
end

if isfield(bem,{'coil_solution'})
    sol = bem.coil_solution';
    V = zeros(3*size(locs,1),length(bem.coils.r));
    disp('Calculating MEG lead fields ');
else
    sol = bem.electrode_solution';
    V = zeros(3*size(locs,1),length(bem.els.elecpos));
    disp('Calculating EEG lead fields ');
end

chunks = 250;
bounds = arange(0,chunks,size(locs,1));
bounds3 = arange(0,3*chunks,3*size(locs,1));


% Lets solve the potential half of Geleowitz's Eqn.
tic
parfor_progress(length(bounds)-1);
for ii = 1:(length(bounds)-1)
    v0 = infinite_potentials(locs((bounds(ii)+1):bounds(ii+1),:),bem_rr);
    v0 = reshape(v0,size(v0,1),3*size(v0,3))';
    V((bounds3(ii)+1):bounds3(ii+1),:) = v0*sol;
    parfor_progress;
end
parfor_progress(0);
toc

if isfield(bem,{'coil_solution'})
    % And now the primary current at the coils if MEG sensors.
    b0 = zeros(size(V,1),length(bem.coils.r));
    for ii = 1:length(bem.coils.r)
        tmp = infinite_Bfields(locs,bem.coils.r(ii,:),bem.coils.o(ii,:));
        b0(:,ii) = reshape(tmp,numel(tmp),1);
    end
    B = V+b0;
    L = B*1e-7;
else
    L = V;
end

end

function V = infinite_potentials(locs,bem_rr)
locs2(1,:,:) = locs;
locs2 = permute(locs2,[1 3 2]);
diff = bsxfun(@minus,bem_rr,locs2);
diff_norm = sum(diff.^2,2);
diff_norm = diff_norm.*sqrt(diff_norm);
diff_norm(diff_norm==0) = 1;
V = bsxfun(@rdivide,diff,diff_norm);
end

function B = infinite_Bfields(locs,cr,co)
diff = bsxfun(@minus,cr,locs);
diff_norm = sum(diff.^2,2);
diff_norm = diff_norm.*sqrt(diff_norm);
diff_norm(diff_norm==0) = 1;
x = fast_cross(diff,co);
B = bsxfun(@rdivide,x,diff_norm)';
end

function vec = arange(lo,step,hi)
vec = lo:step:hi;
if vec(end) ~= hi
    vec = [vec hi];
end
end

function c = fast_cross(a,b)

assert(size(a,2)==3)
assert(size(b,2)==3)

c = [a(:,2).*b(:,3)-a(:,3).*b(:,2) ...
    a(:,3).*b(:,1)-a(:,1).*b(:,3)...
    a(:,1).*b(:,2)-a(:,2).*b(:,1)];
end