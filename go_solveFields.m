function fwd = go_solveFields(bem,sources)

% If the BEM has already been solved for both the boudries and sensor
% locations, then we can just compute the fields.

% Do the field clculations!
L = solve_fields(sources,bem);

% Pack up results
if isfield(L,'meg')
    if isfield(bem.coils,'w')
        fwd.meg = bem.coils.w*L.meg';
    else
        fwd.meg = L.meg';
    end
end

if isfield(L,'eeg')
    if isfield(bem.els,'w')
        fwd.eeg = bem.els.w*L.eeg';
    else
        fwd.eeg = L.eeg';
    end
end

end
