function [fwd bem] = go_solveBEM(opts,meshes,sources,varargin)

% check fieldtrip has been initialised (or we will run into trouble later)
ft = which('ft_defaults');
if isempty(ft)
    error('FieldTrip not found! Please add to path and type ft_defaults to initialise')
else
    ft_defaults
end

% do some sanity checks just to make sure options have been correctly set
if ~isfield(opts,'data')
    if ~isfield(opts,'coils') && ~isfield(opts,'els')
        error('Please provide a dataset location using opts.data, or specify a structure containing electrode or coil locations')
    end
end

switch length(meshes)
    case 3
        opts.cond = ft_getopt(opts,'cond',[0.3, 0.006, 0.3]);
    case 1
        opts.cond = ft_getopt(opts,'cond',0.3);
end

opts.meg = ft_getopt(opts,'meg',0);
opts.eeg = ft_getopt(opts,'eeg',0);

if ~opts.meg && ~opts.eeg
    warning('Modality not specified, assuming MEG only')
    opts.meg = 1;
end

if length(meshes) ~= length(opts.cond)
    error('Number of boundaries and condunctivity values do not match.')
end

if opts.eeg && length(opts.cond)==1
    error('You must have 3 boundaries for EEG solution!')
end

% FieldTrip likes to use pos/pnt/tri for mashes rather than faces/vertices
% but all the code I have written expects the latter. This function renames
% all the relevent items in each mesh structre to the correct name
meshes = rename_mesh_assets(meshes);

% Here is where checking the meshing will happen (to be written)
% Tests will include:
%               1) does the solid angle of each mesh come to +4pi sr?
%               2) is mesh 1 inside 2, which are both inside 3?
%               3) more tests can come.
opts.check_meshes = ft_getopt(opts,'check_meshes',1);
if opts.check_meshes
    for ii = 1:length(meshes)
        try
            [~, meshes(ii)] = check_meshes(meshes(ii),1);
        catch
            % Will try to repair meshes in case of failure
            fprintf('Initial sanity check of %s mesh failed, repairing...',meshes(ii).name)
            if isempty(which('meshcheckrepair'))
                ft_hastoolbox('iso2mesh',1)
            end
            [meshes(ii).vertices, meshes(ii).faces] = meshcheckrepair(meshes(ii).vertices, meshes(ii).faces);
            disp('checking again')
            check_meshes(meshes(ii),1);
        end
    end
    
    if length(meshes) > 1
        % check meshes are in the correct order, if not, flip!
        try
            check_meshes(meshes,2);
        catch
            disp('meshes may be nested in the wrong order, reversing...')
            tmp = meshes;
            meshes(1) = tmp(3);
            meshes(2) = tmp(2);
            meshes(3) = tmp(1);
            check_meshes(meshes,2);
        end
    end
    
end

% With all sanity checks passed, lets fill in all the surface information
% based on just the faces and vertices alone. Then solve the geometric part
% of the BEM common to both EEG and MEG (or skip if already precomputed)
surf = complete_surfaces(meshes);
if nargin==4
    bem = varargin{1};
else
    bem = solve_bem(surf,opts.cond);
end

% Add in the coil part of the solution in (for MEG)
if opts.meg
    if isfield(opts,'data') && ~isfield(opts,'coils')
        if isfield(opts,'coilaccuracy')
            sens = ft_read_sens(opts.data,'senstype','meg','coilaccuracy',opts.coilaccuracy);
        else
            sens = ft_read_sens(opts.data,'senstype','meg');
        end
        sens = ft_convert_units(sens,'m');
        coils.r = sens.coilpos;
        coils.o = sens.coilori;
        if isfield(sens,'tra')
            coils.w = sens.tra;
        end
    else
        coils = opts.coils;
    end
    bem = solve_coils(surf,bem,coils);
end

% Add in the electrodes part of the solution in (for EEG)
if opts.eeg
    opts.eegref = ft_getopt(opts,'eegref',0);
    if isfield(opts,'data') && ~isfield(opts,'els')
        elecs = ft_read_sens(opts.data,'senstype','eeg');            %% locations in cm
        elecs = ft_convert_units(elecs,'m');                          %% locations in m
        els.elecpos = elecs.chanpos;
        if opts.eegref
            fprintf('Adding EEG reference channel and balancing matrix\n');
            hdr = ft_read_header(opts.data);
            type = [hdr.orig.chs.kind];
            id = find(type==2);
            loc = hdr.orig.chs(id(1)).eeg_loc(:,2);
            els.elecpos = cat(1,els.elecpos,loc');
            els.w = ([eye(length(elecs.label)) -1*ones(length(elecs.label),1)]);
        end
    else
        els = opts.els;
    end
    bem = solve_electrodes(surf,bem,els);
end

% Some general housekeeping, there is probably a more elegent way to do
% this but it works for the time being.
tmp = bem;
clear bem
bem = struct;
bem.surfs = surf;
if isfield(tmp,'coil_solution')
    bem.coils = coils;
    bem.coil_solution = tmp.coil_solution;
end
if isfield(tmp,'electrode_solution')
    bem.els = els;
    bem.electrode_solution = tmp.electrode_solution;
end

% Do the field clculations!
L = solve_fields(sources,bem);

% Pack up results
if opts.meg
    if isfield(coils,'w')
        fwd.meg = coils.w*L.meg';
    else
        fwd.meg = L.meg';
    end
end

if opts.eeg
    if isfield(els,'w')
        fwd.eeg = els.w*L.eeg';
    else
        fwd.eeg = L.eeg';
    end
end

fwd.opts = opts;