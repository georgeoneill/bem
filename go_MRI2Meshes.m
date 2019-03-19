function meshes = go_MRI2Meshes(opts,mripath)

opts.method = ft_getopt(opts,'method','fieldtrip');
opts.units  = ft_getopt(opts,'units','m');

% slsh = strfind(mripath,'/');
% if isempty(slsh)
%     base = './';
% else
%     base = mripath(1:slsh(end));
% end

switch lower(opts.method)
    case 'fieldtrip'
        warning('method not yet checked, hold tight!')
        
        mri = ft_read_mri(mripath);
        
        cfg           = [];
        cfg.output    = {'brain','skull','scalp'};
        segmentedmri  = ft_volumesegment(cfg, mri);
        
        cfg=[];
        cfg.tissue={'scalp','skull','brain'};
        cfg.numvertices = [2500 2500 2500];
        bnd=ft_prepare_mesh(cfg,segmentedmri);
        
        try
            meshes(1).pnt =  bnd(3).pos;
            meshes(2).pnt =  bnd(2).pos;
            meshes(3).pnt =  bnd(1).pos;
        catch
            meshes(1).pnt =  bnd(3).pnt;
            meshes(2).pnt =  bnd(2).pnt;
            meshes(3).pnt =  bnd(1).pnt;
        end
        meshes(1).tri = bnd(3).tri;
        meshes(2).tri = bnd(2).tri;
        meshes(3).tri = bnd(1).tri;
        
        meshes(1).name = 'iskull';
        meshes(2).name = 'oskull';
        meshes(3).name = 'scalp';
        
        meshes(1).unit = bnd(3).unit;
        meshes(2).unit = bnd(2).unit;
        meshes(3).unit = bnd(1).unit;
        
    case 'spm'
        
        % check a full version of SPM12 has been added to the path,
        % as fieldtrip's light version of SPM12 is not extensive enough to
        % support using SPM's canonical meshes.
        if isempty(which('spm_eeg_inv_mesh'))
            error(['To use the SPM method of boundary exctaction, '...
                'you need a complete version of SPM12 added to your path'])
        end
        
        % Sanity checks of the anatomical in quesstion
        mri = ft_read_mri(mripath);
        mri = ft_checkdata(mri, 'datatype', 'volume', 'feedback', 'yes', 'hasunit', 'yes', 'hascoordsys', 'yes');
        mri = ft_convert_units(mri, 'mm');
        mri = ft_convert_coordsys(mri, 'spm');
        % flip and permute the 3D volume itself, so that the voxel and
        % headcoordinates approximately correspond this improves the convergence
        % of the segmentation algorithm
        [mri] = align_ijk2xyz(mri);
        nom = tempname;
        Va = ft_write_mri([nom '.img'], mri.anatomy, 'transform', mri.transform);
        
        % We need a full version of SPM12 for this to work, however, fieldtrip
        % comes with elements of spm5/8/12/* to serves its own purposes, we
        % need to removed these from the path for the time being so we do
        % not end up in compatibillity hell.
        spm5 = fullfile(fileparts(which('ft_defaults')), 'external','spm5');
        rmpath(genpath(spm5))
        spm8 = fullfile(fileparts(which('ft_defaults')), 'external','spm8');
        rmpath(genpath(spm8))
        spm12 = fullfile(fileparts(which('ft_defaults')), 'external','spm12');
        rmpath(genpath(spm12))
        
        % SHOWTIME!
        spmesh = spm_eeg_inv_mesh([nom '.img'],2);
        % import meshes
        mesh_iskull_tmp = gifti(spmesh.tess_iskull);
        mesh_oskull_tmp = gifti(spmesh.tess_oskull);
        mesh_scalp_tmp = gifti(spmesh.tess_scalp);
        
        % warp to the coordinate frame of the original MRI
        tmp = mesh_iskull_tmp.vertices;
        tmp1 = ft_warp_apply(mri.head2headOrig,tmp);
        
        meshes(1).tri = mesh_iskull_tmp.faces;
        meshes(1).pnt = tmp1;
        meshes(1).name = 'iskull';
        meshes(1).unit = 'mm';
        
        tmp = mesh_oskull_tmp.vertices;
        tmp1 = ft_warp_apply(mri.head2headOrig,tmp);
        
        meshes(2).tri = mesh_oskull_tmp.faces;
        meshes(2).pnt = tmp1;
        meshes(2).name = 'oskull';
        meshes(2).unit = 'mm';
        
        tmp = mesh_scalp_tmp.vertices;
        tmp1 = ft_warp_apply(mri.head2headOrig,tmp);
        
        meshes(3).tri = mesh_scalp_tmp.faces;
        meshes(3).pnt = tmp1;
        meshes(3).name = 'scalp';
        meshes(3).unit = 'mm';
        
        % flush all SPM temp files when done
        if ispc
            system(['del ' nom '*']);
        else
            system(['rm ' nom '*']);
        end
        
    otherwise
        error('boudnary generation method not supported, select either fieldtrip or spm')
end

% Change units into m (default) or user specified.
meshes(1) = ft_convert_units(meshes(1),opts.units);
meshes(2) = ft_convert_units(meshes(2),opts.units);
meshes(3) = ft_convert_units(meshes(3),opts.units);

end