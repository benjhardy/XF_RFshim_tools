function proc_shim(path,vv,betas,Nc,stopWhen,phase)
% WRAPPER to perform shims when the simulations have been performed!
% Inputs:   path: the fullpath to the Simulations and data folder where the
%                 b1plus_folder will be made along with the voxelizedMesh/mask
%                 information
%           vv: the voxelation function version to be run. v0 is 
%               interactive in terms of making the mask, v1(voxelation) 
%               works for batch_1 unmodfified Yurui Gao data sets.
%           betas: the regularization parameter weighing the value of
%                 the power term in the cost function [1 x n]
% If any statement in a try block generates an error, program control
    % goes immediately to the catch block, which contains your error handling statements.
    try
        home = path;
        % Step 1: Mesh information, extracted to find mask, conductivity and
        % permittivity
        if (~isfile(fullfile(home,'voxelizedMesh.mat')))
            
            mesh = matfile(fullfile(home,'mesh_1.mat'));

            % voxelation requires some interaction, maybe do this for the first head
            % and then generalize somehow?
            % I don't like this as I need to be able to run this 3000 times without
            % having to designate the mask as hands on as I do presently in voxelation
            if ( vv == 0 )
                % this first version lets the user select the roi values to
                % base the brain mask
                voxelation_v0(mesh,home)
            else
                voxelation(mesh,home)
            end
        else
            fprintf('voxelizedMesh.mat already exists in %s\n',home)
        end

        % Some B-Field grid and B-Field Extraction
        % 
        
        proj = home;
        
        % Step 0: (optional) Renumber any abnormal numbering scheme in the simulations folder
        % number order.
        % renumber_sims('H:\256_ellipse_staggered_r.xf')

        grids = cell(Nc,1);
        for j=1:Nc
            % Step 1a: Read B field results with export_B_Field
            % change to an if statement ....
            if ( isfile(fullfile(proj,'B1plus_folder',sprintf('B1plus_%d.mat',j))))
                fprintf('Bplus_%d already extracted ...\n',j)
            else
                export_B_Field(proj,j,path);
            end
            % Step 1b: Extract the grid for each B_Field result
            % to be passed to consolidate B_field function
            [grid_X,grid_Y,grid_Z] = extractGrid(proj,j);
            grids{j} = {grid_X,grid_Y,grid_Z}; 
        end
        %
        % Step 2: Consolidate the B_Field values
        [B1plus_master,master_z,master_y,master_x,iz,iy,ix] = ConsolidateFields(grids,path,Nc);
        indexMesh_z = iz ; indexMesh_y = iy; indexMesh_x = ix;
        save(fullfile(path,'Indices.mat'),'indexMesh_z','indexMesh_y','indexMesh_x')

        %
        % Step 3: Reform the field???
        [B1plus_m,FOX] = B1plusReform(B1plus_master,2,ix,iy,iz);
        % TikhonovShim

        

        
        
        d = 1;
        % Begin the shimming process...
        fprintf('Shimming the B1plus fields for %d element array.\n',Nc)
        [wfull_m, RSD_v, maxI_v, absps_v,cpuTime] = ...
                 tikhonovShim(B1plus_m, betas, stopWhen, FOX, phase, d, path);

        fprintf('Finished Shim for %d element array.\n',Nc)

        % Export/save the data . . . .
        filename = sprintf('L-Curve_%del_Shim_%s_betas_%.0d.mat',Nc,phase,betas(1));         
        save(fullfile(home,filename), 'wfull_m', 'RSD_v', 'maxI_v','absps_v','betas','stopWhen','phase','d','cpuTime');
        fprintf('Results Saved to %s\n',path)
        cd(home)
        
    catch me
        fprintf('...............\nSomething went wrong with iteration %s\n.......\n...............\n',path)
        me
    end








end