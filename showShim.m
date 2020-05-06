function showShim(dir,w_full)

    % just a subset of proc_shim code but keeps enough to show the shim
    Nc = size(w_full,1);
    
    grids = cell(Nc,1);
    for j=1:Nc
             
        % Check the fields have already been extracted
        if ( ~isfile(fullfile(dir,'B1plus_folder',sprintf('B1plus_%d.mat',j))))
            fprintf('Bplus_%d is missing ...\n',j)
            return
        end
        % Extract the grid for each B_Field result
        % to be passed to consolidate B_field function
        [grid_X,grid_Y,grid_Z] = extractGrid(dir,j);
        grids{j} = {grid_X,grid_Y,grid_Z}; 
    end
        
    % Consolidate the B_Field values
    [B1plus_master,master_z,master_y,master_x,iz,iy,ix] = ConsolidateFields(grids,dir,Nc);        

    %
    % Reform
    [B1plus_m,FOX] = B1plusReform(B1plus_master,2,ix,iy,iz);
    
    % show the shim now
    tmp = 0*FOX;
    tmp(FOX) = abs(B1plus_m*w_full);
    zMovie(tmp,100,2)
    





end