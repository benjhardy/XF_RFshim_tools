function [B1plus_master,master_z,master_y,master_x,iz,iy,ix] = ConsolidateFields(grids,workingDir,Nc)
% Intent: takes the field maps from multiple different coils and only
%         pushes through the elements points that align with the head
% 
%         Input: 
%         grids - cell(Nc,1){1,3}(grid_X,grid_Y,grid_Z)
%         B1plus_folder - contains the B1plus_# for each coil
%         Nc - Number of coils
%         Output:
%         B1plus_master - as it seems
%         master_z, master_y, master_x - master cartesian locations, this
%         is also the order of the linear index of the B1plus_master

    cd(workingDir)
    % first load mesh_1
    load('mesh_1.mat','grid_X','grid_Y','grid_Z');
    path_B = sprintf('B1plus_folder/B1plus_%d',1);
    load(path_B,'xpts','ypts','zpts')
    % master grid
    master_x = grid_X(xpts);
    master_y = grid_Y(ypts);
    master_z = grid_Z(zpts);
    iz = zpts;
    iy = ypts;
    ix = xpts;
    
    % initiate to zero
    B1plus_master = zeros(numel(master_x)*numel(master_y)*numel(master_z),Nc);
    % go to each coil
    for i=1:Nc
        
        
        path_B = sprintf('B1plus_folder/B1plus_%d',i);
        % load the points in the grid and the index in these grids . . .
        grid_X = grids{i,1}{1,1};
        grid_Y = grids{i,1}{1,2};
        grid_Z = grids{i,1}{1,3};
        load(path_B,'xpts','ypts','zpts')
        
        % Begin Interpolation
        fprintf('Interpolating Coil map #%d...\n',i)
        [x1,x2,x3] = ndgrid(grid_Z(zpts),grid_Y(ypts),grid_X(xpts));
        load(path_B,'B1plus_m')
        Bi = griddedInterpolant(x1,x2,x3,B1plus_m);
        [x1,x2,x3] = ndgrid(master_z,master_y,master_x);
        B1plus_i = Bi(x1,x2,x3);
        
        % add to the master file
        B1plus_master(:,i) = B1plus_i(:);
            
    end    

end