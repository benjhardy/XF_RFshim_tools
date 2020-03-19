function voxelation(fullheadMesh,path)
% Input: fullheadMesh matlab object that contains various mesh information
%        path - where to save the voxelized mesh
% Output: save a file to the workspace that contains the following
%       rho_m - the density matrix giving the density of the materials at
%       the mesh locations found in x_mesh, y_mesh, z_mesh
%       sigma_m epsilon_m - conductivity and permittivity in mesh
%       frankMask, brainMask, - total face mask, and total brain mask in
%       mesh
%       x_mesh, y_mesh, z_mesh hold the mesh locations

    
%     if isfile(fullfile(path,'voxelizedMesh.mat'))
%         fprintf('VoxelizedMesh already exists...\n')
%         return
%     end
    MeshExDensity = fullheadMesh.MeshExDensity;                             
    MeshExEpsilon_r = fullheadMesh.MeshExEpsilon_r;                           
    MeshExSigma = fullheadMesh.MeshExSigma;                           
    MeshEyDensity = fullheadMesh.MeshEyDensity;                          
    MeshEyEpsilon_r = fullheadMesh.MeshEyEpsilon_r;                             
    MeshEySigma = fullheadMesh.MeshEySigma;                     
    MeshEzDensity = fullheadMesh.MeshEzDensity;                       
    MeshEzEpsilon_r = fullheadMesh.MeshEzEpsilon_r;                        
    MeshEzSigma = fullheadMesh.MeshEzSigma;                   
    grid_X = fullheadMesh.grid_X;                
    grid_Y = fullheadMesh.grid_Y;                       
    grid_Z = fullheadMesh.grid_Z;
    
    % This is for any sort of averaged value like rho_m, sigma, epsilon
    ave_grid_x = grid_X(1:(end-1)) + diff(grid_X)/2;
    ave_grid_y = grid_Y(1:(end-1)) + diff(grid_Y)/2;
    ave_grid_z = grid_Z(1:(end-1)) + diff(grid_Z)/2;

    MeshExDensity = fullheadMesh.MeshExDensity;
    MeshEyDensity = fullheadMesh.MeshEyDensity;
    MeshEzDensity = fullheadMesh.MeshEzDensity;
    
    % Brain Mask
    % Density Mesh creation...
    load('dvb.mat','dvb_X','dvb_Y','dvb_Z')
    
    xD = ismember(MeshExDensity,dvb_X);
    yD = ismember(MeshEyDensity,dvb_Y);
    zD = ismember(MeshEzDensity,dvb_Z);

    rho_m = zeros(length(grid_Z)-1,length(grid_Y)-1,length(grid_X)-1);

    for i=1:length(grid_Z)-1
       for j = 1:length(grid_Y)-1
           for k = 1:length(grid_X)-1

           % for every combination
           rho_m(i,j,k) = ...
               (xD(i,j,k) + xD(i,j+1,k)+ xD(i,j,k+1)+ xD(i,j+1,k+1)...
                + yD(i,j,k) + yD(i+1,j,k) + yD(i,j,k+1)+ yD(i+1,j,k+1)...
                + zD(i,j,k) + zD(i+1,j,k) + zD(i,j+1,k)+ zD(i+1,j+1,k))/12;



           end
       end
    end
    
    rho_m(isnan(rho_m)) = 0;

      rho_ma = bwareaopen(rho_m,100000);
      brainMask = logical(rho_ma);
      brainMask = bwareaopen(brainMask,100000);
      brainMask = imfill(brainMask,'holes');
      
      % clean up the bottom part of the brain
      for i = 1:58
          slice = bwareaopen(squeeze(brainMask(i,:,:)),500);
          brainMask(i,:,:) = slice;
          
      end
     % volumeViewer(brainMask)
     % just to check if it shows up at least!
     %figure
     %imagesc(squeeze(brainMask(100,:,:)))
     %ti = sprintf('%s',path(end-16:end));
     %title(ti)
     %input('do you like?')
     
    % interp the brainMask onto the original grid with NN interp
    [x1,x2,x3] = ndgrid(ave_grid_z,ave_grid_y,ave_grid_x);
    F = griddedInterpolant(x1,x2,x3,double(brainMask),'nearest');
    [x1,x2,x3] = ndgrid(grid_Z,grid_Y,grid_X);
    brainMask_interp = F(x1,x2,x3);
    
    
    
    % frankMask
    % Density Mesh creation...
    xD = MeshExDensity;
    yD = MeshEyDensity;
    zD = MeshEzDensity;

    rho_m = zeros(length(grid_Z)-1,length(grid_Y)-1,length(grid_X)-1);

    for i=1:length(grid_Z)-1
       for j = 1:length(grid_Y)-1
           for k = 1:length(grid_X)-1

           % for every combination
           rho_m(i,j,k) = ...
               (xD(i,j,k) + xD(i,j+1,k)+ xD(i,j,k+1)+ xD(i,j+1,k+1)...
                + yD(i,j,k) + yD(i+1,j,k) + yD(i,j,k+1)+ yD(i+1,j,k+1)...
                + zD(i,j,k) + zD(i+1,j,k) + zD(i,j+1,k)+ zD(i+1,j+1,k))/12;



           end
       end
    end

    
    % frankMask is anything but freespace
    frankMask = ~isnan(rho_m);
    frankMask = imfill(frankMask,'holes');
    frankMask = bwareaopen(frankMask,1000);
    
    % interp frankMask
    [x1,x2,x3] = ndgrid(ave_grid_z,ave_grid_y,ave_grid_x);
    F = griddedInterpolant(x1,x2,x3,double(frankMask),'nearest');
    [x1,x2,x3] = ndgrid(grid_Z,grid_Y,grid_X);
    frankMask_i = F(x1,x2,x3);

    %zMovie(frankMask_i,100);


    % Sigma and Epsilon...
    % epsilon
    xE = MeshExEpsilon_r;
    yE = MeshEyEpsilon_r;
    zE = MeshEzEpsilon_r;
    % sigma
    xS = MeshExSigma;
    yS = MeshEySigma;
    zS = MeshEzSigma;

    sigma_m = zeros(length(grid_Z)-1,length(grid_Y)-1,length(grid_X)-1);
    epsilonr_m = zeros(length(grid_Z)-1,length(grid_Y)-1,length(grid_X)-1);

    for i=1:length(grid_Z)-1
       for j = 1:length(grid_Y)-1
           for k = 1:length(grid_X)-1

           % for every combination
           % sigma
           sigma_m(i,j,k) = ...
               (xS(i,j,k) + xS(i,j+1,k)+ xS(i,j,k+1)+ xS(i,j+1,k+1)...
                + yS(i,j,k) + yS(i+1,j,k) + yS(i,j,k+1)+ yS(i+1,j,k+1)...
                + zS(i,j,k) + zS(i+1,j,k) + zS(i,j+1,k)+ zS(i+1,j+1,k))/12;
           % for every combination
           % epsilonr
           epsilonr_m(i,j,k) = ...
               (xE(i,j,k) + xE(i,j+1,k)+ xE(i,j,k+1)+ xE(i,j+1,k+1)...
                + yE(i,j,k) + yE(i+1,j,k) + yE(i,j,k+1)+ yE(i+1,j,k+1)...
                + zE(i,j,k) + zE(i+1,j,k) + zE(i,j+1,k)+ zE(i+1,j+1,k))/12;




           end
       end
    end
    %
    epsilonr_m = frankMask.*epsilonr_m;
    %zMovie(epsilonr_m,150)
    
    %
    sigma_m = sigma_m.*frankMask;
    %zMovie(sigma_m,250)

    x_mesh = grid_X;
    y_mesh = grid_Y;
    z_mesh = grid_Z;
    %close all
    %clc
    %brainMask = brainMask_interp;
    
save(fullfile(path,'voxelizedMesh.mat'),'rho_m','sigma_m','epsilonr_m','frankMask','frankMask_i','brainMask','brainMask_interp','x_mesh','y_mesh','z_mesh')