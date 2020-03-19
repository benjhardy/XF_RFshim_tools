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
    % load('voxelizedMesh.mat','rho_m')
    % vals = nonzeros(rho_m);
    
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
    % 
    ave_grid_x = grid_X(1:(end-1)) + diff(grid_X)/2;
    ave_grid_y = grid_Y(1:(end-1)) + diff(grid_Y)/2;
    ave_grid_z = grid_Z(1:(end-1)) + diff(grid_Z)/2;
    %load('C:\Users\benja\OneDrive - Vanderbilt\Documents\MATLAB\RFshimTools\densityValsBrain.mat', 'densityValsBrain')
    
    % load('F:\5634\yuruiDensityVals.mat','yuruiDensityVals')
    %densityValsBrain = yuruiDensityVals;
    %densityValsBrain = 911;
    MeshExDensity = fullheadMesh.MeshExDensity;
    MeshEyDensity = fullheadMesh.MeshEyDensity;
    MeshEzDensity = fullheadMesh.MeshEzDensity;
    % Density Mesh creation...
    
    xD = MeshExDensity;
    yD = MeshEyDensity;
    zD = MeshEzDensity;
    xD(isnan(xD))=0;
    yD(isnan(yD))=0;
    zD(isnan(zD))=0;

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
    %volumeViewer(rho_m)
    % select the density values of interest to create the brainMask
    figure
    imagesc(squeeze(rho_m(100,:,:)))
    % find the vals in densityValsBrain
    a = roipoly;
    colorbar
    %
    vals = unique(a.*squeeze(rho_m(100,:,:)));
    vals(isnan(vals)) = 0;
    vals1 = unique(vals);
    
    figure
    imagesc(squeeze(rho_m(70,:,:)))
    % find the vals in densityValsBrain
    a = roipoly;
    colorbar
    %
    vals = unique(a.*squeeze(rho_m(70,:,:)));
    vals(isnan(vals)) = 0;
    vals2 = unique(vals);
    
    vals = [vals1;vals2];

    % ****** LOAD this first
%     figure
%     title('Select brain and cerebellum matter density vals.')
%     imagesc(squeeze(rho_m(65,:,:)))
%     dvb = roipoly;
%     d = rho_m(65,dvb);
%       vals = d(:);
       i = ismember(rho_m,vals);
       brainMask = rho_m.*i;
       brainMask(isnan(brainMask)) = 0;
%     %
      %brainMask = logical(rho_m);
      brainMask = bwareaopen(brainMask,1000);
      brainMask = imfill(brainMask,'holes');
      volumeViewer(brainMask)
%    
%     
    %cancel the bottom z slices
%     inputz = 81;
%     slice = 81;
%     prompt = 'Insert a slice number to observe (0 to exit) plz remember to look for eye top, mid, and bottom: \n';
%     
%     %enter = 0;
% %     while ent.Callback == false
% %         f = figure;
% %         ax = axes('Parent',f,'position',[0.13 0.39  0.77 0.54]);
% %         h = imagesc(squeeze(brainMask(slice,:,:)));
% %         b = uicontrol('Parent',f,'Style','slider','Position',[81,54,419,23],...
% %                   'value',slice, 'min',1, 'max',size(brainMask,1));
% %         bgcolor = f.Color;
% % 
% %         ent = uicontrol;
% %         ent.String = 'This Slice!';
% %         ent.Callback = 
% %             
% %         bl1 = uicontrol('Parent',f,'Style','text','Position',[50,54,23,23],...
% %                             'String','1','BackgroundColor',bgcolor);
% %         bl2 = uicontrol('Parent',f,'Style','text','Position',[500,54,23,23],...
% %                             'String','280','BackgroundColor',bgcolor);
% %         bl3 = uicontrol('Parent',f,'Style','text','Position',[240,25,100,23],...
% %                             'String','slice number','BackgroundColor',bgcolor);
% %         b.Callback = @(es,ed) imagesc(squeeze(brainMask(round(b.Value),:,:))); 
% %     end
%     
%     
%     while inputz > 0
%        
%         slice = inputz;
%         t = sprintf('Z slice: %d',slice);
%         figure; imagesc(squeeze(brainMask(slice,:,:))); title(t)
%         inputz = input(prompt);
%     end
%     
%     eloc = input('please input [bottom,mid,top] slices of eye (0 for no eye: \n');
%     
%     
%     
%     if (eloc(1) ~=0)
%         figure; imagesc(squeeze(brainMask(eloc(2),:,:))); title('selectEye1')
%         selectEye1 = roipoly;
%         selectEye2 = permute(repmat(selectEye1,1,1,eloc(3)-eloc(1)+1),[3,1,2]);
%         brainMask(eloc(1):eloc(3),:,:) = ~selectEye2.*brainMask(eloc(1):eloc(3),:,:);
%     end
%     if (eloc(1) ~=0)
%         figure; imagesc(squeeze(brainMask(eloc(2),:,:))); title('selectEye2')
%         selectEye1 = roipoly;
%         selectEye2 = permute(repmat(selectEye1,1,1,eloc(3)-eloc(1)+1),[3,1,2]);
%         brainMask(eloc(1):eloc(3),:,:) = ~selectEye2.*brainMask(eloc(1):eloc(3),:,:);
%     end
%     
%     
%     
%     fprintf('Keeping slices from %d and above.\n',slice)
%     temp = 0*brainMask;
%     temp(slice:end,:,:) = 1;
%     brainMask = brainMask.*temp;
%     brainMask = bwareaopen(brainMask,1000);

     volumeViewer(brainMask)
     input('do you like?')
%     
    % interp the brainMask onto the original grid with NN interp
    [x1,x2,x3] = ndgrid(ave_grid_z,ave_grid_y,ave_grid_x);
    F = griddedInterpolant(x1,x2,x3,double(brainMask),'nearest');
    [x1,x2,x3] = ndgrid(grid_Z,grid_Y,grid_X);
    brainMask_interp = F(x1,x2,x3);
    
    
    
    %zMovie(brainMask,250)

    %
    % frankMask is anything but freespace minus the shield?
    frankMask = ~isnan(rho_m);
    frankMask = imfill(frankMask,'holes');
    frankMask = bwareaopen(frankMask,1000);
    
    % interp frankMask
    [x1,x2,x3] = ndgrid(ave_grid_z,ave_grid_y,ave_grid_x);
    F = griddedInterpolant(x1,x2,x3,double(frankMask),'nearest');
    [x1,x2,x3] = ndgrid(grid_Z,grid_Y,grid_X);
    frankMask_i = F(x1,x2,x3);

    %zMovie(frankMask,250);


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
    
    % Recreate the Grid...
    % 
%     [gx,gy,gz] = meshgrid(grid_X,grid_Y,grid_Z);
%     gx = permute(gx,[3,2,1]);
%     gy = permute(gy,[3,2,1]);
%     gz = permute(gz,[3,2,1]);
%     gridTot = zeros(length(grid_Z)-1,length(grid_Y)-1,length(grid_X)-1);
%     for i=1:(length(grid_Z)-1)
%        for j = 1:(length(grid_Y)-1)
%            for k = 1:(length(grid_X)-1)
%                
%            % for every combination
%            % same format as E field, new points for the revoxelized data.
% %            gridx(i,j,k) = ...
% %                (gx(i,j,k) + gx(i,j+1,k)+ gx(i,j,k+1)+ gx(i,j+1,k+1))/4;
% %            gridy(i,j,k) = ...
% %                (gy(i,j,k) + gy(i+1,j,k) + gy(i,j,k+1)+ gy(i+1,j,k+1))/4;
% %            gridz(i,j,k) = ...
% %                (gz(i,j,k) + gz(i+1,j,k) + gz(i,j+1,k)+ gz(i+1,j+1,k))/4;
%                
%             gridTot(i,j,k) = ...
%                (gz(i,j,k) + gz(i,j+1,k)+ gz(i,j,k+1)+ gz(i,j+1,k+1)...
%                 + gy(i,j,k) + gy(i+1,j,k) + gy(i,j,k+1)+ gy(i+1,j,k+1)...
%                 + gx(i,j,k) + gx(i+1,j,k) + gx(i,j+1,k)+ gx(i+1,j+1,k))/12;
%            
%             
%            
%            
%            
%            end
%        end
%     end
%     x_mesh = unique(gridy);
%     y_mesh = unique(gridz);
%     z_mesh = unique(gridx);

    % % % recreate the Sensor Grid
    % % load SolidSensorDef.mat
    % gz = repmat(Z_Dimension_1,[1,length(Y_Dimension_2),length(X_Dimension_3)]);
    % gy = repmat(Y_Dimension_2,[1,length(X_Dimension_3),length(Z_Dimension_1)]);
    % gy = permute(gy,[3,1,2]);
    % gx = repmat(X_Dimension_3,[1,length(Y_Dimension_2),length(Z_Dimension_1)]);
    % gx = permute(gx,[3,2,1]);
    % 
    % for i=1:length(Z_Dimension_1)-1
    %    for j = 1:length(Y_Dimension_2)-1
    %        for k = 1:length(X_Dimension_3)-1
    %            
    %        % for every combination
    %        % same format as E field, new points for the revoxelized data.
    %        gridx(i,j,k) = ...
    %            (gx(i,j,k) + gx(i,j+1,k)+ gx(i,j,k+1)+ gx(i,j+1,k+1))/4;
    %        gridy(i,j,k) = ...
    %            (gy(i,j,k) + gy(i+1,j,k) + gy(i,j,k+1)+ gy(i+1,j,k+1))/4;
    %        gridz(i,j,k) = ...
    %            (gz(i,j,k) + gz(i+1,j,k) + gz(i,j+1,k)+ gz(i+1,j+1,k))/4;
    %             
    %        
    %         
    %        
    %        
    %        
    %        end
    %    end
    % end
    % y = unique(gridy);
    % z = unique(gridz);
    % x = unique(gridx);
    x_mesh = grid_X;
    y_mesh = grid_Y;
    z_mesh = grid_Z;
    close all
    clc
    %brainMask = brainMask_interp;
    
save(fullfile(path,'voxelizedMesh.mat'),'rho_m','sigma_m','epsilonr_m','frankMask','frankMask_i','brainMask','brainMask_interp','x_mesh','y_mesh','z_mesh')