function exportFields(dir,voxelizedMesh,savedir,Nc)
% exportFields:
% purpose: to export B1+ field values from the binary files in the .xf dir
% at specific points in the brain mask and eventually the E-field values 
% as well.
% Input: dir -  (string) the name of the .xf directory, DO NOT INCLUDE THE
%               .xf file notation
%        voxelizedMesh - matfile output of the voxelation function that
%        contains the frankMask
%        savedir - (string) location of where files should be saved
%
% Output: saves the B1+ matrix and Indices file to the savedir

xfdir = sprintf('%s.xf/Simulations/',dir);
x_mesh = voxelizedMesh.x_mesh;
y_mesh = voxelizedMesh.y_mesh;
z_mesh = voxelizedMesh.z_mesh;
brainMask = voxelizedMesh.brainMask;
%oad('voxelizedMesh','x_mesh','y_mesh','z_mesh','brainMask')

% number of rf coil loops

    
    fprintf('Begin exporting B1+ map for %s\n',dir)
    for i=1:Nc
        
        j = i;
        %if dir == '5729_pos'
        %    j = i+1;
        %end
        %if dir == '5688_neg'
        %    j = i+1;
        %end
        % navigate to the correct loop
        sim = sprintf('%06d',j);
        outputDir = fullfile(xfdir,sim,'Run0001','output');
        cd(outputDir)
        fprintf('processing simulation %s\n',sim) 
        
        
        % pre-liminary data
        fileID = fopen('MultiPoint_Solid_Sensor_0_info.bin', 'r', 'ieee-le');
        data.rmpt = fread(fileID,[1 4],'*char');
        data.ver = fread(fileID,1, 'uint8');
        if (data.ver ~= 2)
           fprintf('Warning version is not equal to 2\n'); 
        end
        data.fields = fread(fileID,1, 'uint32');
        % this 7 and 10 indicate that steady state discrete freq, E and B are
        % calculated in the sensor ...
        % b = bitget(data.fields,32:-1:1,'uint32')
        data.num = fread(fileID,1,'uint64');
        fclose(fileID);
        
        % read the geometry file . . . 
        fileID = fopen('MultiPoint_Solid_Sensor_0\geom.bin');
        grid = fread(fileID,Inf,'uint32');
        fclose(fileID);
        %
        if ( i~=1 )
           t = [xpts(1),ypts(1),zpts(1)]; 
        end 
        x_d = grid(1:3:end,1);
        y_d = grid(2:3:end,1);
        z_d = grid(3:3:end,1);
        xpts = unique(x_d,'stable');
        ypts = unique(y_d,'stable');
        zpts = unique(z_d,'stable');
        
        dim = [numel(zpts),numel(ypts),numel(xpts)];
        if (i==1)
           masterDim = [numel(zpts),numel(ypts),numel(xpts)];
        else
           fprintf('Old:%d,%d,%d', length(xpts),length(ypts),length(zpts))
           xpts = xpts(1:masterDim(3));
           ypts = ypts(1:masterDim(2));
           zpts = zpts(1:masterDim(1));
           fprintf('- New:%d,%d,%d\n',length(xpts),length(ypts),length(zpts))
        end
        
        
        if ( i~=1)
            if (t(1) ~= xpts(1) || t(2) ~= ypts(1) || t(3) ~= zpts(1))
               fprintf('Warning, sensor starting points are different.\n')
               
            end
            
        end
        
        % ensure the order of change is correct
        xc = find(diff(x_d),1);
        yc = find(diff(y_d),1);
        zc = find(diff(z_d),1);
        
        % Ensure the direction of change of the dimensions is zyx
        if (xc < yc || yc < zc || xc < zc)
            fprintf('Warning, change of indices is unexpected...\n')
        end
        % ensure the number of points is the same as expected
        if (data.num ~= (dim(1)*dim(2)*dim(3)))
           fprintf('Warning number of points saved ~= number expected\n')
           
        end
        % ensure the pts are increasing
        if (sum(diff(xpts))~=numel(xpts)-1 || sum(diff(ypts))~=numel(ypts)-1 || sum(diff(zpts))~=numel(zpts)-1)
            fprintf('Warning number of points are not increasing\n')
            
        end
        %dim = [numel(zpts),numel(ypts),numel(xpts)];
       

        % frequency
        fileID = fopen('MultiPoint_Solid_Sensor_0\frequencies.bin');
        freq = fread(fileID,1,'float32');
        fclose(fileID);
        
        % B-Field x part
        fileID = fopen('MultiPoint_Solid_Sensor_0\ss_Bxrt\0.bin');
        Bx_r = fread(fileID,Inf,'float32');
        fclose(fileID);
        fileID = fopen('MultiPoint_Solid_Sensor_0\ss_Bxit\0.bin');
        Bx_i = fread(fileID,Inf,'float32');
        fclose(fileID);
        Bx_Tot = reshape(Bx_r,dim) + 1i*reshape(Bx_i,dim);
        Bx_Tot = Bx_Tot(1:masterDim(1),1:masterDim(2),1:masterDim(3));
        
        % B-Field y part
        fileID = fopen('MultiPoint_Solid_Sensor_0\ss_Byrt\0.bin');
        By_r = fread(fileID,Inf,'float32');
        fclose(fileID);
        fileID = fopen('MultiPoint_Solid_Sensor_0\ss_Byit\0.bin');
        By_i = fread(fileID,Inf,'float32');
        fclose(fileID);
        By_Tot = reshape(By_r,dim) + 1i*reshape(By_i,dim);
        By_Tot = By_Tot(1:masterDim(1),1:masterDim(2),1:masterDim(3));
        
        
        nonZB_v = find(brainMask(zpts+1,ypts+1,xpts+1));
        
        B1plus_m(:,i) = 1e6*(Bx_Tot(nonZB_v)+1i*By_Tot(nonZB_v))/2;
        fprintf('Coil %d loaded and placed in B1plus_m for %s...\n', i,dir)
        
        % cd('F:/CompletedFields')
    
    end

    % save the file
    saveto = sprintf('%s/B1plus_m.mat',savedir);
    save(saveto,'B1plus_m','-v7.3')
    indexMesh_z = zpts+1;
    indexMesh_y = ypts+1;
    indexMesh_x = xpts+1;
    saveto = sprintf('%s/Indices.mat',savedir);
    save(saveto,'indexMesh_x','indexMesh_y','indexMesh_z');


end