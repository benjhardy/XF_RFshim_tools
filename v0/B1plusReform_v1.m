function [B1plus_m, FOX] = B1plusReform(B1plus_m,ver,ix,iy,iz)
    % Purpose: Restructures the B1plus_m array into an array where only 
    % the brain is included
    %
    % Input: the original B1plus_m
    %        voxelizedMesh -  struct containing the brain and frankMask
    %        Indices - the indices of mesh where the sensor is located
    %        flag to indicate whether or not the operation should be done
    %        or not, may be redundant, yet feels like a safety feature to
    %        me
    % Output: the B1plus_m only present for the brainMask
    
    % first check to see if the B1plus_m has already been restructured
    if ver == 1
        load('Indices.mat','indexMesh_x','indexMesh_y','indexMesh_z')
        load('voxelizedMesh','frankMask')
        sensorfrankMask = frankMask(indexMesh_z,indexMesh_y,indexMesh_x);
        load('voxelizedMesh', 'brainMask');
        FOX = logical(brainMask(indexMesh_z,indexMesh_y,indexMesh_x));
        ptsfm = length(nonzeros(sensorfrankMask));
        pts = size(B1plus_m,1);
        flag = true;
    
            if (pts ~= ptsfm)
            fprintf('B1plus_m has already been restructed and has length %d\n',pts) 
            flag = false;
            end
    
    
        if (flag == true)
            fprintf('Restructuring B1plus_m...\n') 
            temp = sensorfrankMask + FOX;
            temp = nonzeros(temp(:));
            brainIn = find(temp == 2);
            B1plus_m = B1plus_m(brainIn,:);
            fprintf('B1plus_m has been restructured...\n') 
        end
    end
    if ver == 2
        load('voxelizedMesh.mat','frankMask_i')
        sensorfrankMask = frankMask_i(iz,iy,ix);
        load('voxelizedMesh.mat', 'brainMask_interp');
        FOX = logical(brainMask_interp(iz,iy,ix));
        fprintf('Restructuring B1plus_m...\n') 
        temp = (sensorfrankMask + FOX)-1;
        % temp = nonzeros(temp(:));
        brainIn = temp > 0;
        B1plus_m = B1plus_m(brainIn,:);
        fprintf('B1plus_m has been restructured...\n')
    end
    
%     ptsfm = length(nonzeros(sensorfrankMask));
%     pts = size(B1plus_m,1);
%     flag = true;
%     
%     %if (pts ~= ptsfm)
%     %   fprintf('B1plus_m has already been restructed and has length %d\n',pts) 
%     %   flag = false;
%     %end
%     
%     if (flag == true)
%         fprintf('Restructuring B1plus_m...\n') 
%         temp = (sensorfrankMask + FOX)-1;
%         % temp = nonzeros(temp(:));
%         brainIn = temp > 0;
%         B1plus_m = B1plus_m(brainIn,:);
%         fprintf('B1plus_m has been restructured...\n') 
%     end
    
end