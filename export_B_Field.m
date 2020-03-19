function export_B_Field(dir,SimNumber,savedir)
% export_B_Field:
% purpose: to export B1+ field values from the binary files in the .xf dir
%
% Input: dir -  (string) the name of the .xf directory, INCLUDE THE
%               .xf file notation
%        SimNumber - integer number representing the simulation number of
%        interest
%        savedir - (string) location of where files should be saved
%
% Output: saves the B1+ matrix and Indices file to the savedir
h = cd;
sim = sprintf('%06d',SimNumber);
% may need to add Simulations back into this or take it out?
xfdir = sprintf('%s/Simulations/%s',dir,sim);

% fprintf('Begin exporting B1+ map for %s\n',xfdir)


outputDir = fullfile(xfdir,'Run0001','output');
cd(outputDir)

% pre-liminary data
fileID = fopen('MultiPoint_Solid_Sensor_for_strip_line_0_info.bin', 'r', 'ieee-le');
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
fileID = fopen('MultiPoint_Solid_Sensor_for_strip_line_0\geom.bin');
grid = fread(fileID,Inf,'uint32');
fclose(fileID);
% ...

% the index of the grid
x_d = grid(1:3:end,1);
y_d = grid(2:3:end,1);
z_d = grid(3:3:end,1);
xpts = unique(x_d,'stable');
ypts = unique(y_d,'stable');
zpts = unique(z_d,'stable');
dim = [numel(zpts),numel(ypts),numel(xpts)];
% ...

% ensure the grid is ascending and ordered in the typical way...
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
% ...

% frequency
fileID = fopen('MultiPoint_Solid_Sensor_for_strip_line_0\frequencies.bin');
freq = fread(fileID,1,'float32');
fclose(fileID);

% B-Field x part
fileID = fopen('MultiPoint_Solid_Sensor_for_strip_line_0\ss_Bxrt\0.bin');
Bx_r = fread(fileID,Inf,'float32');
fclose(fileID);
fileID = fopen('MultiPoint_Solid_Sensor_for_strip_line_0\ss_Bxit\0.bin');
Bx_i = fread(fileID,Inf,'float32');
fclose(fileID);
Bx_Tot = reshape(Bx_r,dim) + 1i*reshape(Bx_i,dim);

% B-Field y part
fileID = fopen('MultiPoint_Solid_Sensor_for_strip_line_0\ss_Byrt\0.bin');
By_r = fread(fileID,Inf,'float32');
fclose(fileID);
fileID = fopen('MultiPoint_Solid_Sensor_for_strip_line_0\ss_Byit\0.bin');
By_i = fread(fileID,Inf,'float32');
fclose(fileID);
By_Tot = reshape(By_r,dim) + 1i*reshape(By_i,dim);

B1plus_m = 1e6*(Bx_Tot+1i*By_Tot)/2;

% save the file
if ~exist(sprintf('%s/B1plus_folder',savedir),'dir')
   mkdir(sprintf('%s/B1plus_folder',savedir))
end
saveto = sprintf('%s/B1plus_folder/B1plus_%d.mat',savedir,SimNumber);
save(saveto,'B1plus_m','xpts','ypts','zpts','-v7.3')
cd(h)


fprintf('B1+ map exported for %d\n',SimNumber)

end