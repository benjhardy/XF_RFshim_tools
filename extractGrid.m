function [grid_Xr,grid_Yr,grid_Zr] = extractGrid(dir,sim)
% input: dir, including the .xf and sim number to extract grid from
% read the geometry file
% contains everything like material definitions and numbers and values of
% permittivity and everything
% for now I just want the specific grid indices, later will adjust to get
% the entire makeup of the mesh

ff = fullfile(dir,'Simulations',sprintf('%06d',sim),'Run0001','geometry.input');
% ff = fullfile(dir,sprintf('%06d',sim),'Run0001','geometry.input');
fid = fopen(ff,'r');
text = fscanf(fid,'%c');
fclose(fid);

% find the number of Cells in Y and X a
Nxs = strfind(text,'NumberOfCellsInX');
Nys = strfind(text,'NumberOfCellsInY');
Nzs = strfind(text,'NumberOfCellsInZ');
en = strfind(text,'end_<GridDefinition>');

% find origin
os = strfind(text,'GridOriginInMeters');
O = text(os+18:Nxs-1);
origin = str2num(O);

Nx = str2double(text(Nxs+16:Nys-1));
Ny = str2double(text(Nys+16:Nzs-1));
Nz = str2double(text(Nzs+16:en-1));

xi = strfind(text,'begin_<DelX>');
xf = strfind(text,'end_<DelX>');
yi = strfind(text,'begin_<DelY>');
yf = strfind(text,'end_<DelY>');
zi = strfind(text,'begin_<DelZ>');
zf = strfind(text,'end_<DelZ>');

% might need to edit the precision of these numbers
x_step = str2num(text(xi+13:xf-1));
y_step = str2num(text(yi+13:yf-1));
z_step = str2num(text(zi+13:zf-1));

% confirm these were read correctly
if (sum(isnan([x_step(:);y_step(:);z_step(:);])) > 0 )
    fprintf('Warning, step sizes were not read correctly\n')
    fprintf('%s...\n',cd)
    %return
end



% make step array
% X
step = zeros(Nx,1);
x_step = cat(1,x_step,[Nx,x_step(size(x_step,1),2)]);
step(1:x_step(2,1)+1) = x_step(1,2);
for i=2:size(x_step,1)
   step(x_step(i,1)+1:x_step(i+1)) = x_step(i,2);
end

x = zeros(Nx,1);
x(1) = origin(1);
for i = 1:Nx-1
    x(i+1) = x(i) + step(i);
end
x = 1000*x;

% Y
step = zeros(Ny,1);
y_step = cat(1,y_step,[Ny,y_step(size(y_step,1),2)]);
step(1:y_step(2,1)+1) = y_step(1,2);
for i=2:size(y_step,1)-1
   step(y_step(i,1)+1:y_step(i+1)) = y_step(i,2);
end

y = zeros(Ny,1);
y(1) = origin(2);
for i = 1:Ny-1
    y(i+1) = y(i) + step(i);
end
y = 1000*y;

% Z
step = zeros(Nz,1);
z_step = cat(1,z_step,[Nz,z_step(size(z_step,1),2)]);
step(1:z_step(2,1)+1) = z_step(1,2);
for i=2:size(z_step,1)-1
   step(z_step(i,1)+1:z_step(i+1)) = z_step(i,2);
end

z = zeros(Nz,1);
z(1) = origin(3);
for i = 1:Nz-1
    z(i+1) = z(i) + step(i);
end
z = 1000*z;

grid_Xr = x;
grid_Yr = y;
grid_Zr = z;

