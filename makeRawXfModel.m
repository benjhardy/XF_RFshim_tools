function makeRawXfModel(label_m,list,list_nums,dir,name)
% 
% Takes a label matrix and a list of labels to create an xf
% label file that can be imported into an XF project using the xf gui...
% the following code is based on already existing scripts that can be found
% at the following locations:
% E:/CoordinateSystemConfirmation/batch_1/get_heads.m
%
%
% ---------- > 4-21-2020 ===== Benjamin M. Hardy < -----------

% determine orientation
% may need to set the label_m file into XF space orientation
% set into XF space
% label_1 = permute(label_m,[2,1,3]);

% pass the nifti

% Create the RAW file!
raw_name = sprintf('%s.raw',name);
fid = fopen(fullfile(dir,raw_name),'w+');
fwrite(fid,label_m,'uint8');
fclose(fid);
        
if (class(list) ~= 'cell')
    % the text file
    text = sprintf(['1	0.411765	0	0	CSF General\n',...
    '2	0.521569	0.521569	0.533333	Brain Gray Matter\n',...
    '3	1	1	1	Brain White Matter\n',...
    '4	1	1	0.588235	Skull\n',...
    '8	1	0	0	Fat\n',...
    '10	1	0.72549	0.564706	skin\n',...
    '20	0	0	1	Electrode\n',...
    '0	0	0	0	Background\n',...
    '\n',...
    '\n',...
    'Grid extent (number of cells)\n',...
    'nx	%d\n',...
    'ny	%d\n',...
    'nz	%d\n',...
    '\n',...
    'Spatial steps [m]\n',...
    'dx	0.001\n',...
    'dy	0.001\n',...
    'dz	0.001'],size(label_m,1),size(label_m,2),size(label_m,3));
    text_name = sprintf('%s.txt',name);
    path2file = fullfile(dir,text_name);
    fid = fopen(path2file,'w+');
    fwrite(fid,text);
    fclose(fid);
end

if (class(list) == 'cell')
    % the text file
    
    % somewhat random color distributions
    r_v = linspace(.01,.75,numel(list));
    g_v = flip(r_v);
    b_v = [r_v(1:2:end),r_v(2:2:end)];
    
    text = '';
    % add the labels to the text file
    for i=1:numel(list)
       key = sprintf(['%d\t%0.5f\t%0.5f\t%0.5f\t%s\n'],list_nums(i),r_v(i),g_v(i),b_v(i),list{i});
       text = sprintf('%s%s',text,key);
    end
    
    % add the tail end of the text file
    tail = sprintf(['\n',...
    '\n',...
    'Grid extent (number of cells)\n',...
    'nx	%d\n',...
    'ny	%d\n',...
    'nz	%d\n',...
    '\n',...
    'Spatial steps [m]\n',...
    'dx	0.001\n',...
    'dy	0.001\n',...
    'dz	0.001'],size(label_m,1),size(label_m,2),size(label_m,3));
    text = strcat(text,tail);

    text_name = sprintf('%s.txt',name);
    path2file = fullfile(dir,text_name);
    fid = fopen(path2file,'w+');
    fwrite(fid,text);
    fclose(fid);
end











end