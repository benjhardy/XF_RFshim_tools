function renumber_sims(xfdir)
    % takes a group of abnormally ordered simulations and changes the file names to linear ordered sims
    % Input: the xfdir full path in a string, including the '.xf'
    % Output: the simulations folder reflects -> 000001, 000002, 000003 ...
    h = pwd;

    % Jan's Code: 
    % https://www.mathworks.com/matlabcentral/answers/16283-renaming-a-lot-of-folders-automatically-by-matlab
    PathName = [xfdir,'/Simulations'];
    ADir = dir(PathName);
    ADir = ADir([ADir.isdir]);
    AName = {ADir.name};
    AName(strncmp(AName, '.', 1)) = [];
    
    % check if the simulations need renamed or not
    simnumbers = arrayfun(@(x) str2double(x), AName, 'UniformOutput',false);
    simnumberss = cell2mat(simnumbers);
    if (simnumberss(1) == 1 && sum(diff(simnumberss)) == numel(AName) - 1)
       fprintf('Simulations folders are already in order...\n')
       return 
    end
    
    for iFolder = 1:numel(AName)
        
        try
            APath   = fullfile(PathName, AName{iFolder});
            % BDir    = dir(fullfile(APath, '*'));
            newName = sprintf('%06d',iFolder);
            movefile(APath, fullfile(PathName, newName));
        catch
            fprintf('%s already exists and is in correct location.\n',newName);
        end
        fprintf('%.2f %% complete...\n',100*iFolder/numel(AName));
    end
    
    cd(h)




end