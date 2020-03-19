function names_array = getHeadNames()
    % label_%s_%d_%s.raw',vector,a-1,simIDs{i}),'w+'
	simIDs = {'5634','5639','5652','5657','5662','5671','5674','5688','5693','5698','5708','5719','5729','5757','5762','5767','5772','5780','5785','5790'};
    names_array = {};
    for i=1:numel(simIDs)
	% s[0] = 'label_100_0_5634';
		% 100
		for j=1:50
            names_array{end+1} = sprintf('label_100_%d_%s', j-1, simIDs{i});	
        end
		% 010
		for j=1:50
            names_array{end+1} = sprintf('label_010_%d_%s',j-1,simIDs{i});	
        end
		% 001
		for j=1:50
            names_array{end+1} = sprintf('label_001_%d_%s',j-1,simIDs{i});
        end
    end
end



        