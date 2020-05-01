function benView(matrix)
% I'm always looking at a matrix somewhere in the middle to get a feel for
% the data so I created this function.
% Benjamin M Hardy 4-21-2020

middle = round(size(matrix)/2);
try
    figure;
    imagesc(squeeze(matrix(:,:,middle(3))))
catch
    sprintf('matrix is complex')
    figure;
    imagesc(squeeze(abs(matrix(:,:,middle(3)))))
end



end