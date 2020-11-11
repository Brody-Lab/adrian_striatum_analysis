function [fits1,fits2] = match_up_fits(fits1,fits2)
    % takes a table of glm fits across cells and matches the indices so
    % that fits1(1,:) corresponds to the same cell as fits2(i,:)
    idx1 = fits1(:,{'cellno','sess_date','rat'});
    idx2 = fits2(:,{'cellno','sess_date','rat'});
    [~,i1,i2] = intersect(idx1,idx2);
    fits1 = fits1(i1,:);
    fits2 = fits2(i2,:);
    fprintf('Found %g matching rows.\n',height(fits1));
end
    
