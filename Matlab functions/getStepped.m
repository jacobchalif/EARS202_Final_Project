function [stepX,stepY] = getStepped(minX,maxX,y)
    % getStepped  returns stepped x and y values for y data with top and
    %             bottom values. Useful for discrete ice core data.
    %   
    %   [stepX,stepY] = getStepped(minX,maxX,y)
    %

    [~,si] = sort(minX,'ascend');
    minX = minX(si);
    maxX = maxX(si);
    y = y(si);
    N = length(y);
    stepX = nan(N*2,1);
    stepY = nan(N*2,1);
    for i = 1:N
        stepX((i-1)*2+1) = minX(i);
        stepX((i-1)*2+2) = maxX(i);
        stepY((i-1)*2+1) = y(i);
        stepY((i-1)*2+2) = y(i);
    end
end
