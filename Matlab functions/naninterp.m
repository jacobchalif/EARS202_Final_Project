function X = naninterp(X)
    % naninterp  interpolates over NaNs
    %   X = naninterp(X)
    %   See INTERP1 for more info
    
    X(isnan(X)) = interp1(find(~isnan(X)), X(~isnan(X)), find(isnan(X)),'cubic');
return