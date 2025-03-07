function leap = isLeap(year)
    % isLeap  returns TRUE if year is a leap year, otherwise FALSE
    %   leap = isLeap(year)

    if mod(year,4) ~= 0
        leap = false;
    elseif mod(year,100) ~= 0
        leap = true;
    elseif mod(year,400) ~= 0
        leap = false;
    else
        leap = true;
    end
end

