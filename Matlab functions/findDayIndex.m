function ind = findDayIndex(year,month,day,startYear)
    % findDayIndex  returns index of date given any start year
    %   ind = findDayIndex(year,month,day,startYear)
    %
    %   Includes leap days.

    if year < startYear
        error("year must be greater than or equal to startYear")
    end
    
    yearsCompleted = (startYear:year-1)'; % total number of completed years
    completedLeaps = sum(arrayfun(@(x) isLeap(x),yearsCompleted)); % number of completed leap years
    completedNonLeaps = (year - startYear) - completedLeaps; % number of completed non-leap years
    
    ind = (366*completedLeaps) + (365*completedNonLeaps) + date2yearDay(year,month,day);
end