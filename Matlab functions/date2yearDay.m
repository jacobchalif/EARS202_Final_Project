function ind = date2yearDay(year,month,day)
    % date2yearDay  returns index of date in a given year
    %   ind = date2yearDay(year,month,day)
    %
    %   IF LEAP YEAR: includes Feb 29.
    %       For example, date2yearDay(2020,12,31) will return 366
    %                but date2yearDay(2019,12,31) will return 365


    days = monthDays(isLeap(year)); % stores number of days in each month
    ind = day; % adds days passed in current month to index;
    ind = ind + sum(days(1:month-1)); % adds number of days that have passed before current month
end
