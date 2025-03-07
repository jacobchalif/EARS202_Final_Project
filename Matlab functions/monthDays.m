function months = monthDays(LY)
    % monthDays  returns vector with number of days in each month
    %   months = monthDays(LY)
    %   LY should be boolean - TRUE if leap year, otherwise FALSE
    %
    % also see isLeap
    if LY
        months = [31,29,31,30,31,30,31,31,30,31,30,31]';  % leap year days
    else
        months = [31,28,31,30,31,30,31,31,30,31,30,31]';  % non-leap year days
    end
end