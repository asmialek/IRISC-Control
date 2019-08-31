function [ dayJ2000 ] = date2J2000( year, month, day, hours, minutes, seconds )
%date2J2000 date conversion
%   converts date from format yy, mm, dd, time ot J2000 format (days from
%   1.1.2000 12:00 incl. time)

% days_month_normal = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334];
% days_month_leap = [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335];
% days_year_list = [2018, 2019, 2020; 6573.5, 6938.5, 7303.5];
% 
% %check if year is leap year
% if year == 2020
%     index = find(days_year_list(1,:)==year);
%     days_month = days_month_leap;
% else
%     days_month = days_month_normal;
% end
% 
% index = find(days_year_list(1,:)==year);
% days_year = days_year_list(2,index);

dayJ2000 = juliandate(year, month, day, hours, minutes, seconds) - juliandate(2000, 1, 1, 12, 0, 0);

end

