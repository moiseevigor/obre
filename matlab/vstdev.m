function [vstdevresult] = vstdev(xs,n)
% VSTDDEV Calculate the standard deviation of an array of data
%
%       n   size of the array
%       xs  array of data (Stored in Single Precision)
%
%OBRE project, 2011.

%    OBRE for Matlab
%    Copyright (C) 2011, Igor Moiseev
%
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program; if not, see http://www.gnu.org/licenses
%    or write to the Free Software Foundation,Inc., 51 Franklin Street,
%    Fifth Floor, Boston, MA 02110-1301  USA

mean = vmean(xs,n);
summlv = 0;

for  i=1:n;
    summlv = summlv +((xs(i))-mean).*((xs(i))-mean);
end;

vstdevresult = sqrt(summlv./(n-1));

return;
end % END FUNCTION vstdev

