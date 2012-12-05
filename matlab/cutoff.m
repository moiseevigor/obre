function ret = cutoff(x, minval, maxval)
%CUTOFF the function cuts off the values out of interval
%

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

if nargin<2
  error('MATLAB:cutoff','Not enough input arguments.'); 
end

if nargin == 2
    x(x<minval) = minval;
end

if nargin == 3
    x(x>maxval) = maxval;
    x(x<minval) = minval;
end

ret = x;

return;