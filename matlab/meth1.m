function [afin,bfin] = meth1(n,p)
% METH1 Estimate the parameters of a GP distribution via the Method of Moments
%
%        n  number of data points
%        p  array of data
%     afin  final alpha value
%     bfin  final beta value 
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

mean = vmean(p,n);
var  = vstdev(p,n).^2;
afin = 0.5d0.*mean.*((mean.*mean)./var + 1.0d0);
bfin = 0.5d0.*((mean.*mean)./var - 1.0d0);

return;
end % subroutine meth1
