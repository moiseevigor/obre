function [quangpresult] = quangp(thre,alpha,k,p)
%QUANGP Calculates an estimate of the pth quantile of a GP based
%   on the current threshold and estimates of alpha and k.
%
%     thre  value of the threshold used in the analysis
%    alpha  estimate of parameter alpha
%        k  estimate of parameter k
%        p  percentile of the quantile of interest
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

quangpresult = thre + alpha./k.*(1. - p.^k);

return;
end % END FUNCTION quangp

