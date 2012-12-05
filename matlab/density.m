function [densityresult] = density(x, alpha, beta)
%DENSITY returns the density for a given x corresponding to a PWM distribution
%
%     x  real value at which the function is to be evaluated
% alpha  alpha parameter of GP
%  beta  beta  paremeter of GP
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

if(beta == 0)
    densityresult = 1.0/alpha*exp(-x/alpha);
elseif( (beta > 0 && x > alpha/beta) || x < 0 ) ;
    densityresult = 0.0;
else
    densityresult = 1.0/alpha*(1.0 - beta*x/alpha)^(1.0/beta - 1.0);
end

return;
end % END FUNCTION density

