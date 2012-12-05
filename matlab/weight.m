function [weightresult] = weight(x, alpha, beta, avec, amat)
%WIEGHT returns the weight for a given c, am, score, and av
%
%     x  real value at which the function is to be evaluated
% alpha  the alpha parameter of GP
%  beta  the beta parameter of GP
%  amat  A matrix
%  avec  a vector
% dat.c  bound on the IF
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


global dat;

[sco1,feas1] = score1(alpha, beta, x, false);
[sco2,feas2] = score1(alpha, beta, x, false);

diff(1) = sco1-avec(1);
diff(2) = sco2-avec(2);

vector(1) = amat(1,1).*diff(1)+amat(1,2).*diff(2);
vector(2) = amat(2,1).*diff(1)+amat(2,2).*diff(2);

norm = sqrt(vector(1).^2+vector(2).^2);

weightresult = min( 1.0d0, dat.c./norm );

if((~feas1)||(~feas2))
    weightresult = 0.0d0;
end

return;
end % END FUNCTION weight

