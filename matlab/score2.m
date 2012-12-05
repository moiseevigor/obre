function [score2result, feas] = score2(alpha, beta, x, feas)
%SCORE1 is to return the second element (score2) of the score function corresponding to a GP distribution
%
%  alpha alpha value
%   beta  beta value
%      x  real value at which the function is to be evaluated
%   feas  boolean, the result of successfull execution
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

if( (beta > 0 && x >(alpha/beta)) || x < 0 )
    score2result = 0;
    feas_sco2 = false;
else
    score2result =((alpha-beta.*x).*log(alpha-beta.*x) +(beta.*x-alpha).*log(alpha) + beta.*x.*(1-beta) )./((beta.^2).*(beta.*x-alpha) );
    feas_sco2 = true;
end

feas = feas_sco2;

return;
end % END FUNCTION score2

