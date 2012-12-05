function j = fishinfo(par1,par2,y,n)
%FISHINFO Calculates the fisher information matrix corresponding to
%    a vector of n observations from a distribution with
%    para(2) parameters whose scores are stored in score1 and score2.
%    This essentially performs step two of the OBRE fitting algorithm.
%
%    Outputs:
%        j   fisher information matrix
%
%    Inputs:
%     par1   alpha value
%     par2   beta value
%        y   vector containing the simulated observations (OUTPUT)
%        n   number of observations (INPUT)
%
%OBRE project, 2011.

% Ver. 0.1 by Moiseev Igor, 27/01/2011
% Converted to Matlab from fishinfo.f

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

% initializations
j = zeros(2,2);
feas1=false; feas2=false;

% Have replaced dH with its empirical distribution function.
% This means replacing the integral with averages over the sample.

for  i=1:n
    [sco1,feas1] = score1(par1,par2,y(i),feas1);
    [sco2,feas2] = score2(par1,par2,y(i),feas2);

    j(1,1) = j(1,1) + sco1.*sco1;
    j(1,2) = j(1,2) + sco1.*sco2;
    j(2,2) = j(2,2) + sco2.*sco2;
end

j(1,1)=j(1,1)./n;
j(1,2)=j(1,2)./n;
j(2,1)=j(1,2);
j(2,2)=j(2,2)./n;

return;
end %END FUNCTION fishinfo

