function done = compare(avec,oldav,amat,oldam)
%COMPARE Compares new values of a and A with previously calculated values to see if we have converged to a solution
%
%        p   dimension of a (px1) and A (pxp)
%       avec   current a vector
%       amat   current A matrix
%        s   logical value indicating if solved satisfactorily
%      oldav   a vector from previous step
%      oldam   A matrix from previous step
%       df   definition of 'significantly different'
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

% Maximum relative difference (rdf) is 5%
% Maximum absolute difference (df) used in cases of division by zero
rdf = .05d0; df = .1d0;
done = true;

% Check 'a' vector
% No division by zero expected if oldav(1) or oldav(2) not equal to zero
if(oldav(1) ~= 0 && oldav(2) ~= 0)
    if(abs((avec(1)-oldav(1))./oldav(1)) > rdf ||abs((avec(2)-oldav(2))./oldav(2)) > rdf)
        done = false;
    end
else
    if(abs(avec(1)-oldav(1)) > df || abs(avec(2)-oldav(2)) > df)
        done = false;
    end
end

% Check 'A' matrix
if(oldam(1,1) ~= 0 && oldam(1,2) ~= 0 && oldam(2,2) ~= 0)
    if(abs((amat(1,1)-oldam(1,1))./oldam(1,1)) > rdf ||abs((amat(1,2)-oldam(1,2))./oldam(1,2)) > rdf ||abs((amat(2,2)-oldam(2,2))./oldam(2,2)) > rdf )
        done = false;
    end
    
else
    if(abs(amat(1,1)-oldam(1,1)) > df || abs(amat(1,2)-oldam(1,2)) > df ||abs(amat(2,2)-oldam(2,2)) > df)
        done = false;
    end
end

return;
end % END FUNCTION compare

