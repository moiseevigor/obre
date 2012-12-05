function acm = covmat(n,x,alpha,beta,av,we)
%COVMAT To compute the Asymptotic Covariance Matrix (4.2.11, p.231).
%   MQmat computes the M and Q matrices (4.2.13, p.230) needed
%   by CovMat (via averaging over the integral).
%
%      n  the number of values above the threshold
%      x  an array of the data points above the threshold
%  alpha  the current OBRE estimate of alpha
%      beta  the current OBRE estimate of beta
%    acm  the return variable (Asymptotic Covariance Matrix)
%     av  the current a vector
%     we  the array of weights for the array of data points
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

% Evaluate the M and Q matrices
[m,q] = mqmat(n,x,alpha,beta,av,we);

% Compute the inverse of M (Minv) and get (Minv)^t
% ---> Since M is symmetric so is Minv thus: (Minv)^t = Minv
temp = m(1,1)*m(2,2) - m(1,2)*m(2,1);

minv(1,1) = m(2,2)/temp;
minv(2,2) = m(1,1)/temp;
minv(1,2) = -m(1,2)/temp;
minv(2,1) = -m(2,1)/temp;

% Compute MIQ = Minv * Q
miq(1,1) = minv(1,1)*q(1,1) + minv(1,2)*q(2,1);
miq(1,2) = minv(1,1)*q(1,2) + minv(1,2)*q(2,2);
miq(2,1) = minv(2,1)*q(1,1) + minv(2,2)*q(2,1);
miq(2,2) = minv(2,1)*q(1,2) + minv(2,2)*q(2,2);

% Multiply MIQ by M to get ACM
acm(1,1) = miq(1,1)*minv(1,1) + miq(1,2)*minv(2,1);
acm(1,2) = miq(1,1)*minv(1,2) + miq(1,2)*minv(2,2);
acm(2,1) = miq(2,1)*minv(1,1) + miq(2,2)*minv(2,1);
acm(2,2) = miq(2,1)*minv(1,2) + miq(2,2)*minv(2,2);

return;
end % END FUCNTION covmat


function [m,q] = mqmat(n,x,alpha,beta,av,we)
%MQMAT

m = zeros(2,2); q = m;

for i=1:n
    sc1 = score1(alpha, beta, x(i), false);
    sc2 = score2(alpha, beta, x(i), false);
    m(1,1) = m(1,1) + (sc1-av(1))*we(i)*sc1;
    m(1,2) = m(1,2) + (sc1-av(1))*we(i)*sc2;
    m(2,2) = m(2,2) + (sc2-av(2))*we(i)*sc2;
    q(1,1) = q(1,1) + (sc1-av(1))*we(i)*(sc1-av(1))*we(i);
    q(1,2) = q(1,2) + (sc1-av(1))*we(i)*(sc2-av(2))*we(i);
    q(2,2) = q(2,2) + (sc2-av(2))*we(i)*(sc2-av(2))*we(i);
end

m(1,1) = m(1,1)/n;
m(1,2) = m(1,2)/n;
m(2,2) = m(2,2)/n;
m(2,1) = m(1,2);

q(1,1) = q(1,1)/n;
q(1,2) = q(1,2)/n;
q(2,2) = q(2,2)/n;
q(2,1) = q(1,2);

return;
end % END FUCNTION mqmat

