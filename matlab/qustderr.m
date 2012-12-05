function qustderrresult = qustderr(alpha, k, p, acm, n)
%QUSTDERR Calculate an estimate of the standard error of the quantile 
%	estimate based on the Delta method.
%
%    alpha  OBRE of parameter alpha
%        k  OBRE of parameter k
%        p  percentile of the quantile of interest
%      acm  asymptotic covariance matrix
%        n  sample size
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

% First-order partials of the quantile estimate
% dqda = dq_p/dalpha
dqda =(1.0d0 - p.^k)./k;

% dqdk = dq_p/dk
dqdk = -alpha*(1.0d0 - p.^k + p.^k*log(p)*k)./(k.^2);

qustderrresult = sqrt( (dqda.^2*acm(1,1) + dqdk.^2*acm(2,2) +2.0d0*dqda*dqdk*acm(1,2))/n );

return;
end % END FUNCTION qustderr

