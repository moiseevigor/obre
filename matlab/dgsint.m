function [dgsintresult] = dgsint(fnc, a, b, varargin)
%DGSINT estimate the integral of a user supplied function using 20 pt Gaussian quadrature.
%
%  Returns the integral of a user supplied function via Gaussian Quadrature.
%  Arguments are described as follows;
%
%    fnc   externally supplied double precision function name which returns
%          the function values (this is the function to be integrated).
%          `fnc' is assumed to not blow up at any of the gauss points.
%          (input)
%
%    a     lower integration bound (input)
%
%    b     upper integration bound (input)
%
%    supp  contains a parameter set which is passed to the user-supplied
%          function fnc.  supp may be used to pass any auxiliary 
%          parameters necessary for the computation of the function
%          fnc.  (input)
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

% check the function fnc
f = fcnchk(fnc);

w20 = [0.017614007139152118312d0,0.040601429800386941331d0, ...
       0.062672048334109063570d0,0.083276741576704748725d0, ...
       0.101930119817240435037d0,0.118194531961518417312d0, ... 
       0.131688638449176626898d0,0.142096109318382051329d0, ... 
       0.149172986472603746788d0,0.152753387130725850698d0];
   
z20 = [0.993128599185094924786d0,0.963971927277913791268d0, ...
       0.912234428251325905868d0,0.839116971822218823395d0, ...
       0.746331906460150792614d0,0.636053680726515025453d0, ...
       0.510867001950827098004d0,0.373706088715419560673d0, ...
       0.227785851141645078080d0,0.076526521133497333755d0];

rr = 0.5.*(b - a);
ss = 0.5.*(b + a);

f1  = f(ss+rr.*z20(1), varargin{:});
f2  = f(ss-rr.*z20(1), varargin{:});
f3  = f(ss+rr.*z20(2), varargin{:});
f4  = f(ss-rr.*z20(2), varargin{:});
f5  = f(ss+rr.*z20(3), varargin{:});
f6  = f(ss-rr.*z20(3), varargin{:});
f7  = f(ss+rr.*z20(4), varargin{:});
f8  = f(ss-rr.*z20(4), varargin{:});
f9  = f(ss+rr.*z20(5), varargin{:});
f10 = f(ss-rr.*z20(5), varargin{:});
f11 = f(ss+rr.*z20(6), varargin{:});
f12 = f(ss-rr.*z20(6), varargin{:});
f13 = f(ss+rr.*z20(7), varargin{:});
f14 = f(ss-rr.*z20(7), varargin{:});
f15 = f(ss+rr.*z20(8), varargin{:});
f16 = f(ss-rr.*z20(8), varargin{:});
f17 = f(ss+rr.*z20(9), varargin{:});
f18 = f(ss-rr.*z20(9), varargin{:});
f19 = f(ss+rr.*z20(10), varargin{:});
f20 = f(ss-rr.*z20(10), varargin{:});

dd = w20(1).*(f1+f2) + ... 
     w20(2).*(f3+f4) + ... 
     w20(3).*(f5+f6) + ...
     w20(4).*(f7+f8) + ... 
     w20(5).*(f9+f10) + ...
     w20(6).*(f11+f12) + ...
     w20(7).*(f13+f14) + ... 
     w20(8).*(f15+f16) + ...
     w20(9).*(f17+f18) + ... 
     w20(10).*(f19+f20);

dgsintresult = rr.*dd;

return;
end %END FUNCTION dgsint

