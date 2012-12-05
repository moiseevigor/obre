function [thre,x,n,iouts,p,finala,finalk,we]=goodoffit(thre,x,n,iouts,p,finala,finalk,we);
%GOODOFFIT

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

if isempty(i), i=0; end;
if isempty(header), header=false; end;
if isempty(prognm), prognm=repmat(' ',1,32); end;
if isempty(z), z=zeros(1,n); end;
if isempty(alphahat), alphahat=0; end;
if isempty(khat), khat=0; end;
if isempty(w2), w2=0; end;
if isempty(a2), a2=0; end;
if isempty(al), al=0; end;
if isempty(be), be=0; end;
if isempty(dz), dz=zeros(1,n); end;
if isempty(finalz), finalz=zeros(1,n); end;
if firstCall,   header=[true];  end;
if firstCall,   prognm =['ThisisGoodofFit(Version2.0):'];  end;
firstCall=0;
% *** Formatting statements ***
%format(a,a,a);
%format(a20,i9);
%format(a20,f10.5);
%format(a30,f10.5);
% *** Print introduction ***
if( header )
    writef(iouts,['%s','%s','%s' ' \n'],,' ', prognm);
    writef(iouts,['%0.15g \n']);
    writef(iouts,['%30s','%10.5f' ' \n'],'Starting Value of thre: ',thre);
    writef(iouts,['%0.15g \n']);
    writef(iouts,['%s','%s','%s' ' \n'],'*****************************************');
    header = false;
end;
%************************** Main Body of Program ***************************
% *** Compute parameter estimates using PWM as starting values.  The PWM ***
% ***  estimates are fed to MLE and the MLE estimates are fed to gof.    ***
[n,x,al,be]=meth2(n,x,al,be);
[x,n,al,be]=meth3(x,n,al,be,(al./be));
alphahat = (al);
khat     = (be);

% *** Cramer-von Mises Statistics ***
w2 = 0.0;
for  i=1:n;
    z(i) = 1.0 -( 1.0 - khat.*x(i)./alphahat ).^(1.0./khat);
    finalz(i)=1.0d0-(1.0d0-finalk.*(x(i))./finala).^(1.0d0./finalk);
    w2 = w2 +(z(i) -(2.0.*(i)-1.0)./(2.0.*(n))).^2;
end;  i=fix(n+1);
w2 = w2 + 1.0./(12.0.*(n));

% *** Anderson-Darling Statistic ***
for  i=1:n;
    dz(i) = (z(i));
end;  i=fix(n+1);
a2 = (dad(dz,n,p));

% *** Output results ***
writef(iouts,['%0.15g \n']);
writef(iouts,['%0.15g \n']);
writef(iouts,['%20s','%10.5f' ' \n'],'Value of thre: ',thre);
writef(iouts,['%20s','%9i' ' \n'],'n    :',n);
writef(iouts,['%20s','%10.5f' ' \n'],'alphahat: ',alphahat);
writef(iouts,['%20s','%10.5f' ' \n'],'khat: ',khat);
writef(iouts,['%20s','%10.5f' ' \n'],'Cramer-von Mises =',w2);
writef(iouts,['%20s','%10.5f' ' \n'],'Anderson-Darling =',a2);
writef(iouts,['%20s','%10.5f' ' \n'],'OBRE A-D statistic =', dad(finalz,n,p));
writef(iouts,['%20s','%10.5f' ' \n'],'OBRE A-D weighted  =', dadwght(finalz,n,p,we));
return;
end %subroutine goodoffit
