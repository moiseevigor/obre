function [d,v]=jacobi(a,n,np,d,v,nrot);
%JACOBI
% Automatic coinversion from Fortran
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

persistent b c g h i ip iq j nmax s sm t tau theta tresh z ;

a_orig=a;a_shape=[n,n];a=reshape([a_orig(1:min(prod(a_shape),numel(a_orig))),zeros(1,max(0,prod(a_shape)-numel(a_orig)))],a_shape);
v_orig=v;v_shape=[n,n];v=reshape([v_orig(1:min(prod(v_shape),numel(v_orig))),zeros(1,max(0,prod(v_shape)-numel(v_orig)))],v_shape);

if isempty(nmax), nmax=500; end;
if isempty(i), i=0; end;
if isempty(ip), ip=0; end;
if isempty(iq), iq=0; end;
if isempty(j), j=0; end;
if isempty(c), c=0; end;
if isempty(g), g=0; end;
if isempty(h), h=0; end;
if isempty(s), s=0; end;
if isempty(sm), sm=0; end;
if isempty(t), t=0; end;
if isempty(tau), tau=0; end;
if isempty(theta), theta=0; end;
if isempty(tresh), tresh=0; end;
if isempty(b), b=zeros(1,nmax); end;
if isempty(z), z=zeros(1,nmax); end;

for  ip=1:n;
    for  iq=1:n;
        v(ip,iq)=0.0d0;
    end;  iq=fix(n+1);
    v(ip,ip)=1.0d0;
end;  ip=fix(n+1);
for  ip=1:n;
    b(ip)=a(ip,ip);
    d(ip)=b(ip);
    z(ip)=0.0d0;
end;  ip=fix(n+1);
nrot=0;
for  i=1:50;
    sm=0.0d0;
    for  ip=1:n-1;
        for  iq=ip+1:n;
            sm=sm+abs(a(ip,iq));
        end;  iq=fix(n+1);
    end;  ip=fix(n-1+1);
    if(sm == 0.)
        return;
    end;
    if(i < 4)
        tresh=(0.2d0.*sm./n.^2);
    else;
        tresh=0.0d0;
    end;
    for  ip=1:n-1;
        for  iq=ip+1:n;
            g=(100.0d0.*abs(a(ip,iq)));
            if((i > 4)&&(abs(d(ip))+g == abs(d(ip)))&&(abs(d(iq))+g == abs(d(iq))))
                a(ip,iq)=0.0d0;
            elseif(abs(a(ip,iq)) > tresh);
                h=d(iq)-d(ip);
                if(abs(h)+g == abs(h))
                    t=a(ip,iq)./h;
                else;
                    theta=0.5d0.*h./a(ip,iq);
                    t=1.0d0./(abs(theta)+sqrt(1.0d0+theta.^2));
                    if(theta < 0.0d0)
                        t=-t;
                    end;
                end;
                c=1.0d0./sqrt(1+t.^2);
                s=t.*c;
                tau=s./(1.0d0+c);
                h=t.*a(ip,iq);
                z(ip)=z(ip)-h;
                z(iq)=z(iq)+h;
                d(ip)=d(ip)-h;
                d(iq)=d(iq)+h;
                a(ip,iq)=0.0d0;
                for  j=1:ip-1;
                    g=a(j,ip);
                    h=a(j,iq);
                    a(j,ip)=g-s.*(h+g.*tau);
                    a(j,iq)=h+s.*(g-h.*tau);
                end;  j=fix(ip-1+1);
                for  j=ip+1:iq-1;
                    g=a(ip,j);
                    h=a(j,iq);
                    a(ip,j)=g-s.*(h+g.*tau);
                    a(j,iq)=h+s.*(g-h.*tau);
                end;  j=fix(iq-1+1);
                for  j=iq+1:n;
                    g=a(ip,j);
                    h=a(iq,j);
                    a(ip,j)=g-s.*(h+g.*tau);
                    a(iq,j)=h+s.*(g-h.*tau);
                end;  j=fix(n+1);
                for  j=1:n;
                    g=v(j,ip);
                    h=v(j,iq);
                    v(j,ip)=g-s.*(h+g.*tau);
                    v(j,iq)=h+s.*(g-h.*tau);
                end;  j=fix(n+1);
                nrot=fix(nrot+1);
            end;
        end;  iq=fix(n+1);
    end;  ip=fix(n-1+1);
    for  ip=1:n;
        b(ip)=b(ip)+z(ip);
        d(ip)=b(ip);
        z(ip)=0.0d0;
    end;  ip=fix(n+1);
end;  i=fix(50+1);

pause 'too many iterations in jacobi';

return;
end %subroutine jacobi
%  (C) Copr. 1986-92 Numerical Recipes Software !00,]v45.3.