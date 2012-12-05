function [avec, amat, fin, weights] = step345(n, y, alpha_init, beta_init, avec, amat, simul)
%STEP345 Performs Steps three, four and five of Ronchetti's algorithm
%
%    n            number of observations
%    y            vector containing the simulated observations
%    alpha_init   initial alpha value of parameter vector
%    beta_init    initial beta value of parameter vector
%    avec         value of 'a' vector
%    amat         value of A matrix
%    simul        boolean to simulate
%
%    fin          vector of parameter (alpha, beta) estimates obtained (OUTPUT)
%    weights      vector of weights corresponding to OBRE estimates (OUTPUT)
%
%OBRE project, 2011.

% Ver. 0.1 by Moiseev Igor, 27/01/2011
% Converted to Matlab from step.f
%  Written by Joanna E. Mills, TUNS, October 1996
%  Modified by D.J. Dupuis, 24 July 2000

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

% *********************Initializations*********************************

% The DAT structure refs
global dat;

% Parameters
p=2; l1max=20; l2max=15; 

negmax=10;
quant=0.0001;

% maxobs= 5000;
% err=0;
% m1=zeros(p,p); m2=zeros(p,p);
% del=zeros(1,p);
% s1=zeros(1,maxobs);
% scalemlv=zeros(p,p);
%
% s2=zeros(1,maxobs);
rdel1=0; rdel2=0; lm=0;
% top1=0; top2=0; bottom=0;
% limit=0; limpr=0;
%
rsmall = 1.0d-10;
done = false; 
loop1 = 0; loop2 = 0;

alpha = alpha_init;
beta = beta_init;

% *******************Iteration Starting Point***************************
% LOOP B
while( (loop2 > l2max) || rdel1 < dat.epsiln && rdel2 < dat.epsiln)
    
    % LOOP A
    while(~done && loop1 < l1max)
        loop1 = loop1+1;
        
        % Store 'old' values of a and A
        oldav = avec;
        oldam = amat;
        
        % Step three:  Calculating the A matrix and the a vector
        % The integral for m2 is calculated via numerical integration
        % If k is presently less than -1.0, use only k = -1.0
        % when computing the quantile
        
        if(beta < -1)
            lm = quangp(0, alpha, -1, quant*10);
        else
            lm = quangp(0, alpha, beta, quant);
        end
        
        % The maximum number of loops equals the limit divided by dat.factor plus one for good measure
        maxlp = fix(lm/dat.factor) + 1;
        m2 = zeros(p,p);
        limit = 0;
        
        if(beta <= 0)
            for i=1:maxlp
                limpr = limit;
                limit = limit + dat.factor;
                
                m2(1,1) = m2(1,1) + dgsint(@(x)m211(x, alpha, beta, avec, amat),limpr,limit);
                m2(1,2) = m2(1,2) + dgsint(@(x)m212(x, alpha, beta, avec, amat),limpr,limit);
                m2(2,1) = m2(2,1) + dgsint(@(x)m221(x, alpha, beta, avec, amat),limpr,limit);
                m2(2,2) = m2(2,2) + dgsint(@(x)m222(x, alpha, beta, avec, amat),limpr,limit);
            end
        else
            for i=1:maxlp
                limpr=limit;
                limit=limit + dat.factor;
                
                % when k is positive, function is bounded at alpha/beta
                if(limit >(alpha/beta))
                    limit = alpha/beta;
                end
                
                m2(1,1) = m2(1,1) + dgsint(@(x)m211(x, alpha, beta, avec, amat),limpr,limit);
                m2(1,2) = m2(1,2) + dgsint(@(x)m212(x, alpha, beta, avec, amat),limpr,limit);
                m2(2,1) = m2(2,1) + dgsint(@(x)m221(x, alpha, beta, avec, amat),limpr,limit);
                m2(2,2) = m2(2,2) + dgsint(@(x)m222(x, alpha, beta, avec, amat),limpr,limit);
            end
        end
        
        % Calculate the inverse of m2
        m2inv = inv(m2);
        
        % Calculate new amat (A^t*A = m2^-1)
        amat(1,1) = sqrt(m2inv(1,1));
        amat(2,1) = 0;
        amat(1,2) = m2inv(2,1)/amat(1,1);
        amat(2,2) = sqrt(m2inv(2,2)-(amat(1,2)*amat(1,2)));
        
        % Calculate new avec (Using the equation on top of pg.260)
        % Computing integrals of the 'a' vector equation
        bottom = 0;
        top1 = 0;
        top2 = 0;
        limit = 0;
        
        if(beta <= 0)
            for i=1:maxlp;
                limpr = limit;
                limit = limit + dat.factor;
                bottom = bottom + dgsint(@(x)denom(x, alpha, beta, avec, amat),  limpr,limit);
                top1   = top1   + dgsint(@(x)numer1(x, alpha, beta, avec, amat), limpr,limit);
                top2   = top2   + dgsint(@(x)numer2(x, alpha, beta, avec, amat), limpr,limit);
            end
        else
            for i=1:maxlp;
                limpr=limit;
                limit=limit + dat.factor;
                
                % when k is positive, function is bounded at alpha/beta
                
                if(limit >(alpha/beta))
                    limit = alpha./beta;
                end;
                bottom = bottom + dgsint(@(x)denom(x, alpha, beta, avec, amat),  limpr,limit);
                top1   = top1   + dgsint(@(x)numer1(x, alpha, beta, avec, amat), limpr,limit);
                top2   = top2   + dgsint(@(x)numer2(x, alpha, beta, avec, amat), limpr,limit);
            end
        end
        
        if(abs(bottom) > rsmall)
            avec(1) = top1./bottom;
            avec(2) = top2./bottom;
        else
            avec(1) = 0;
            avec(2) = 0;
        end
        
        % 'Check' to see if a and A close to previous values
        % (i.e. they have converged).  If not we loop back to the start
        done = compare(avec,oldav,amat,oldam);
        
    end % END OF LOOP A
    
    % At this point we have converged to values for a and A
    % Reset loop1 so it can be performed again in loop2
    loop2 = loop2+1;
    
    % Calculate new weights corresponding to new a, and A
    % Calculate m1 via numerical integration and deltapar
    % via averaging
    % Step four:  finding deltheta
    deltapar = zeros(1,p);
    s1 = zeros(1,n); s2 = s1; weights = s1;
    
    for i=1:n
        s1(i) = score1(alpha, beta, y(i), false);
        s2(i) = score2(alpha, beta, y(i), false);
        weights(i)  = weight(y(i), alpha, beta, avec, amat);
        deltapar(1) = deltapar(1)+(s1(i)-avec(1)) * weights(i);
        deltapar(2) = deltapar(2)+(s2(i)-avec(2)) * weights(i);
    end
    
    % write(iout,*),'delpar',deltapar(1),deltapar(2)
    
    deltapar(1) = deltapar(1)/n;
    deltapar(2) = deltapar(2)/n;
    
    if(beta < -1.0d0)
        lm = quangp(0.0,(alpha),-1.0,(quant*10.0));
    else
        lm = quangp(0.0,(alpha),(beta),quant);
    end
    
    % The maximum number of loops equals the limit divided by dat.factor
    % plus one for good measure
    maxlp = fix(lm./dat.factor) + 1;
    m1    = zeros(p,p);
    limit = 0;
    
    if(beta <= 0)
        for i=1:maxlp
            limpr = limit;
            limit = limit + dat.factor;
            
            m1(1,1) = m1(1,1) + dgsint(@(x)m111(x, alpha, beta, avec, amat),limpr,limit);
            m1(1,2) = m1(1,2) + dgsint(@(x)m112(x, alpha, beta, avec, amat),limpr,limit);
            m1(2,1) = m1(2,1) + dgsint(@(x)m121(x, alpha, beta, avec, amat),limpr,limit);
            m1(2,2) = m1(2,2) + dgsint(@(x)m122(x, alpha, beta, avec, amat),limpr,limit);
        end
    else
        for i=1:maxlp
            limpr = limit;
            limit = limit + dat.factor;
            
            % when k is positive, function is bounded at alpha/beta
            if(limit >(alpha/beta))
                limit = alpha/beta;
            end;
            
            m1(1,1) = m1(1,1) + dgsint(@(x)m111(x, alpha, beta, avec, amat),limpr,limit);
            m1(1,2) = m1(1,2) + dgsint(@(x)m112(x, alpha, beta, avec, amat),limpr,limit);
            m1(2,1) = m1(2,1) + dgsint(@(x)m121(x, alpha, beta, avec, amat),limpr,limit);
            m1(2,2) = m1(2,2) + dgsint(@(x)m122(x, alpha, beta, avec, amat),limpr,limit);
        end
    end
    
    % Calculate the inverse of M1
    m1inv = inv(m1);
    
    del(1) = m1inv(1,1)*deltapar(1)+m1inv(1,2)*deltapar(2);
    del(2) = m1inv(2,1)*deltapar(1)+m1inv(2,2)*deltapar(2);
    
    % These lines throw us out of the Loop2 if we have been in it for too
    % long or if there is convergence (based on the epsilon parameter)
    
    %if( loop2 > l2max )
    %    go to 110;
    %end;
    
    rdel1 = abs(del(1)/alpha);
    rdel2 = abs(del(2)/beta);
    
    %if(rdel1 < dat.epsiln && rdel2 < dat.epsiln)
    %    go to 110;
    %end
    
    % If not, check to see if negative alpha estimate has been obtained
    % First check to see if we have reasonable estimates
    % If not, return to main program (redo simulation)
    negcount=0;
    
    % label 90
    while( (beta+del(2)) > 1.0 && negcount < negmax )
        del(1) = del(1)/2.0;
        
        % Half estimates to find convergence
        del(2) = del(2)/2.0;
        negcount = fix(negcount + 1);
    end % END WHILE label 90
    
    if(negcount >= negmax)
        if(simul)
            fprintf(1,'%s \n','   Redoing Simulation');
            fprintf(1,'%s \n','   -->Impossible K<--');
            return;
        else
            fprintf(1,'%s \n','The OBRE algorithm could not converge.');
            fprintf(1,'%s \n','       Using latest estimates');
            break;
        end
    end
    
    negcount=0;
    % label 91
    while( (alpha+del(1)) < 0 && negcount < negmax)
        % Half estimates to find convergence
        del(1) = del(1)/2.0;
        del(2) = del(2)/2.0;
        negcount = fix(negcount + 1);
    end
    
    if(negcount >= negmax)
        if(simul)
            fprintf(1,'%s \n','   Redoing Simulation');
            fprintf(1,'%s \n','  -->Negative Alpha<--');
            return;
        else
            fprintf(1,'%s \n','The OBRE algorithm could not converge.');
            fprintf(1,'%s \n','       Using latest estimates');
            break;
        end
    end
    
    negcount=0;
    % label 92
    while( beta_init >= 0.2 && (beta+del(2)) < -0.5 && negcount < negmax )
        % Half estimates to find convergence
        del(1) = del(1)/2.0;
        del(2) = del(2)/2.0;
        negcount = fix(negcount + 1);
    end
    
    if(negcount >= negmax)
        if(simul)
            fprintf(1,'%s \n','   Redoing Simulation');
            fprintf(1,'%s \n','  -->Bad k estimate<--');
            return;
        else
            fprintf(1,'%s \n','The OBRE algorithm could not converge.');
            fprintf(1,'%s \n','       Using latest estimates');
            break;
        end
    end
    
    negcount = 0;
    % label 93
    while( (500.0*alpha_init)<(alpha+del(1)) && negcount < negmax)
        % Half estimates to find convergence
        del(1) = del(1)/2.0;
        del(2) = del(2)/2.0;
        negcount=fix(negcount + 1);
    end

    if(negcount >= negmax)
        if(simul)
            fprintf(1,'%s \n','   Redoing Simulation');
            fprintf(1,'%s \n','   -->Bad alpha estimate<--');
            return;
        else
            fprintf(1,'%s \n','The OBRE algorithm could not converge.');
            fprintf(1,'%s \n','       Using latest estimates');
            break;
        end
    end
    
    % Step five: theta = theta + deltheta
    alpha = alpha + del(1);
    beta  = beta + del(2);

end % END WHILE, END OF LOOP B

% LABEL 110
% Return the values of alpha, beta, the a vector and the A matrix
fin(1)  = alpha;
fin(2)  = beta;

return;
end % END FUNCTION step


%  ************************************************************************
%  *                                                                      *
%  *    function m111, m112, m121, m122, m211, m212, m221, m222           *
%  *                                                                      *
%  ************************************************************************
%  Double Precision Version 1.0
%  Written by Andrew P. Wilson, June 1999
%
%  PURPOSE  These functions implement equations to be integrated via another
%           program.  These equations are ones seen as the second equation
%           found on page 260.  Each of the following programs implement
%           one cell each of the 2 by 2, m1 and m2 matrices.  The numbers
%           after the m indicate which matrix and cell the function
%           corresponds to. (m121 is matrix m1, row 2, column 1)
%
%      x    real value at which each function is to be evaluated
%
%--------------------------------------------------------------------------

function [m211result,feas] = m211(x, alpha, beta, avec, amat)
%M211 moment
feas = false;

[sco1,feas] = score1(alpha, beta, x,feas);

m211result = (sco1-avec(1))*(sco1-avec(1))*weight(x, alpha, beta, avec, amat)*density(x, alpha, beta)*weight(x, alpha, beta, avec, amat);

return;
end % END FUNCTION m211


function [m212result,feas] = m212(x, alpha, beta, avec, amat)
%M212 moment
feas=false;

[sco1,feas1] = score1(alpha, beta, x,feas);
[sco2,feas2] = score2(alpha, beta, x,feas);

m212result =(sco2-avec(2))*(sco1-avec(1))*weight(x, alpha, beta, avec, amat)*density(x, alpha, beta)*weight(x, alpha, beta, avec, amat);

feas = feas1 * feas2;

return;
end % END FUNCTION m212


function [m221result,feas] = m221(x, alpha, beta, avec, amat)
%M221 moment
feas=false;

[sco1,feas1] = score1(alpha, beta, x,feas);
[sco2,feas2] = score2(alpha, beta, x,feas);

m221result =(sco1-avec(1))*(sco2-avec(2))*weight(x, alpha, beta, avec, amat)*density(x, alpha, beta)*weight(x, alpha, beta, avec, amat);

feas = feas1 * feas2;

return;
end % END FUNCTION m221


function [m222result,feas] = m222(x, alpha, beta, avec, amat)
%M222 moment
feas = false;

[sco2,feas] = score2(alpha, beta, x,feas);

m222result =(sco2-avec(2))*(sco2-avec(2))*weight(x, alpha, beta, avec, amat)*density(x, alpha, beta)*weight(x, alpha, beta, avec, amat);

return;
end % END FUNCTION m222


function [m111result,feas] = m111(x, alpha, beta, avec, amat)
%M111 moment
feas=false;

[sco1,feas] = score1(alpha, beta, x,feas);

m111result =(sco1-avec(1))*(sco1-avec(1))*weight(x, alpha, beta, avec, amat)*density(x, alpha, beta);

return;
end % END FUNCTION m111


function [m112result,feas] = m112(x, alpha, beta, avec, amat)
%M112 moment
feas=false;

[sco1,feas1] = score1(alpha, beta, x,feas);
[sco2,feas2] = score2(alpha, beta, x,feas);

m112result =(sco2-avec(2))*(sco1-avec(1))*weight(x, alpha, beta, avec, amat)*density(x, alpha, beta);

feas = feas1 * feas2;

return;
end % END FUNCTION m112


function [m121result,feas] = m121(x, alpha, beta, avec, amat)
%M121 moment
feas=false;

[sco1,feas1] = score1(alpha, beta, x,feas);
[sco2,feas2] = score2(alpha, beta, x,feas);

m121result =(sco1-avec(1))*(sco2-avec(2))*weight(x, alpha, beta, avec, amat)*density(x, alpha, beta);

feas = feas1 * feas2;

return;
end % END FUNCTION m121


function [m122result,feas] = m122(x, alpha, beta, avec, amat)
%M122 moment
feas=false;

[sco2,feas] = score2(alpha, beta, x, feas);

m122result =(sco2-avec(2))*(sco2-avec(2))*weight(x, alpha, beta, avec, amat)*density(x, alpha, beta);

return;
end % END FUNCTION m122


%  *******************************************************************
%  *                                                                 *
%  *            function numer1, numer2, denom                       *
%  *                                                                 *
%  *******************************************************************
%  Double Precision Version 1.0
%  Written by Joanna E. Mills, September 1996
%
%  PURPOSE  These function perform the equations of the a vector as
%           seen on the top of pg.260.  These functions are used to
%           find the limits of integration in the step subroutine.
%
%     x  real value at which each function is to be evaluated
%--------------------------------------------------------------------------

function [numer1result] = numer1(x, alpha, beta, avec, amat)
    numer1result = score1(alpha, beta, x, false)*weight(x, alpha, beta, avec, amat)*density(x, alpha, beta);
return;
end % END FUNCTION numer1

function [numer2result] = numer2(x, alpha, beta, avec, amat)
    numer2result = score2(alpha, beta, x, false)*weight(x, alpha, beta, avec, amat)*density(x, alpha, beta);
return;
end % END FUNCTION numer2

function [denomresult,x] = denom(x, alpha, beta, avec, amat)
    denomresult = weight(x, alpha, beta, avec, amat)*density(x, alpha, beta);
return;
end % END FUNCTION denom
