function [alfa, beta, nex, finalthreshold, vref_obre, vref_mom, vref_pwm, weights_out, errors, dat] = obre(initialthreshold, returnperiod, DATINIT)
%OBRE Compute the MMs, PWMs, and OBREs of the parameters of a    
%   GP distribution for a given set of data.  These parameters
%   are then used to determine the weights and p-values of the
%   ten highest data points.  All data points are then subjected
%   to two goodness of fit equations to examine how well they fit
%   the GP distribution.
%
%   Generalized Pareto Distribution:
%
%       GP = THRE + ALPHA/BETA*(1 - P^BETA);
%   
%   THRE    value of the threshold used in the analysis
%   ALPHA   dispersion or scale parameter
%   BETA    mode of the extreme value distribution or shape parameter
%   P       percentile of the quantile of interest
%
%   VREF speed, the wind speed estimation with a GP distribution parameters
%   
%       VREF = THRE + ALPHA/BETA*(1 - (N/M*T)^(-BETA));
%
%   THRE    value of the threshold used in the analysis
%   ALPHA   dispersion or scale parameter
%   BETA    mode of the extreme value distribution or shape parameter
%   NEX     number of exceedances over the selected threshold
%   M       M is the number of years of record
%   T       Return period
%
%   [ALPHA, BETA, N, FINALTHRESHOLD, VREF_OBRE, VREF_MOM, VREF_PWM, WEIGHTS_OUT] 
%       = OBRE(INITIALTHRESOLD, RETURNPERIOD, DATINIT)
%
%   Obtains the OBREs estimates using Ronchetti's algorithm with 
%   PWM estimates as it's starting values.
%
%   RETURN VALUES
%
%       APLHA           dispersion estimation. Vector of DAT.rep elements.
%       BETA            shape estimation.  Vector of DAT.rep elements.
%       NEX             number of exceedances. Vector of DAT.rep elements.
%       FINALTHRESHOLD  the estimated threshold. Vector of DAT.rep elements.
%       VREF_OBRE       VREF speed estimation by OBRE. Vector of DAT.rep elements.
%       VREF_MOM        VREF speed estimation by MoM. Vector of DAT.rep elements.
%       VREF_PWM        VREF speed estimation by PWM. Vector of DAT.rep elements.
%       WEIGHTS_OUT     Weights vector
%       ERRORS          The errors vector
%       DAT             The last in the iteration the DAT structure
%
%   INPUT PARAMETERS
%       
%       The array with data will be requested upon execution. It should be
%       the MATLAB complatible format MAT file. 
%
%       INITIALTHRESOLD is the initial threshold for the OBRE algorithm
%       RETURNPERIOD is the return period in the VREF formula
%       DATINIT structure (optional)
%           dim     Dimention to be elaborated
%           NN      Desired Number of Data Points 
%           incr    Increment for Threshold
%           rep     Number of increases of Threshold
%           c       Desired bound on the influence function
%           epsiln  The perscribed epsilon needed to determine closeness
%           factor  The integration factor used in the step function to 
%                   speed up the program.  (1.0 is slow; 50.0 is fast)
%           sim     Desired Number of simulations performed
%           answer  Has the data been preprocessed for independance? (yes/no)
%           quant   What quantile do you wish to make the threshold?
%           dist    Minimum distance you wish there to be between peaks?
%           p       On what upper proportion do you want the Anderson
%                   Darling statistic to be performed(0.0 to 1.0)?
%           PROC    Boolean, corresponds to DAT.ANSWER
%
%   EXAMPLE:
%       [alfa, beta, nex, finalthreshold, vref] = obre(25.8, 50, struct('incr', 0.2, 'rep', 5));
%
%
%OBRE project, 2011.

% Ver. 0.1 by Moiseev Igor, 27/01/2011
% Converted to Matlab from realdata.f 
%  Single Precision Version 3.0
%  Written by Joanna E. Mills, TUNS, October 1996
%  Modified by D.J. Dupuis, February 1997
%  Modified by Andrew P. Wilson, May 1999

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

clc;
fprintf(1,'%s \n','This is OBRE for Matlab, ver. 0.1');
fprintf(1,'%s \n','Usage: obre(initialthreshold, returnperiod, DATINIT)');
fprintf(1,'%s \n','');

if nargin<2
  error('MATLAB:obre:NotEnoughInputs','Not enough input arguments.'); 
end

if nargin==3 && ~isstruct(DATINIT)
    error('MATLAB:obre', 'The DATINIT should be the structure to be processed.');
else
    DATINIT.thre = initialthreshold;
end

%% Read data file
[f1,pathload] = uigetfile('*.mat','Load data');
file = strcat(pathload, f1);
clear f1;

fprintf(1,'File to be ellaborated: %s \n', file);
fprintf(1,'%s \n','');

% Init DAT structure from the datas FILE
global dat;
dat = readinfo(file, DATINIT);

fprintf(1,'%s \n','Basic DAT structure info:');
fprintf(1,'   %s %0.15g \n','Number of observations: ', dat.NN);
fprintf(1,'   %s %0.15g \n','Number of simulations:  ', dat.sim);
fprintf(1,'   %s %0.15g \n','Starting Value of thre: ', dat.thre);
fprintf(1,'   %s %0.15g \n','Increment amount:       ', dat.incr);
fprintf(1,'   %s %0.15g \n','Number of repeats:      ', dat.rep);
fprintf(1,'   %s %0.15g \n','Bound on the IF (c):    ', dat.c);
fprintf(1,'   %s %0.15g \n','Prescribed epsilon:     ', dat.epsiln);
fprintf(1,'   %s %0.15g \n','Integration factor:     ', dat.factor);


%% Compute Mean Excess

y = sort(dat.peak);

thre2 = dat.thre;
numx = fix(y(dat.NN)-thre2);
dat.mean_excess = [];

for k=1:numx
    i=1;
    
    % if point is below the threshold, check the next point
    while( y(i) <= thre2 )
        i = i+1;
    end
    
    % when points above threshold are found, find the excess
    for j=i:dat.NN
        xmrp(j-i+1) = y(j) - thre2;
    end
    
    n = dat.NN - i + 1;
    
    dat.mean_excess = [dat.mean_excess; thre2,vmean(xmrp,n)];   %! *** find mean excess ***

    thre2 = thre2 + 1.0;

end % END FOR  k=1:numx

%% Begin repetitions

alfa = zeros(1,dat.rep);
beta = alfa;
nex = alfa;
finalthreshold = alfa;
vref_obre = alfa;
errors = struct;

for k=1:dat.rep
    i=1;
    while y(i) <= dat.thre
        i = i+1;
    end
    
    for j=i:dat.NN
        x(j-i+1) = y(j) - dat.thre;
    end
    
    n = dat.NN-i+1;
    
    % showing these exceedances so that we can do cross validation
    % and or construction with other data sets.
    fprintf(1,'\n%s \n','*****************************************');
    fprintf(1,'%s %0.15g \n','   Value of thre: ',dat.thre);
    fprintf(1,'%s %0.15g \n','   n: ',n);
    
    % Parameter Estimation via Method 1:  MoM
    [avals2(k),bvals2(k)] = meth1(n,x);
    fprintf(1,'\n%s \n',     '   MoM  Estimates:');
    fprintf(1,'%s %0.15g \n','   alpha: ',avals2(k));
    fprintf(1,'%s %0.15g \n','   beta:  ',bvals2(k));
    %fprintf(1,'%s %0.15g \n','   halving s.d. alpha: ',MoM_al);
    %fprintf(1,'%s %0.15g \n','   halving s.d. beta: ',MoM_be);
    
    % Step one:  finding initial parameter values
    % Parameter Estimation via Method 2:  PWM
    
    [avals3(k),bvals3(k)] = meth2(n,x);
    fprintf(1,'\n%s \n',     '   PWM  Estimates:');
    fprintf(1,'%s %0.15g \n','   alpha: ',avals3(k));
    fprintf(1,'%s %0.15g \n','   beta:  ',bvals3(k));
    %fprintf(1,'%s %0.15g \n','   halving s.d. alpha: ',PWM_al);
    %fprintf(1,'%s %0.15g \n','   halving s.d. beta: ',PWM_be);
    
    % ************************* Begin OBRE algorithm **************************
    % *-----------------------------------------------------------------------*
    % *   PWM estimates are the starting values for this OBRE algorithm.      *
    % *-----------------------------------------------------------------------*
    % *************************************************************************

    % Step Two:  finding the A matrix 
    % Calculation of the fisher information matrix 
    % and its the eigenvalues and eigenvectors

    %call fishinfo(avals3(k),bvals3(k),x,n,fish)
    dat.fish = fishinfo(avals3(k),bvals3(k),x,n);
    
    %call jacobi(fish,np,np,ei,ev,nrot)
    
    [ei,ev] = jacobi(dat.fish,2,2,[],[],0);

    % Take the reciprocal of the squareroot of each eigenvalue ***
                        
    ei(1)=1.0/sqrt(ei(1));
    ei(2)=1.0/sqrt(ei(2));
    
    % Transpose the fish matrix to arrive at amat                   ***
    % Calculation of initial values for avec and amat               ***
    
    avec(1) = 0;
    avec(2) = 0;
    amat(1,1)=(ei(1)*ev(1,1)*ev(1,1))+(ei(2)*ev(1,2)*ev(1,2));
    amat(1,2)=(ei(1)*ev(2,1)*ev(1,1))+(ei(2)*ev(2,2)*ev(1,2));
    amat(2,1)=(ei(1)*ev(1,1)*ev(2,1))+(ei(2)*ev(1,2)*ev(2,2));
    amat(2,2)=(ei(1)*ev(2,1)*ev(2,1))+(ei(2)*ev(2,2)*ev(2,2));

    % Steps three, four, and five of algorithm
    
    [avec, amat, fin, weights] = step345(n,x,avals3(k),bvals3(k),avec,amat,false);

    % ***************************End OBRE algorithm*****************************
    
    % Store final values in estimate array
    if( fin(1)<= 0 || fin(2)>= 1 )
        fprintf(1,'\n%s \n', '   PWM  Estimates:');
        fprintf(1,'\n%s \n', '       OBRE algorithm failed for this threshold.');
        fprintf(1,'\n%s \n', '       Poor estimate obtained.');
        fprintf(1,'\n%s \n', '       Continuing to next threshold.');
        %goto 299
    end
    
    avals4(k) = fin(1);
    bvals4(k) = fin(2);

    fprintf(1,'\n%s ',         '   OBRE Estimates:');
    fprintf(1,'\n%s %0.15g',   '   alpha:', avals4(k));
    fprintf(1,'\n%s %0.15g\n', '   beta: ', bvals4(k));

    % ************************Computing Standard Errors*************************

    % compute the Asymptotic Covariance Matrix
    acm = covmat(n,x,avals4(k),bvals4(k),avec,weights);
    
    fprintf(1,'\n%s ',         '   OBRE Stderr Estimates:');
    fprintf(1,'\n%s %0.15g',   '   s.d. alpha: ', sqrt(acm(1,1)/n));
    fprintf(1,'\n%s %0.15g',   '   s.d. beta:  ', sqrt(acm(2,2)/n));

    fprintf(1,'\n');
    
    fprintf(1,'\n%s %0.15g',   '   95th quantile: ', quangp(dat.thre, avals4(k), bvals4(k), 0.05));
    fprintf(1,'\n%s %0.15g',   '   98th quantile: ', quangp(dat.thre, avals4(k), bvals4(k), 0.02));
    fprintf(1,'\n%s %0.15g',   '      s.d. q_.05: ', qustderr(avals4(k), bvals4(k), 0.05, acm, n));
    fprintf(1,'\n%s %0.15g',   '      s.d. q_.02: ', qustderr(avals4(k), bvals4(k), 0.02, acm, n));
    
    errors.stdalfa(k) = sqrt(acm(1,1)/n);
    errors.stdbeta(k) = sqrt(acm(2,2)/n);
    errors.quant95(k) = quangp(dat.thre, avals4(k), bvals4(k), 0.05);
    errors.quant98(k) = quangp(dat.thre, avals4(k), bvals4(k), 0.02);
    errors.stdquant95(k) = qustderr(avals4(k), bvals4(k), 0.05, acm, n);
    errors.stdquant98(k) = qustderr(avals4(k), bvals4(k), 0.02, acm, n);
    
    
    if(bvals4(k) >= 0)
        fprintf(1,'\n%s %0.15g',   '     Upper bound: ', dat.thre+(avals4(k)/bvals4(k)) );
        fprintf(1,'\n%s %0.15g',   'Highest observed point: ', x(n)+dat.thre );
        
        errors.uppper_bound(k) = dat.thre+(avals4(k)/bvals4(k));
        errors.highest_pnt(k)  = x(n)+dat.thre;
                
        if( x(n) >= (avals4(k)/bvals4(k)) )
            fprintf(1,'\n');
            fprintf(1,'\nX~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~X');
            fprintf(1,'\n| WARNING:  Non-feasable estimates! |');
            fprintf(1,'\nX___________________________________X');
            %go to 299 
        end
    end
    fprintf(1,'\n');
    fprintf(1,'\n---------------WEIGHTS---------------');

    if(length(weights) > 1)
        weights_out(k,1:length(weights)) = weights;
    end
    
    for i=1:n
        if(weights(i) ~= 1.0)
            fprintf(1,'\n%s %0.15g %s %0.15g',   'Observation ', i, '  -weight=', weights(i) );
        end
    end
        
    fprintf(1,'\n*****************************************\n');
    
    % *** OBRE and goodness of fit algorithms for this threshold ***
    
    % A review of methods to calculate extreme wind speeds, page 6
    % If the exceedance process is assumed to be Poisson with rate
    % lambda per year, an unbiased estimate of lambda is n/M where n is
    % the total number of exceedances over the selected
    % threshold x, and M is the number of years of record.

    M = dat.NN/52560;
    
    vref_mom(k)  = dat.thre + avals2(k)/bvals2(k)*(1 - (n/M*returnperiod)^(-bvals2(k)));
    vref_pwm(k)  = dat.thre + avals3(k)/bvals3(k)*(1 - (n/M*returnperiod)^(-bvals3(k)));
    vref_obre(k) = dat.thre + avals4(k)/bvals4(k)*(1 - (n/M*returnperiod)^(-bvals4(k)));

    fprintf(1,'\n%s %0.15g',   'Estimated VREF by OBRE: ', vref_obre(k));
    fprintf(1,'\n%s %0.15g',   'Estimated VREF by MoM:  ', vref_mom(k));
    fprintf(1,'\n%s %0.15g\n',   'Estimated VREF by PWM:  ', vref_pwm(k));
    
    % TODO have not been converted
    %call goodoffit(thre,x,n,iout4,p,avals4(k),bvals4(k),we)
    %call gpobre(n,sim,beta,we,iout3)

    alfa(k) = avals4(k);
    beta(k) = bvals4(k);
    nex(k) = n;
    finalthreshold(k) = dat.thre;
    
    % Increase threshold for the next repetition
    dat.thre = dat.thre + dat.incr;
    
end % END FOR k=1:dat.rep

return;
end
