function [dat] = readinfo(file, DATINIT)
%READINFO creates the config data structure for further processing
%   READINFO(FILE, DATINIT) returns DAT structure
%   FILE variable is the file with a path to be processed
%
%   The DAT structure
%
%       serie   The array to be loaded
%       dim     Dimention to be elaborated
%       NN      Desired Number of Data Points
%       thre    Threshold
%       incr    Increment for Threshold
%       rep     Number of increases of Threshold
%       c       Desired bound on the influence function
%       epsiln  The perscribed epsilon needed to determine closeness
%       factor  The integration factor used in the step function to 
%               speed up the program.  (1.0 is slow; 50.0 is fast)
%       sim     Desired Number of simulations performed
%       answer  Has the data been preprocessed for independance? (yes/no)
%       quant   What quantile do you wish to make the threshold?
%       dist    Minimum distance you wish there to be between peaks?
%       p       On what upper proportion do you want the Anderson
%               Darling statistic to be performed(0.0 to 1.0)?
%       PROC    Boolean, corresponds to DAT.ANSWER
%
%OBRE project, 2011.

% Ver. 0.1 by Moiseev Igor, 27/01/2011
% Converted to Matlab from readinfo.f
% Written by Joanna E. Mills, TUNS, October 1996
% Modified by D.J. Dupuis, February 1997
% Modified by Andrew P. Wilson, May 1999

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

% check if file exist
if(~exist(file,'file'))
    error('MATLAB:OBRE:readinfo', 'File does not exist.');
end

dat = load(file);

if( ~isstruct(dat) )
    error('MATLAB:OBRE:readinfo', 'Could not read the file %s or occured the error on data loading', file);
end

sizeofSerie = size(dat.serie);

% Dimension of serie to perform the stat. tests on
dat.dim    = 6;
if(isfield(DATINIT,'dim')), dat.dim = DATINIT.dim; end

dat.NN     = sizeofSerie(1,1);
if(isfield(DATINIT,'NN')), dat.NN = DATINIT.NN; end
    
dat.thre   = 25.80;
if(isfield(DATINIT,'thre')), dat.thre = DATINIT.thre; end

dat.incr   = 0.2;
if(isfield(DATINIT,'incr')), dat.incr = DATINIT.incr; end

dat.rep    = 10;
if(isfield(DATINIT,'rep')), dat.rep = DATINIT.rep; end

dat.c      = 4.0;
if(isfield(DATINIT,'c')), dat.c = DATINIT.c; end

dat.epsiln = 0.001;
if(isfield(DATINIT,'epsiln')), dat.epsiln = DATINIT.epsiln; end

dat.factor = 5.0;
if(isfield(DATINIT,'factor')), dat.factor = DATINIT.factor; end

dat.sim    = 200;
if(isfield(DATINIT,'sim')), dat.sim = DATINIT.sim; end

dat.answer = 'no';
if(isfield(DATINIT,'answer')), dat.answer = DATINIT.answer; end

dat.quant  = 95;
if(isfield(DATINIT,'quant')), dat.quant = DATINIT.quant; end

dat.dist   = 576;
if(isfield(DATINIT,'dist')), dat.dist = DATINIT.dist; end

dat.p      = 1.0;
if(isfield(DATINIT,'p')), dat.p = DATINIT.p; end

dat.PROC = false;
if( strcmp(dat.answer,'yes') )
    dat.PROC = true;
end

% lets check for default values
if(dat.p == 0.0d0)
    dat.p = 1.0d0;
end

if(dat.c == 0.0)
    dat.c = 4.0d0;
end;

if(dat.epsiln == 0.0)
    dat.epsiln = 0.001d0;
end

if(dat.factor == 0.0)
    dat.factor = 1.0d0;
end

if(dat.quant == 0)
    dat.quant = 90;
end;

% now will do the prcoess for finding threshold
if( ~dat.PROC )
    dat = process(dat);
else
    dat.peak = sort(dat.serie(:,dat.dim), 1, 'ascend');
end;

return;
end % END FUCNTION readinfo