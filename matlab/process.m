function dat = process(dat)
%PROCESS is processing an array of datapoints for independance. 
%   This process uses the quantile as a threshold and dist 
%   as the distance required between peaks in order to gain independance. 
%   It also checks to be sure the user's starting threshold is not too far below the quantile threshold.
%
%   DAT = PROCESS(DAT) accepts the the DAT structure defined in the
%   READINFO and returns the same object with the threshold info.
%
%OBRE project, 2011.

% Ver. 0.1 by Moiseev Igor, 27/01/2011
% Converted to Matlab from process.f, written by Andrew P. Wilson, Daltech, June 1999

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

% initialization 
%n=0; i=0; j=0; 
%qtint=0; lastp=0; qthres=0; maxy=0;
%x       = zeros(1,dat.NN);
%peak    = zeros(1,dat.NN); 
counter=0; 

% getting data to analyse
y = dat.serie(:,dat.dim);

% sorting the serie
sorted = sort(y,1,'ascend');

qtint = fix(dat.NN*dat.quant/100.0);
qthres = sorted(qtint);

j=1;
lastp = -1000;
while( j<=dat.NN )
    if( y(j) >= qthres && (j-lastp) > dat.dist )
        peak(1) = y(j);
        maxy  = y(j);
        lastp = j;
        maxy_j = j;
        
        i = 1; j = j+1;
        while( (j <= dat.NN && y(j) >= qthres) || ((j-lastp) < dat.dist && j <= dat.NN) )
            if(y(j) >= qthres)
                i = i+1;
                peak(i) = y(j);
                lastp = j;
                if(peak(i) >= maxy)
                    maxy = peak(i);
                    maxy_j = j;
                end
            end
            j = j+1;
        end % END WHILE
        
        counter = counter + 1;
        x(counter) = maxy;
        peak_j(counter) = maxy_j;
        j = j+1;
    else
        j = j+1;
    end
end

% Saving new independant data sets back into original variables to be returned to the main program
dat.peak = x;
dat.peak_j = peak_j;
% Number of peaks extracted from datas
dat.NN = counter;

% dialog on the correct thresold selection
if( (qthres-dat.thre) > 1.0 )
    fprintf(1,'%s %0.15g \n','All of the processed points are above:   ', qthres);
    fprintf(1,'%s %0.15g \n','You have chosen a starting threshold of: ', dat.thre);
    fprintf(1,'%0.15g \n','');
    answer = input('Do you wish to change your starting threshold? (yes/no): ', 's');
    if(strcmp(deblank(answer),deblank('yes')))
        dat.thre = input('What would you like the new threshold to be? (number): ');
    end;
end;

return;

end % END FUNCTION process
