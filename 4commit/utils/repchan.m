% Replace bad channel(s) using spherical spline interpolation
% X  = repchan(X, chans,chanlocs)

function X  = repchan(X, chans,chanlocs)
% Optional inputs:
%   'nterms'  - scalar int > 0 number of terms {default 50}
%   'm'       - scalar int > 1 m {default 4}
%   'nframes' - scalar int > 0 number of frames to interpolate
%               vectorized (trade-off between speed and memory usage)
%               {default 1000}
%
%   Channel coordinates should be located on the surface of a (unit) sphere.


if nargin < 1
    help repchan;
    return
end
if isempty(X)
    error('Cannot process empty dataset.');
end



nterms = 50;
mm = 4;
nframes = size(X,2);
ntr=size(X,3);
nbchan=size(X,1);


% Find channels with location coordinates
chanArray = find(~sum(cellfun('isempty', {chanlocs.X; chanlocs.Y; chanlocs.Z})));
% Unit sphere
E = unitsph([chanlocs(chanArray).X; chanlocs(chanArray).Y; chanlocs(chanArray).Z]');
% Channel location coordinate matrices E and F
F = E(ismember(chanArray, chans), :);
if size(F, 1) ~= length(chans)
    error('No channel location coordinates for channel(s) to interpolate.');
end
E = E(~ismember(chanArray, chans), :);
chanArray = setdiff(chanArray, chans);
% G matrix
[Ginv, g] = sserpgfcn(E, F, 'sp', 0, nterms, mm);
% Reshape X
if ntr > 1
    X = reshape(X, size(X,1), []);
end
blockArray = [1 : nframes : size(X, 2) size(X, 2) + 1];
for iBlock = 1 : length(blockArray) - 1
    % C matrix
    C = sserpweights(X(chanArray, blockArray(iBlock) : blockArray(iBlock + 1) - 1), Ginv);
    % Interpolation
    X(chans, blockArray(iBlock) : blockArray(iBlock + 1) - 1) = sserp(C, g, 'sp');
end

% Reshape X
if ntr > 1
    X = reshape(X, nbchan, nframes, ntr);
end













% arg2str() - Convert function argument cell array or structure to
%             EEGLAB compatible history string
%
% Usage:
%   >> str = arg2str(arg);
%
% Inputs:
%   arg       - cell array arguments with parameter value pairs or
%               structure arguments
%
% Outputs:
%   str       - history string
%
% Author: Andreas Widmann, University of Leipzig, 2006

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2006 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Id: arg2str.m 8 2006-04-03 09:05:23Z widmann $

function str = arg2str(arg)

if isstruct(arg)
    str = '';
    fieldArray = fieldnames(arg);
    arg = struct2cell(arg);
    for iArg = 1:length(arg)
        if ischar(arg{iArg})
            str = [str '''' fieldArray{iArg} ''', ''' arg{iArg} ''''];
        elseif isnumeric(arg{iArg}) | islogical(arg{iArg})
            str = [str '''' fieldArray{iArg} ''', ' mat2str(arg{iArg})];
        elseif iscell(arg{iArg})
            str = [str '''' fieldArray{iArg} ''', ' arg2str(arg{iArg})];
        end
        if iArg < length(arg)
            str = [str ', '];
        end
    end
elseif iscell(arg)
    str = '{';
    for iArg = 1:length(arg)
        if ischar(arg{iArg})
            str = [str '''' arg{iArg} ''''];
        elseif isnumeric(arg{iArg}) | islogical(arg{iArg})
            str = [str mat2str(arg{iArg})];
        elseif iscell(arg{iArg})
            str = [str arg2str(arg{iArg})];
        end
        if iArg < length(arg)
            str = [str ' '];
        else
            str = [str '}'];
        end
    end
end




% sserp() - Perform spherical spline interpolation
%
% Usage:
%   >> erpData = sserp(C, g, type)
%
% Inputs:
%   C         - nbchan + 1 by pnts matrix weights
%   g         - locs by nbchan matrix g(cos(E, F))
%
% Optional inputs:
%   type      - string type of interpolation. 'sp' (scalp potential,
%               �V), 'scd' (scalp current density, scaled to mA/m^3), or
%               'lap' (surface Laplacian, scaled to mV/m^2) {default:
%               'sp'}
%
% Outputs:
%   erpData   - locs by pnts matrix interpolated data
%
% References:
%   [1] Perrin, F., Pernier, J., Bertrand, O., & Echallier, J. F.
%       (1989). Spherical splines for scalp potential and current
%       density mapping. Electroencephalography and Clinical
%       Neurophysiology, 72, 184-187
%   [2] Perrin, F., Pernier, J., Bertrand, O., & Echallier, J. F.
%       (1990). Corrigenda EEG 02274. Electroencephalography and
%       Clinical Neurophysiology, 76, 565
%   [3] Kayser, J., & Tenke, C. E. (2006). Principal components analysis
%       of Laplacian waveforms as a generic method for identifying ERP
%       generator patterns: I. Evaluation with auditory oddball tasks.
%       Clinical Neurophysiology, 117, 348-368
%   [4] Weber, D. L. (2001). Scalp current density and source current
%       modelling. Retrieved March 26, 2006, from
%       dnl.ucsf.edu/users/dweber/dweber_docs/eeg_scd.html
%   [5] Ferree, T. C. (2000). Spline Interpolation of the Scalp EEG.
%       Retrieved March 26, 2006, from
%       www.egi.com/Technotes/SplineInterpolation.pdf
%   [6] Ferree, T. C., & Srinivasan, R. (2000). Theory and Calculation
%       of the Scalp Surface Laplacian. Retrieved March 26, 2006, from
%       http://www.egi.com/Technotes/SurfaceLaplacian.pdf
%
% Author: Andreas Widmann, University of Leipzig, 2006
%
% See also:
%   sserpgfcn(), sserpweights()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2006 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Id: sserp.m 10 2006-04-12 15:57:01Z widmann $

function erpData = sserp(C, g, type)

% Arguments
if nargin < 3 || isempty(type)
    type = 'sp';
end
if nargin < 2 || isempty(g) || isempty(C)
    error('Not enough input arguments.')
end

% Constants
CONDUCTIVITY = -0.45; % Siemens/meter
HEAD_RADIUS = 0.1; % meter
SCALE = 1e-3; % �A/m^3 to mA/m^3 and �V/m^2 to mV/m^2

% Interpolation
switch type
    case 'sp'
        % Perrin et al., 1989, eqn. (1)
        erpData = C(ones(1, size(g, 1)), :) + g * C(2:end, :);
    case 'scd'
        % Perrin et al., 1989, eqn. (5), negative 2D spherical laplacian
        erpData = -g * C(2:end, :) / HEAD_RADIUS ^ 2 * CONDUCTIVITY * SCALE;
    case 'lap'
        % Perrin et al., 1989, eqn. (5), negative 2D spherical laplacian
        erpData = -g * C(2:end, :) / HEAD_RADIUS ^ 2 * SCALE;
    otherwise
        error('Unrecognized or ambiguous interpolation type specified.');
end


% sserpgfcn() - Compute G matrices used for spherical spline
%               interpolation
%
% Usage:
%   >> [Ginv, g, G] = sserpgfcn(E, F, type, lambda, nterms, m)
%
% Inputs:
%   E         - nbchan by 3 matrix with cartesian channel location
%               coordinates x, y, z of measured channels
%   F         - locs by 3 matrix with cartesian channel location
%               coordinates x, y, z to interpolate in columns
%
% Optional inputs:
%   type      - string type of interpolation. 'sp' (scalp potential),
%               'scd' (scalp current density), or 'lap' (surface
%               Laplacian) {default 'sp'}
%   lambda    - scalar smoothing factor (commonly used values are 1e-7
%               for sp and 1e-5 for scd) {default 0}
%   nterms    - scalar int > 0 number of terms {default 50}
%   m         - scalar int > 1 m {default 4}
%
% Outputs:
%   Ginv      - nbchan + 1 by nbchan + 1 matrix padded inverse of
%               g(cos(E, E))
%   g         - locs by nbchan matrix g(cos(E, F))
%   G         - nbchan by nbchan matrix g(cos(E, E))
%
% Note:
%   Recursive version of gfcn is tens of times faster than using MATLAB
%   legendre function. Sserpgfcn implicitly understands that channel
%   coordinates in E and F are located on the surface of a unit sphere.
%   It will return wrong results if they are not. Use unitsph() to
%   project channel coordinates to a unit sphere surface.
%
% References:
%   [1] Perrin, F., Pernier, J., Bertrand, O., & Echallier, J. F.
%       (1989). Spherical splines for scalp potential and current
%       density mapping. Electroencephalography and Clinical
%       Neurophysiology, 72, 184-187
%   [2] Perrin, F., Pernier, J., Bertrand, O., & Echallier, J. F.
%       (1990). Corrigenda EEG 02274. Electroencephalography and
%       Clinical Neurophysiology, 76, 565
%   [3] Kayser, J., & Tenke, C. E. (2006). Principal components analysis
%       of Laplacian waveforms as a generic method for identifying ERP
%       generator patterns: I. Evaluation with auditory oddball tasks.
%       Clinical Neurophysiology, 117, 348-368
%   [4] Weber, D. L. (2001). Scalp current density and source current
%       modelling. Retrieved March 26, 2006, from
%       dnl.ucsf.edu/users/dweber/dweber_docs/eeg_scd.html
%   [5] Ferree, T. C. (2000). Spline Interpolation of the Scalp EEG.
%       Retrieved March 26, 2006, from
%       www.egi.com/Technotes/SplineInterpolation.pdf
%   [6] Ferree, T. C., & Srinivasan, R. (2000). Theory and Calculation
%       of the Scalp Surface Laplacian. Retrieved March 26, 2006, from
%       http://www.egi.com/Technotes/SurfaceLaplacian.pdf
%
% Author: Andreas Widmann, University of Leipzig, 2006
%
% See also:
%   sserpweights(), sserp(), unitsph()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2006 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Id: sserpgfcn.m 10 2006-04-12 15:57:01Z widmann $

function [Ginv, g, G] = sserpgfcn(E, F, type, lambda, nterms, m)

if nargin < 6 || isempty(m)
    m = 4;
end
if nargin < 5 || isempty(nterms)
    nterms = 50;
end
if nargin < 4 || isempty(lambda)
    lambda = 0;
end
if nargin < 3 || isempty(type)
    type = 'sp';
end
if nargin < 2
    error('Not enough input arguments.')
end

% Cosines, quaternion based analog to Perrin et al., 1989, eqn. (4)
x = E * E';

% G matrix
G = gfcn(x, nterms, m);

% Pad, add smoothing constant to diagonale, and invert G
Ginv = inv([0 ones(1, size(G, 2)); ones(size(G, 1), 1) G + eye(size(G)) * lambda]);

% Cosines, quaternion based analog to Perrin et al., 1989, eqn. (4)
x = F * E';

% g matrix
switch type
    case 'sp'
        g = gfcn(x, nterms, m); % Perrin et al., 1989, eqn. (3)
    case {'scd' 'lap'}
        g = gfcn(x, nterms, m - 1); % Perrin et al., 1990, eqn. after eqn. (5)
    otherwise
        error('Unrecognized or ambiguous interpolation type specified.');
end

function [G] = gfcn(x, nterms, m)
    P = cat(3, ones(size(x)), x);
    G = 3 / 2 ^ m * P(:, :, 2);
    for n = 2:nterms
        % nth degree legendre polynomial; Perrin et al., 1989, eqn. after eqn. (4)
        P(:, :, 3) = ((2 * n - 1) * x .* P(:, :, 2) - (n - 1) * P(:, :, 1)) / n;
        P(:, :, 1) = [];
        % Perrin et al., 1989, eqn. (3)
        G = G + (2 * n + 1) / (n ^ m * (n + 1) ^ m) * P(:, :, 2);
    end
    G = G / (4 * pi);

    
    
    
    
    
    % sserpweights() - Compute C matrix used for spherical spline
%                  interpolation
%
% Usage:
%   >> C = sserpweights(data, Ginv)
%
% Inputs:
%   data  - nbchan by pnts matrix of measured data
%   Ginv  - matrix padded inverse of g(cos(E, E))
%
% Outputs:
%   C     - nbchan + 1 by pnts matrix weights
%
% References:
%   [1] Perrin, F., Pernier, J., Bertrand, O., & Echallier, J. F.
%       (1989). Spherical splines for scalp potential and current
%       density mapping. Electroencephalography and Clinical
%       Neurophysiology, 72, 184-187
%   [2] Perrin, F., Pernier, J., Bertrand, O., & Echallier, J. F.
%       (1990). Corrigenda EEG 02274. Electroencephalography and
%       Clinical Neurophysiology, 76, 565
%   [3] Kayser, J., & Tenke, C. E. (2006). Principal components analysis
%       of Laplacian waveforms as a generic method for identifying ERP
%       generator patterns: I. Evaluation with auditory oddball tasks.
%       Clinical Neurophysiology, 117, 348-368
%   [4] Weber, D. L. (2001). Scalp current density and source current
%       modelling. Retrieved March 26, 2006, from
%       dnl.ucsf.edu/users/dweber/dweber_docs/eeg_scd.html
%   [5] Ferree, T. C. (2000). Spline Interpolation of the Scalp EEG.
%       Retrieved March 26, 2006, from
%       www.egi.com/Technotes/SplineInterpolation.pdf
%   [6] Ferree, T. C., & Srinivasan, R. (2000). Theory and Calculation
%       of the Scalp Surface Laplacian. Retrieved March 26, 2006, from
%       http://www.egi.com/Technotes/SurfaceLaplacian.pdf
%
% Author: Andreas Widmann, University of Leipzig, 2006
%
% See also:
%   sserpgfcn(), sserp()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2006 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Id: sserpweights.m 10 2006-04-12 15:57:01Z widmann $

function C = sserpweights(data, Ginv)

% Arguments
if nargin < 2
    error('Not enough input arguments.')
end

% C matrix
C = Ginv * [zeros(1, size(data, 2)); double(data)]; % Perrin et al., 1989, eqn. (2)








% unitsph() - Re-center channel location coordinates and project to
%             unit sphere surface
%
% Usage:
%   >> E = unitsph(E, thresh)
%
% Inputs:
%   E     - nbchan by 3 matrix with cartesian channel location
%           coordinates x, y, z
%
% Optional inputs:
%   thres - scalar threshold < abs(radius - 1) {default 1e-6}
%
% Outputs:
%   E     - nbchan by 3 matrix with cartesian channel location
%           coordinates x, y, z re-centered to best fitting sphere and
%           projected to unit sphere surface
%
% Author: Andreas Widmann, University of Leipzig, 2006
%
% See also:
%   chancenter()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2006 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Id: unitsph.m 10 2006-04-12 15:57:01Z widmann $

function E = unitsph(E, thresh)

if nargin < 2
    thresh = 1e-6;
end

[th, phi, r] = cart2sph(E(:, 1), E(:, 2), E(:, 3));

if any(abs(r - 1) > thresh)
%     warning('Channel coordinates not located on unit sphere surface. Recentering and projecting.')
    [E(:, 1), E(:, 2), E(:, 3)] = chancenter(E(:, 1), E(:, 2), E(:, 3), []); % Re-center
    [th, phi] = cart2sph(E(:, 1), E(:, 2), E(:, 3));
    [E(:, 1), E(:, 2), E(:, 3)] = sph2cart(th, phi, 1); % Project to unit sphere
end



% chancenter() - recenter cartesian X,Y,Z channel coordinates
%
% Usage:  >> [x y z newcenter] = chancenter(x,y,z,center); 
%
% Optional inputs:
%    x,y,z     = 3D coordintates of the channels
%    center    = [X Y Z] known center different from [0 0 0]
%                [] will optimize the center location according
%                to the best sphere. Default is [0 0 0].
%
% Note: 4th input gui is obsolete. Use pop_chancenter instead.
%
% Authors: Arnaud Delorme, Luca Finelli & Scott Makeig SCCN/INC/UCSD,
%          La Jolla, 11/1999-03/2002 
%
% See also: spherror(), cart2topo()

% Copyright (C) 11/1999 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% 3-16-00 improved help message -sm
% 1-25-02 put spherror subfunction inside chancenter -ad
% 1-25-02 help pops-up if no arguments -ad
% 01-25-02 reformated help & license -ad 
% 02-13-02 now center fitting works, need spherror outside chancenter -lf
% 02-14-02 radii are squeezed of squeeze in to fit topoplot circle -lf
% 03-31-02 center fitting is optional
% 04-01-02 automatic squeeze calculation -ad & sm
 
function [ x, y, z, newcenter, optim] = chancenter( x, y, z, center, gui)

optim = 0;

if nargin<4
    help chancenter
    return;
end;

if nargin > 4 && gui
    error('Chancenter: 4th input'' gui'' is obsolete. Use pop_chancenter instead');
else 
	if isempty(center)
		optim = 1;
		center = [0 0 0];
	end;
end;

options = {'MaxFunEvals',1000*length(center)};
x = x - center(1);  % center the data
y = y - center(2);
z = z - center(3);
radius = (sqrt(x.^2+y.^2+z.^2));   % assume xyz values are on a sphere
wobble = std(radius);              % test if xyz values are on a sphere
% fprintf('Radius values: %g (mean) +/- %g (std)\n',mean(radius),wobble);
newcenter = center;

if  wobble/mean(radius) > 0.01 && optim==1
	% Find center
	% ----------------------------------------------
% 	fprintf('Optimizing center position...\n');
	kk=0;
	while wobble/mean(radius) > 0.01 && kk<5
		newcenter = fminsearch('spherror',center,options,x,y,z);
		nx = x - newcenter(1);  % re-center the data
		ny = y - newcenter(2);
		nz = z - newcenter(3);
		nradius = (sqrt(nx.^2+ny.^2+nz.^2));   % assume xyz values are on a sphere
		newobble = std(nradius);   
		if newobble<wobble
			center=newcenter;
% 			fprintf('Wobble too strong (%3.2g%%)! Re-centering data on (%g,%g,%g)\n',...
% 					100*wobble/mean(radius),newcenter(1),newcenter(2),newcenter(3))
			x = nx;  % re-center the data
			y = ny;
			z = nz;
			radius=nradius;
			wobble=newobble;
			kk=kk+1;
		else
			newcenter = center;
			kk=5;
		end
	end
% 	fprintf('Wobble (%3.2g%%) after centering data on (%g,%g,%g)\n',...
% 			100*wobble/mean(radius),center(1),center(2), ...
% 			center(3))
	%else
	%  fprintf('Wobble (%3.2g%%) after centering data on (%g,%g,%g)\n',...
	%              100*wobble/mean(radius),center(1),center(2),center(3))
end






function [x,fval,exitflag,output] = fminsearch(funfcn,x,options,varargin)
%FMINSEARCH Multidimensional unconstrained nonlinear minimization (Nelder-Mead).
%   X = FMINSEARCH(FUN,X0) starts at X0 and attempts to find a local minimizer 
%   X of the function FUN.  FUN is a function handle.  FUN accepts input X and 
%   returns a scalar function value F evaluated at X. X0 can be a scalar, vector 
%   or matrix.
%
%   X = FMINSEARCH(FUN,X0,OPTIONS)  minimizes with the default optimization
%   parameters replaced by values in the structure OPTIONS, created
%   with the OPTIMSET function.  See OPTIMSET for details.  FMINSEARCH uses
%   these options: Display, TolX, TolFun, MaxFunEvals, MaxIter, FunValCheck,
%   PlotFcns, and OutputFcn.
%
%   X = FMINSEARCH(PROBLEM) finds the minimum for PROBLEM. PROBLEM is a
%   structure with the function FUN in PROBLEM.objective, the start point
%   in PROBLEM.x0, the options structure in PROBLEM.options, and solver
%   name 'fminsearch' in PROBLEM.solver. The PROBLEM structure must have
%   all the fields.
%
%   [X,FVAL]= FMINSEARCH(...) returns the value of the objective function,
%   described in FUN, at X.
%
%   [X,FVAL,EXITFLAG] = FMINSEARCH(...) returns an EXITFLAG that describes 
%   the exit condition of FMINSEARCH. Possible values of EXITFLAG and the 
%   corresponding exit conditions are
%
%    1  Maximum coordinate difference between current best point and other
%       points in simplex is less than or equal to TolX, and corresponding 
%       difference in function values is less than or equal to TolFun.
%    0  Maximum number of function evaluations or iterations reached.
%   -1  Algorithm terminated by the output function.
%
%   [X,FVAL,EXITFLAG,OUTPUT] = FMINSEARCH(...) returns a structure
%   OUTPUT with the number of iterations taken in OUTPUT.iterations, the
%   number of function evaluations in OUTPUT.funcCount, the algorithm name 
%   in OUTPUT.algorithm, and the exit message in OUTPUT.message.
%
%   Examples
%     FUN can be specified using @:
%        X = fminsearch(@sin,3)
%     finds a minimum of the SIN function near 3.
%     In this case, SIN is a function that returns a scalar function value
%     SIN evaluated at X.
%
%     FUN can also be an anonymous function:
%        X = fminsearch(@(x) norm(x),[1;2;3])
%     returns a point near the minimizer [0;0;0].
%
%   If FUN is parameterized, you can use anonymous functions to capture the 
%   problem-dependent parameters. Suppose you want to optimize the objective     
%   given in the function myfun, which is parameterized by its second argument c. 
%   Here myfun is an M-file function such as
%
%     function f = myfun(x,c)
%     f = x(1)^2 + c*x(2)^2;
%
%   To optimize for a specific value of c, first assign the value to c. Then 
%   create a one-argument anonymous function that captures that value of c 
%   and calls myfun with two arguments. Finally, pass this anonymous function 
%   to FMINSEARCH:
%    
%     c = 1.5; % define parameter first
%     x = fminsearch(@(x) myfun(x,c),[0.3;1])
%
%   FMINSEARCH uses the Nelder-Mead simplex (direct search) method.
%
%   See also OPTIMSET, FMINBND, FUNCTION_HANDLE.

%   Reference: Jeffrey C. Lagarias, James A. Reeds, Margaret H. Wright,
%   Paul E. Wright, "Convergence Properties of the Nelder-Mead Simplex
%   Method in Low Dimensions", SIAM Journal of Optimization, 9(1):
%   p.112-147, 1998.

%   Copyright 1984-2007 The MathWorks, Inc.
%   $Revision: 1.21.4.16 $  $Date: 2008/10/31 06:19:57 $


defaultopt = struct('Display','notify','MaxIter','200*numberOfVariables',...
    'MaxFunEvals','200*numberOfVariables','TolX',1e-4,'TolFun',1e-4, ...
    'FunValCheck','off','OutputFcn',[],'PlotFcns',[]);

% If just 'defaults' passed in, return the default options in X
if nargin==1 && nargout <= 1 && isequal(funfcn,'defaults')
    x = defaultopt;
    return
end

if nargin<3, options = []; end

% Detect problem structure input
if nargin == 1
    if isa(funfcn,'struct') 
        [funfcn,x,options] = separateOptimStruct(funfcn);
    else % Single input and non-structure
        error('MATLAB:fminsearch:InputArg','The input to FMINSEARCH should be either a structure with valid fields or consist of at least two arguments.');
    end
end

if nargin == 0
    error('MATLAB:fminsearch:NotEnoughInputs',...
        'FMINSEARCH requires at least two input arguments');
end


% Check for non-double inputs
if ~isa(x,'double')
  error('MATLAB:fminsearch:NonDoubleInput', ...
         'FMINSEARCH only accepts inputs of data type double.')
end

n = numel(x);
numberOfVariables = n;

printtype = optimget(options,'Display',defaultopt,'fast');
tolx = optimget(options,'TolX',defaultopt,'fast');
tolf = optimget(options,'TolFun',defaultopt,'fast');
maxfun = optimget(options,'MaxFunEvals',defaultopt,'fast');
maxiter = optimget(options,'MaxIter',defaultopt,'fast');
funValCheck = strcmp(optimget(options,'FunValCheck',defaultopt,'fast'),'on');

% In case the defaults were gathered from calling: optimset('fminsearch'):
if ischar(maxfun)
    if isequal(lower(maxfun),'200*numberofvariables')
        maxfun = 200*numberOfVariables;
    else
        error('MATLAB:fminsearch:OptMaxFunEvalsNotInteger',...
            'Option ''MaxFunEvals'' must be an integer value if not the default.')
    end
end
if ischar(maxiter)
    if isequal(lower(maxiter),'200*numberofvariables')
        maxiter = 200*numberOfVariables;
    else
        error('MATLAB:fminsearch:OptMaxIterNotInteger',...
            'Option ''MaxIter'' must be an integer value if not the default.')
    end
end

switch printtype
    case {'notify','notify-detailed'}
        prnt = 1;
    case {'none','off'}
        prnt = 0;
    case {'iter','iter-detailed'}
        prnt = 3;
    case {'final','final-detailed'}
        prnt = 2;
    case 'simplex'
        prnt = 4;
    otherwise
        prnt = 1;
end
% Handle the output
outputfcn = optimget(options,'OutputFcn',defaultopt,'fast');
if isempty(outputfcn)
    haveoutputfcn = false;
else
    haveoutputfcn = true;
    xOutputfcn = x; % Last x passed to outputfcn; has the input x's shape
    % Parse OutputFcn which is needed to support cell array syntax for OutputFcn.
    outputfcn = createCellArrayOfFunctions(outputfcn,'OutputFcn');
end

% Handle the plot
plotfcns = optimget(options,'PlotFcns',defaultopt,'fast');
if isempty(plotfcns)
    haveplotfcn = false;
else
    haveplotfcn = true;
    xOutputfcn = x; % Last x passed to plotfcns; has the input x's shape
    % Parse PlotFcns which is needed to support cell array syntax for PlotFcns.
    plotfcns = createCellArrayOfFunctions(plotfcns,'PlotFcns');
end

header = ' Iteration   Func-count     min f(x)         Procedure';

% Convert to function handle as needed.
funfcn = fcnchk(funfcn,length(varargin));
% Add a wrapper function to check for Inf/NaN/complex values
if funValCheck
    % Add a wrapper function, CHECKFUN, to check for NaN/complex values without
    % having to change the calls that look like this:
    % f = funfcn(x,varargin{:});
    % x is the first argument to CHECKFUN, then the user's function,
    % then the elements of varargin. To accomplish this we need to add the 
    % user's function to the beginning of varargin, and change funfcn to be
    % CHECKFUN.
    varargin = {funfcn, varargin{:}};
    funfcn = @checkfun;
end

n = numel(x);

% Initialize parameters
rho = 1; chi = 2; psi = 0.5; sigma = 0.5;
onesn = ones(1,n);
two2np1 = 2:n+1;
one2n = 1:n;

% Set up a simplex near the initial guess.
xin = x(:); % Force xin to be a column vector
v = zeros(n,n+1); fv = zeros(1,n+1);
v(:,1) = xin;    % Place input guess in the simplex! (credit L.Pfeffer at Stanford)
x(:) = xin;    % Change x to the form expected by funfcn
fv(:,1) = funfcn(x,varargin{:});
func_evals = 1;
itercount = 0;
how = '';
% Initial simplex setup continues later

% Initialize the output and plot functions.
if haveoutputfcn || haveplotfcn
    [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,v(:,1),xOutputfcn,'init',itercount, ...
        func_evals, how, fv(:,1),varargin{:});
    if stop
        [x,fval,exitflag,output] = cleanUpInterrupt(xOutputfcn,optimValues);
        if  prnt > 0
            disp(output.message)
        end
        return;
    end
end

% Print out initial f(x) as 0th iteration
if prnt == 3
    disp(' ')
    disp(header)
    disp(sprintf(' %5.0f        %5.0f     %12.6g         %s', itercount, func_evals, fv(1), how));
elseif prnt == 4
    clc
    formatsave = get(0,{'format','formatspacing'});
    format compact
    format short e
    disp(' ')
    disp(how)
    v
    fv
    func_evals
end
% OutputFcn and PlotFcns call
if haveoutputfcn || haveplotfcn
    [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,v(:,1),xOutputfcn,'iter',itercount, ...
        func_evals, how, fv(:,1),varargin{:});
    if stop  % Stop per user request.
        [x,fval,exitflag,output] = cleanUpInterrupt(xOutputfcn,optimValues);
        if  prnt > 0
            disp(output.message)
        end
        return;
    end
end

% Continue setting up the initial simplex.
% Following improvement suggested by L.Pfeffer at Stanford
usual_delta = 0.05;             % 5 percent deltas for non-zero terms
zero_term_delta = 0.00025;      % Even smaller delta for zero elements of x
for j = 1:n
    y = xin;
    if y(j) ~= 0
        y(j) = (1 + usual_delta)*y(j);
    else
        y(j) = zero_term_delta;
    end
    v(:,j+1) = y;
    x(:) = y; f = funfcn(x,varargin{:});
    fv(1,j+1) = f;
end

% sort so v(1,:) has the lowest function value
[fv,j] = sort(fv);
v = v(:,j);

how = 'initial simplex';
itercount = itercount + 1;
func_evals = n+1;
if prnt == 3
    disp(sprintf(' %5.0f        %5.0f     %12.6g         %s', itercount, func_evals, fv(1), how))
elseif prnt == 4
    disp(' ')
    disp(how)
    v
    fv
    func_evals
end
% OutputFcn and PlotFcns call
if haveoutputfcn || haveplotfcn
    [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,v(:,1),xOutputfcn,'iter',itercount, ...
        func_evals, how, fv(:,1),varargin{:});
    if stop  % Stop per user request.
        [x,fval,exitflag,output] = cleanUpInterrupt(xOutputfcn,optimValues);
        if  prnt > 0
            disp(output.message)
        end
        return;
    end
end
exitflag = 1;

% Main algorithm: iterate until 
% (a) the maximum coordinate difference between the current best point and the 
% other points in the simplex is less than or equal to TolX. Specifically,
% until max(||v2-v1||,||v2-v1||,...,||v(n+1)-v1||) <= TolX,
% where ||.|| is the infinity-norm, and v1 holds the 
% vertex with the current lowest value; AND
% (b) the corresponding difference in function values is less than or equal
% to TolFun. (Cannot use OR instead of AND.)
% The iteration stops if the maximum number of iterations or function evaluations 
% are exceeded
while func_evals < maxfun && itercount < maxiter
    if max(abs(fv(1)-fv(two2np1))) <= max(tolf,10*eps(fv(1))) && ...
            max(max(abs(v(:,two2np1)-v(:,onesn)))) <= max(tolx,10*eps(max(v(:,1))))
        break
    end
    
    % Compute the reflection point
    
    % xbar = average of the n (NOT n+1) best points
    xbar = sum(v(:,one2n), 2)/n;
    xr = (1 + rho)*xbar - rho*v(:,end);
    x(:) = xr; fxr = funfcn(x,varargin{:});
    func_evals = func_evals+1;
    
    if fxr < fv(:,1)
        % Calculate the expansion point
        xe = (1 + rho*chi)*xbar - rho*chi*v(:,end);
        x(:) = xe; fxe = funfcn(x,varargin{:});
        func_evals = func_evals+1;
        if fxe < fxr
            v(:,end) = xe;
            fv(:,end) = fxe;
            how = 'expand';
        else
            v(:,end) = xr;
            fv(:,end) = fxr;
            how = 'reflect';
        end
    else % fv(:,1) <= fxr
        if fxr < fv(:,n)
            v(:,end) = xr;
            fv(:,end) = fxr;
            how = 'reflect';
        else % fxr >= fv(:,n)
            % Perform contraction
            if fxr < fv(:,end)
                % Perform an outside contraction
                xc = (1 + psi*rho)*xbar - psi*rho*v(:,end);
                x(:) = xc; fxc = funfcn(x,varargin{:});
                func_evals = func_evals+1;
                
                if fxc <= fxr
                    v(:,end) = xc;
                    fv(:,end) = fxc;
                    how = 'contract outside';
                else
                    % perform a shrink
                    how = 'shrink';
                end
            else
                % Perform an inside contraction
                xcc = (1-psi)*xbar + psi*v(:,end);
                x(:) = xcc; fxcc = funfcn(x,varargin{:});
                func_evals = func_evals+1;
                
                if fxcc < fv(:,end)
                    v(:,end) = xcc;
                    fv(:,end) = fxcc;
                    how = 'contract inside';
                else
                    % perform a shrink
                    how = 'shrink';
                end
            end
            if strcmp(how,'shrink')
                for j=two2np1
                    v(:,j)=v(:,1)+sigma*(v(:,j) - v(:,1));
                    x(:) = v(:,j); fv(:,j) = funfcn(x,varargin{:});
                end
                func_evals = func_evals + n;
            end
        end
    end
    [fv,j] = sort(fv);
    v = v(:,j);
    itercount = itercount + 1;
    if prnt == 3
        disp(sprintf(' %5.0f        %5.0f     %12.6g         %s', itercount, func_evals, fv(1), how))
    elseif prnt == 4
        disp(' ')
        disp(how)
        v
        fv
        func_evals
    end
    % OutputFcn and PlotFcns call
    if haveoutputfcn || haveplotfcn
        [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,v(:,1),xOutputfcn,'iter',itercount, ...
            func_evals, how, fv(:,1),varargin{:});
        if stop  % Stop per user request.
            [x,fval,exitflag,output] = cleanUpInterrupt(xOutputfcn,optimValues);
            if  prnt > 0
                disp(output.message)
            end
            return;
        end
    end
end   % while

x(:) = v(:,1);
fval = fv(:,1);

if prnt == 4,
    % reset format
    set(0,{'format','formatspacing'},formatsave);
end
output.iterations = itercount;
output.funcCount = func_evals;
output.algorithm = 'Nelder-Mead simplex direct search';

% OutputFcn and PlotFcns call
if haveoutputfcn || haveplotfcn
    callOutputAndPlotFcns(outputfcn,plotfcns,x,xOutputfcn,'done',itercount, func_evals, how, fval, varargin{:});
end

if func_evals >= maxfun
    msg = sprintf(['Exiting: Maximum number of function evaluations has been exceeded\n' ...
                   '         - increase MaxFunEvals option.\n' ...
                   '         Current function value: %f \n'], fval);
    if prnt > 0
        disp(' ')
        disp(msg)
    end
    exitflag = 0;
elseif itercount >= maxiter
    msg = sprintf(['Exiting: Maximum number of iterations has been exceeded\n' ... 
                   '         - increase MaxIter option.\n' ...
                   '         Current function value: %f \n'], fval);
    if prnt > 0
        disp(' ')
        disp(msg)
    end
    exitflag = 0;
else
    msg = ...
      sprintf(['Optimization terminated:\n', ...
               ' the current x satisfies the termination criteria using OPTIONS.TolX of %e \n' ...
               ' and F(X) satisfies the convergence criteria using OPTIONS.TolFun of %e \n'], ...
               tolx, tolf);
    if prnt > 1
        disp(' ')
        disp(msg)
    end
    exitflag = 1;
end

output.message = msg;

%--------------------------------------------------------------------------
function [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,x,xOutputfcn,state,iter,...
    numf,how,f,varargin)
% CALLOUTPUTANDPLOTFCNS assigns values to the struct OptimValues and then calls the
% outputfcn/plotfcns.
%
% state - can have the values 'init','iter', or 'done'.

% For the 'done' state we do not check the value of 'stop' because the
% optimization is already done.
optimValues.iteration = iter;
optimValues.funccount = numf;
optimValues.fval = f;
optimValues.procedure = how;

xOutputfcn(:) = x;  % Set x to have user expected size
stop = false;
% Call output functions
if ~isempty(outputfcn)
    switch state
        case {'iter','init'}
            stop = callAllOptimOutputFcns(outputfcn,xOutputfcn,optimValues,state,varargin{:}) || stop;
        case 'done'
            callAllOptimOutputFcns(outputfcn,xOutputfcn,optimValues,state,varargin{:});
        otherwise
            error('MATLAB:fminsearch:InvalidState', ...
                'Unknown state in CALLOUTPUTANDPLOTFCNS.')
    end
end
% Call plot functions
if ~isempty(plotfcns)
    switch state
        case {'iter','init'}
            stop = callAllOptimPlotFcns(plotfcns,xOutputfcn,optimValues,state,varargin{:}) || stop;
        case 'done'
            callAllOptimPlotFcns(plotfcns,xOutputfcn,optimValues,state,varargin{:});
        otherwise
            error('MATLAB:fminsearch:InvalidState', ...
                'Unknown state in CALLOUTPUTANDPLOTFCNS.')
    end
end

%--------------------------------------------------------------------------
function [x,FVAL,EXITFLAG,OUTPUT] = cleanUpInterrupt(xOutputfcn,optimValues)
% CLEANUPINTERRUPT updates or sets all the output arguments of FMINBND when the optimization
% is interrupted.

x = xOutputfcn;
FVAL = optimValues.fval;
EXITFLAG = -1;
OUTPUT.iterations = optimValues.iteration;
OUTPUT.funcCount = optimValues.funccount;
OUTPUT.algorithm = 'golden section search, parabolic interpolation';
OUTPUT.message = 'Optimization terminated prematurely by user.';

%--------------------------------------------------------------------------
function f = checkfun(x,userfcn,varargin)
% CHECKFUN checks for complex or NaN results from userfcn.

f = userfcn(x,varargin{:});
% Note: we do not check for Inf as FMINSEARCH handles it naturally.
if isnan(f)
    error('MATLAB:fminsearch:checkfun:NaNFval', ...
        'User function ''%s'' returned NaN when evaluated;\n FMINSEARCH cannot continue.', ...
        localChar(userfcn));  
elseif ~isreal(f)
    error('MATLAB:fminsearch:checkfun:ComplexFval', ...
        'User function ''%s'' returned a complex value when evaluated;\n FMINSEARCH cannot continue.', ...
        localChar(userfcn));  
end

%--------------------------------------------------------------------------
function strfcn = localChar(fcn)
% Convert the fcn to a string for printing

if ischar(fcn)
    strfcn = fcn;
elseif isa(fcn,'inline')
    strfcn = char(fcn);
elseif isa(fcn,'function_handle')
    strfcn = func2str(fcn);
else
    try
        strfcn = char(fcn);
    catch
        strfcn = '(name not printable)';
    end
end































