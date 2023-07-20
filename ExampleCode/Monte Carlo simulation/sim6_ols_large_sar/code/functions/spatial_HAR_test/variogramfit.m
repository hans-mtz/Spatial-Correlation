function [a,c,n,S] = variogramfit(h,gammaexp,a0,c0,numobs,varargin)

% fit a theoretical variogram to an experimental variogram
%
% Syntax:
%
%     [a,c,n] = variogramfit(h,gammaexp,a0,c0)
%     [a,c,n] = variogramfit(h,gammaexp,a0,c0,numobs)
%     [a,c,n] = variogramfit(h,gammaexp,a0,c0,numobs,'pn','pv',...)
%     [a,c,n,S] = variogramfit(...)
%
% Description:
%
%     variogramfit performs a least squares fit of various theoretical 
%     variograms to an experimental, isotropic variogram. The user can
%     choose between various bounded (e.g. spherical) and unbounded (e.g.
%     exponential) models. A nugget variance can be modelled as well, but
%     higher nested models are not supported.
%
%     The function works best with the function fminsearchbnd available on
%     the FEX. You should download it from the File Exchange (File ID:
%     #8277). If you don't have fminsearchbnd, variogramfit uses
%     fminsearch. The problem with fminsearch is, that it might return 
%     negative variances or ranges.
%
%     The variogram fitting algorithm is in particular sensitive to initial
%     values below the optimal solution. In case you have no idea of
%     initial values variogramfit calculates initial values for you
%     (c0 = max(gammaexp); a0 = max(h)*2/3;). If this is a reasonable
%     guess remains to be answered. Hence, visually inspecting your data
%     and estimating a theoretical variogram by hand should always be
%     your first choice.
%
%     Note that for unbounded models, the supplied parameter a0 (range) is
%     the distance where gamma equals 95% of the sill variance. The
%     returned parameter a0, however, is the parameter r in the model. The
%     range at 95% of the sill variance is then approximately 3*r.
%
% Input arguments:
%
%     h         lag distance of the experimental variogram
%     gammaexp  experimental variogram values (gamma)
%     a0        initial value (scalar) for range
%     c0        initial value (scalar) for sill variance
%     numobs    number of observations per lag distance (used for weight
%               function)
%
% Output arguments:
%
%     a         range
%     c         sill
%     n         nugget (empty if nugget variance is not applied)
%     S         structure array with additional information
%               .model - theoretical variogram 
%               .h  - distance
%               .gamma  - experimental variogram values
%               .gammahat - estimated variogram values
%               .residuals - residuals
%               .Rs - R-square of goodness of fit
%               .weights - weights
%               .exitflag - see fminsearch
%               .algorithm - see fminsearch
%               .funcCount - see fminsearch
%               .iterations - see fminsearch
%               .message - see fminsearch
%
% Property name/property values:
% 
%     'model'   a string that defines the function that can be fitted 
%               to the experimental variogram. 
% 
%               Supported bounded functions are:
%               'blinear' (bounded linear) 
%               'circular' (circular model)
%               'spherical' (spherical model, =default)
%               'pentaspherical' (pentaspherical model)
% 
%               Supported unbounded functions are:
%               'exponential' (exponential model)
%               'gaussian' (gaussian variogram)
%               'whittle' Whittle's elementary correlation (involves a
%                         modified Bessel function of the second kind.
%               'stable' (stable models sensu Wackernagel 1995). Same as
%                         gaussian, but with different exponents. Supply 
%                         the exponent alpha (<2) in an additional pn,pv 
%                         pair: 
%                        'stablealpha',alpha (default = 1.5).
%               'matern' Matern model. Requires an additional pn,pv pair. 
%                        'nu',nu (shape parameter > 0, default = 1)
%                        Note that for particular values of nu the matern 
%                        model reduces to other authorized variogram models.
%                        nu = 0.5 : exponential model
%                        nu = 1 : Whittles model
%                        nu -> inf : Gaussian model
%               
%               See Webster and Oliver (2001) for an overview on variogram 
%               models. See Minasny & McBratney (2005) for an
%               introduction to the Matern variogram.
%           
%     'nugget'  initial value (scalar) for nugget variance. The default
%               value is []. In this case variogramfit doesn't fit a nugget
%               variance. 
% 
%     'plotit'  true (default), false: plot experimental and theoretical 
%               variogram together.
% 
%     'solver'  'fminsearchbnd' (default) same as fminsearch , but with  
%               bound constraints by transformation (function by John 
%               D'Errico, File ID: #8277 on the FEX). The advantage in 
%               applying fminsearchbnd is that upper and lower bound 
%               constraints can be applied. That prevents that nugget 
%               variance or range may become negative.           
%               'fminsearch'
%
%     'weightfun' 'none' (default). 'cressie85' and 'mcbratney86' require
%               you to include the number of observations per experimental
%               gamma value (as returned by VARIOGRAM). 
%               'cressie85' uses m(hi)/gammahat(hi)^2 as weights
%               'mcbratney86' uses m(hi)*gammaexp(hi)/gammahat(hi)^3
%               
%
% Example: fit a variogram to experimental data
%
%     load variogramexample
%     a0 = 15; % initial value: range 
%     c0 = 0.1; % initial value: sill 
%     [a,c,n] = variogramfit(h,gammaexp,a0,c0,[],...
%                            'solver','fminsearchbnd',...
%                            'nugget',0,...
%                            'plotit',true);
%
%           
% See also: VARIOGRAM, FMINSEARCH, FMINSEARCHBND
%           
%
% References: Wackernagel, H. (1995): Multivariate Geostatistics, Springer.
%             Webster, R., Oliver, M. (2001): Geostatistics for
%             Environmental Scientists. Wiley & Sons.
%             Minsasny, B., McBratney, A. B. (2005): The Mat�rn function as
%             general model for soil variograms. Geoderma, 3-4, 192-207.
% 
% Date: 8. February, 2010
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)




% check input arguments

if nargin == 0
    help variogramfit
    return
elseif nargin>0 && nargin < 2;
    error('Variogramfit:inputargs',...
          'wrong number of input arguments');
end
if ~exist('a0','var') || isempty(a0)
    a0 = max(h)*2/3;
end
if ~exist('c0','var') || isempty(c0)
    c0 = max(gammaexp);
end
if ~exist('numobs','var') || isempty(a0)
    numobs = [];
end
      

% check input parameters
params.model       = 'spherical';
params.nugget      = [];
params.plotit      = true;
params.solver      = 'fminsearchbnd';
params.stablealpha = 1.5;
params.weightfun   = 'cressie85';
params.nu          = 1;
params = parseargs(params,varargin{:});

% check if fminsearchbnd is in the search path
switch lower(params.solver)
    case 'fminsearchbnd'
        if ~exist('fminsearchbnd.m','file')==2
            params.solver = 'fminsearch';
            warning('Variogramfit:fminsearchbnd',...
            'fminsearchbnd was not found. fminsearch is used instead')
        end
end

% check if h and gammaexp are vectors and have the same size
if ~isvector(h) || ~isvector(gammaexp)
    error('Variogramfit:inputargs',...
          'h and gammaexp must be vectors');
end

% force column vectors
h = h(:);
gammaexp = gammaexp(:);

% check size of supplied vectors 
if numel(h) ~= numel(gammaexp)
    error('Variogramfit:inputargs',...
          'h and gammaexp must have same size');
end

% remove nans;
nans = isnan(h) | isnan(gammaexp);
if any(nans);
    h(nans) = [];
    gammaexp(nans) = [];
    if ~isempty(numobs)
        numobs(nans) = [];
    end
end

% create options for fminsearch
% options = optimset('MaxFunEvals',1000000);
options = optimset('MaxFunEvals',10000000,'MaxIter',10000000);

% create vector with initial values
% b(1) range
% b(2) sill
% b(3) nugget if supplied
b0 = [a0 c0 params.nugget];

% variogram function definitions
switch lower(params.model)    
    case 'spherical'
        type = 'bounded';
        func = @(b,h)b(2)*((3*h./(2*b(1)))-1/2*(h./b(1)).^3);
    case 'pentaspherical'
        type = 'bounded';
        func = @(b,h)b(2)*(15*h./(8*b(1))-5/4*(h./b(1)).^3+3/8*(h./b(1)).^5);
    case 'blinear'
        type = 'bounded';
        func = @(b,h)b(2)*(h./b(1));
    case 'circular'
        type = 'bounded';
        func = @(b,h)b(2)*(1-(2./pi)*acos(h./b(1))+2*h/(pi*b(1)).*sqrt(1-(h.^2)/(b(1)^2)));
    case 'exponential'
        type = 'unbounded';
        func = @(b,h)b(2)*(1-exp(-h./b(1)));
    case 'gaussian'
        type = 'unbounded';
        func = @(b,h)b(2)*(1-exp(-(h.^2)/(b(1)^2)));
    case 'stable'
        type = 'unbounded';
        func = @(b,h)b(2)*(1-exp(-(h.^params.stablealpha)/(b(1)^params.stablealpha)));  
    case 'whittle'
        type = 'unbounded';
        func = @(b,h)b(2)*(1-h/b(1).*besselk(1,h/b(1)));
    case 'matern'
        type = 'unbounded';
        func = @(b,h)b(2)*(1-(1/((2^(params.nu-1))*gamma(params.nu))) * (h/b(1)).^params.nu .* besselk(params.nu,h/b(1)));
    otherwise
        error('unknown model')
end


% check if there are zero distances 
% if yes, remove them, since the besselk function returns nan for
% zero
switch lower(params.model) 
    case {'whittle','matern'}
        izero = h==0;
        if any(izero)
            flagzerodistances = true;
        else
            flagzerodistances = false;
        end
    otherwise
        flagzerodistances = false;
end
        

% if model type is unbounded, then the parameter b(1) is r, which is
% approximately range/3. 
switch type
    case 'unbounded'
        b0(1) = b0(1)/3;
end


% nugget variance
if isempty(params.nugget)
    nugget = false;
    funnugget = @(b) 0;
else
    nugget = true;
    funnugget = @(b) b(3);
end

% generate upper and lower bounds when fminsearchbnd is used
switch lower(params.solver)
    case {'fminsearchbnd'};
        % lower bounds
        lb = zeros(size(b0));
        % upper bounds
        if nugget;
            ub = [inf max(gammaexp) max(gammaexp)]; %
        else
            ub = [inf max(gammaexp)];
        end
end

% create weights (see Webster and Oliver)
if isempty(numobs);
    weights = @(b,h) 1;
else
    switch params.weightfun
        case 'cressie85'
            weights = @(b,h) (numobs./variofun(b,h).^2)./sum(numobs./variofun(b,h).^2); 
        case 'mcbratney86'
            weights = @(b,h) (numobs.*gammaexp./variofun(b,h).^3)/sum(numobs.*gammaexp./variofun(b,h).^3);
        otherwise
            weights = @(b,h) 1;
    end
end

% create objective function: weighted least squares
if isempty(numobs);
    objectfun = @(b)sum(((variofun(b,h)-gammaexp).^2).*weights(b,h));
else
    objectfun = @(b)sum((variofun(b,h)-gammaexp).^2);
end

% call solver
switch lower(params.solver)
    case 'fminsearch'                
        % call fminsearch
        [b,fval,exitflag,output] = fminsearch(objectfun,b0,options);
    case 'fminsearchbnd'
        % call fminsearchbnd
        [b,fval,exitflag,output] = fminsearchbnd(objectfun,b0,lb,ub,options);
    otherwise
        error('Variogramfit:Solver','unknown or unsupported solver')
end


% prepare output
a = b(1); %range
c = b(2); %sill
if nugget;
    n = b(3);%nugget
else
    n = [];
end


% Create structure array with results 
if nargout == 4;    
    S.model     = params.model; % model
    S.h         = h; % distance
    S.gamma     = gammaexp; % experimental values
    S.gammahat  = variofun(b,h); % estimated values
    S.residuals = gammaexp-S.gammahat; % residuals
    COVyhaty    = cov(S.gammahat,gammaexp);
    S.Rs        = (COVyhaty(2).^2) ./...
                  (var(S.gammahat).*var(gammaexp)); % Rsquare
    S.weights   = weights(b,h); %weights
    S.exitflag  = exitflag; % exitflag (see doc fminsearch)
    S.algorithm = output.algorithm;
    S.funcCount = output.funcCount;
    S.iterations= output.iterations;
    S.message   = output.message;
end



% if you want to plot the results...
if params.plotit
    switch lower(type)
        case 'bounded'
            plot(h,gammaexp,'rs','MarkerSize',10);
            hold on
            fplot(@(h) funnugget(b) + func(b,h),[0 b(1)])
            fplot(@(h) funnugget(b) + b(2),[b(1) max(h)])
            
        case 'unbounded'
            plot(h,gammaexp,'rs','MarkerSize',10);
            hold on
            fplot(@(h) funnugget(b) + func(b,h),[0 max(h)])
    end
    axis([0 max(h) 0 max(gammaexp)])
    xlabel('lag distance h')
    ylabel('\gamma(h)')
    hold off
end


% fitting functions for  fminsearch/bnd
function gammahat = variofun(b,h)
    
    switch type
        % bounded model
        case 'bounded'
            I = h<=b(1);
            gammahat     = zeros(size(I));
            gammahat(I)  = funnugget(b) + func(b,h(I));
            gammahat(~I) = funnugget(b) + b(2);        
        % unbounded model
        case 'unbounded'
            gammahat = funnugget(b) + func(b,h);
            if flagzerodistances
                gammahat(izero) = funnugget(b);
            end    
    end
end

end

% and that's it...





