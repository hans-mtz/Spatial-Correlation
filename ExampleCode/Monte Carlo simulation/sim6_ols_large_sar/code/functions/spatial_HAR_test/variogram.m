function S = variogram(x,y,varargin)

% isotropic and anisotropic experimental (semi-)variogram
%
% Syntax:
%   d = variogram(x,y)
%   d = variogram(x,y,'propertyname','propertyvalue',...)
%
% Description:
%   variogram calculates the experimental variogram in various 
%   dimensions. It makes extensive use of a function written
%   by John D'Errico (IPDM). IPDM is required to run the function.
%
% Input:
%   x - array with coordinates. Each row is a location in a 
%       size(x,2)-dimensional space (e.g. [x y elevation])
%   y - column vector with values of the locations in x. 
%
% Propertyname/-value pairs:
%   nrbins - number bins the distance should be grouped into
%            (default = 20)
%   maxdist - maximum distance for variogram calculation
%            (default = maximum distance in the dataset / 2)
%   type -   'gamma' returns the variogram value (default)
%            'cloud1' returns the binned variogram cloud
%            'cloud2' returns the variogram cloud
%   plotit - true -> plot variogram
%            false -> don't plot (default)
%   anisotropy - false (default), true (works only in two dimensions)
%   thetastep - if anisotropy is set to true, specifying thetastep 
%            allows you the angle width (default 30�)
%   
%   
% Output:
%   d - structure array with distance and gamma - vector
%   
% Example: Generate a random field with periodic variation in x direction
% 
%     x = rand(1000,1)*4-2;  
%     y = rand(1000,1)*4-2;
%     z = 3*sin(x*15)+ randn(size(x));
%
%     subplot(2,2,1)
%     scatter(x,y,4,z,'filled'); box on;
%     ylabel('y'); xlabel('x')
%     title('data (coloring according to z-value)')
%     subplot(2,2,2)
%     hist(z,20)
%     ylabel('frequency'); xlabel('z')
%     title('histogram of z-values')
%     subplot(2,2,3)
%     d = variogram([x y],z,'plotit',true,'nrbins',50);
%     title('Isotropic variogram')
%     subplot(2,2,4)
%     d2 = variogram([x y],z,'plotit',true,'nrbins',50,'anisotropy',true);
%     title('Anisotropic variogram')
%
% Requirements:
%   The function requires IPDM written by John D'Errico
%   available on the Mathworks FEX (objectId=18937)
%   http://www.mathworks.com/matlabcentral/
%   The function uses parseargs (objectId=10670) 
%   by Malcolm wood as subfunction.
%
% See also: IPDM
%
% Date: 8. February, 2010
% Author: Wolfgang Schwanghart (w.schwanghart[at]unibas.ch)


% error checking
if size(y,1) ~= size(x,1);
    error('x and y must have the same number of rows')
end
% check if consolidator and ipdm are available
if exist('ipdm.m','file') ~= 2;
    error('IPMD is not available. See help variogram for more infos')
end

% check for nans
II = any(isnan(x),2) | isnan(y);
x(II,:) = [];
y(II)   = [];


% extent of dataset
minx = min(x,[],1);
maxx = max(x,[],1);
maxd = sqrt(sum((maxx-minx).^2));
nrdims = size(x,2);

% check input using PARSEARGS
params.nrbins      = 20;
params.maxdist     = maxd/2;
params.type        = {'default','gamma','cloud1','cloud2'};
params.plotit      = false;
params.anisotropy  = false;
params.thetastep   = 30;
params = parseargs(params,varargin{:});

if params.maxdist > maxd;
    warning('Matlab:Variogram',...
            ['Maximum distance exceeds maximum distance \n' ... 
             'in the dataset. maxdist was decreased to ' num2str(maxd) ]);
    params.maxdist  = maxd;
end

if params.anisotropy && nrdims ~= 2 
    params.anisotropy = false;
    warning('Matlab:Variogram',...
            'Anistropy is only supported for 2D data');
end

% calculate bin tolerance
tol      = params.maxdist/params.nrbins;


% calculate euclidean interpoint distances using ipdm
d = ipdm(x,'Result','Structure',...
           'Subset','Maximum',...
           'Limit',params.maxdist);

% remove distances were d.columnindex = d.rowindex
iid = [d.rowindex d.columnindex d.distance];

% clear workspace variable
clear d

% remove double entries in iid
iid(iid(:,1) == iid(:,2),:) = []; 
[m,m] = unique(sort(iid(:,[1 2]),2),'rows');
iid = iid(m,:);

% calculate squared difference between values of coordinate pairs
lam      = (y(iid(:,1))-y(iid(:,2))).^2;

% anisotropy
if params.anisotropy 
    nrthetaedges = floor(180/params.thetastep);
  
    % calculate with radians, not degrees
    params.thetastep = params.thetastep/180*pi;

    % calculate angles, note that angle is calculated clockwise from top
    theta    = atan2(x(iid(:,2),1)-x(iid(:,1),1),...
                     x(iid(:,2),2)-x(iid(:,1),2));
    
    % only the semicircle is necessary for the directions
    I        = theta < 0;
    theta(I) = theta(I)+pi;
    I        = theta >= pi-params.thetastep/2;
    theta(I) = 0;
        
    % create a vector with edges for binning of theta
    % directions go from 0 to 180 degrees;
    thetaedges = linspace(-params.thetastep/2,pi-params.thetastep/2,nrthetaedges);
    
    % bin theta
    [ntheta,ixtheta] = histc(theta,thetaedges);
    
    % bin centers
    thetacents = thetaedges(1:end)+params.thetastep/2;
    thetacents(end) = pi; %[];
end

% calculate variogram
switch params.type
    case {'default','gamma'}
        % variogram anonymous function
        fvar     = @(x) 1./(2*numel(x)) * sum(x);
        
        % distance bins
        edges      = linspace(0,params.maxdist,params.nrbins+1);
        edges(end) = inf;

        [nedge,ixedge] = histc(iid(:,3),edges);
        
        if params.anisotropy
            S.val      = accumarray([ixedge ixtheta],lam,...
                                 [numel(edges) numel(thetaedges)],fvar,nan);
            S.val(:,end)=S.val(:,1); 
            S.theta    = thetacents;
            S.num      = accumarray([ixedge ixtheta],ones(size(lam)),...
                                 [numel(edges) numel(thetaedges)],@sum,nan);
            S.num(:,end)=S.num(:,1);                 
        else
            S.val      = accumarray(ixedge,lam,[numel(edges) 1],fvar,nan);     
            S.num      = accumarray(ixedge,ones(size(lam)),[numel(edges) 1],@sum,nan);
        end
        S.distance = (edges(1:end-1)+tol/2)';
        S.val(end,:) = [];
        S.num(end,:) = [];

    case 'cloud1'
        edges      = linspace(0,params.maxdist,params.nrbins+1);
        edges(end) = inf;
        
        [nedge,ixedge] = histc(iid(:,3),edges);
        
        S.distance = edges(ixedge);
        S.val      = lam;  
        if params.anisotropy            
            S.theta   = thetacents(ixtheta);
        end
    case 'cloud2'
        S.distance = iid(:,3);
        S.val      = lam;
        if params.anisotropy            
            S.theta   = thetacents(ixtheta);
        end
end


% create plot if desired
if params.plotit
    switch params.type
        case {'default','gamma'}
            marker = 'o--';
        otherwise
            marker = '.';
    end
    
    if ~params.anisotropy
        plot(S.distance,S.val,marker);
        axis([0 params.maxdist 0 max(S.val)*1.1]);
        xlabel('h');
        ylabel('\gamma (h)');
        title('(Semi-)Variogram');
    else
        [Xi,Yi] = pol2cart(repmat(S.theta,numel(S.distance),1),repmat(S.distance,1,numel(S.theta)));
        surf(Xi,Yi,S.val)
        xlabel('h y-direction')
        ylabel('h x-direction')
        zlabel('\gamma (h)')
        title('directional variogram')
%         set(gca,'DataAspectRatio',[1 1 1/30])
    end
end
        
end













