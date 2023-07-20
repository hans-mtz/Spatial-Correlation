

clear
tic


%% LOAD SETTINGS
B = 1000; % number of iterations
Bboot = 1000; % number of iterations in parametric bootstrap
modelName = 'iv'; 
sampleScale = 'small'; % small = 205 locations, large = 820 locations
spatialModel = 'baseline'; % true DGP of underlying dependence
settingName = 'iv_small_baseline';
testMode = 0; % set to 1 for faster dugging


%% INITIALIZATION
initialization


%% MAIN LOOP
WaitMessage = parfor_wait(B);
parfor bb = 1 : B

WaitMessage.Send
rng(bb)


%% DGP

% Generate error for this realization
if lower(spatialModel) == "sar"
    [U,V] = dgp_sar_iv(N,T,dis_mat_unique);
else
    err = CSigma'*mvnrnd([0 0],[1 .8 ; .8 1],n);
    U = err(:,1);
    V = err(:,2);
    if lower(spatialModel) == "hetero"
        U = hetero.*U;
        V = hetero.*V;
    end
end

% Generate outcome for this realization
D = Z*pi0 + X*betaD + V;
Y = D*alpha0 + X*betaY + U;  

% First stage estimates (full sample)
fscoef = [Z X]\D;
vhat = D - [Z X]*fscoef;

% IV estimate
mD = D - X*(X\D);
mY = Y - X*(X\Y);
ahat(bb) = (mZ'*mD)\(mZ'*mY);
uhat = mY-mD*ahat(bb); % residual

betaY_hat = ([Z X]'*[D X])\([Z X]'*Y);
betaY_hat(1) = [];


%% FIXED-G GROUP-BASED METHODS

% Clustered by location
seLocation(bb) = cluster_se(mZ,uhat,1/(mZ'*mD),group_location,p+2);

seG_cce_temp = zeros(1,l_G);
ahatG_fm_temp = zeros(1,l_G);
seG_im_temp = zeros(1,l_G);
pValNullG_crs_temp = zeros(1,l_G);
for kk = 1 : l_G
    % full sample - CCE
    seG_cce_temp(kk) = cluster_se(mZ,uhat,1/(mZ'*mD),...
        group_matrix_km(:,kk),p+2);
    
    % Fama-MacBeth - IM & CRS
    [bfmtmp,sfmtmp,beta_g] = FamaMacbeth(D,X,Y,Z,group_matrix_km(:,kk)); 
    ahatG_fm_temp(kk) = bfmtmp(1);
    seG_im_temp(kk) = sfmtmp(1); % IM se
    pValNullG_crs_temp(kk) = CRSpVal(beta_g(~isnan(...
        beta_g(:,1)),1)); 
    
end
seG_cce(bb,:) = seG_cce_temp;
ahatG_fm(bb,:) = ahatG_fm_temp;
seG_im(bb,:) = seG_im_temp;
pValNullG_crs(bb,:) = pValNullG_crs_temp;


%% COVARIANCE ESTIMATION

Sigma_func = @(w) M(useQML,:)*(exp(w(1))*exp(-dis_mat/exp(w(2))-...
    time_mat/exp(w(3))))*(M(useQML,:)');
QU = @(w) .5*logdet(Sigma_func(w))+.5*uhat(useQML,1)'*...
    (Sigma_func(w)\uhat(useQML,1));
QV = @(w) .5*logdet(Sigma_func(w))+.5*vhat(useQML,1)'*...
    (Sigma_func(w)\vhat(useQML,1));
options = optimoptions('fminunc','display','off');
alpha_init = zeros(3,1);
alphaHatU = fminunc(QU,alpha_init,options);
alphaHatV = fminunc(QV,alpha_init,options);

Sigma_func_DGP = @(w) exp(w(1))*exp(-dis_mat/exp(w(2))-...
    time_mat/exp(w(3)));
SigmaHatU = Sigma_func_DGP(alphaHatU);
SigmaHatV = Sigma_func_DGP(alphaHatV);
    
if min(eig(SigmaHatU)) < 0 || min(eig(SigmaHatV)) < 0
    error('Non psd covariance matrix');
end
        

%% PARAMETRIC BOOTSTRAP

CSHatU = sqrtm(SigmaHatU);
CSHatV = sqrtm(SigmaHatV);
rhohat = corr(CSHatU\uhat,CSHatV\vhat);
    
% check criterion
qmle(bb,:) = [exp(alphaHatU)' exp(alphaHatV)' rhohat] ;

aboot_cce = zeros(Bboot,1);
sloc = zeros(Bboot,1);
sboot_cce = zeros(Bboot,l_G);
aboot_fm = zeros(Bboot,l_G);
sboot_fm = zeros(Bboot,l_G);
tboot_fm = zeros(Bboot,l_G);

pValTHboot = zeros(Bboot,l_G);
pValTHboot_alt = zeros(Bboot,nAltPowerSim,l_G);
pValTHboot_test = zeros(Bboot,l_G);

for rr = 1 : Bboot
    
    % DGP
    errboot = mvnrnd([0 0],[1 rhohat ; rhohat 1],n);
    Uboot = CSHatU'*errboot(:,1);
    Vboot = CSHatV'*errboot(:,2);        
    Dboot = [Z X]*fscoef + Vboot;
    Yboot = X*betaY_hat+Uboot;
    mDboot = Dboot - X*(X\Dboot);
	mYboot = Yboot - X*(X\Yboot);
    
    % full-sample estimation
    aboot_cce(rr) = (mZ'*mDboot)\(mZ'*mYboot);
    residboot = mYboot-mDboot*aboot_cce(rr);

    % CCE by location
    sloc(rr,1) = cluster_se(mZ,residboot,1/(mZ'*mDboot),...
        group_location,p+2);
 
    for kk = 1:l_G
        
        group = group_matrix_km(:,kk);
        
        % CCE by k-medoids groupings
        sboot_cce(rr,kk) = cluster_se(mZ,residboot,1/(mZ'*mDboot),...
            group_matrix_km(:,kk),p+2);
        
        % IM
        [bfmtmp,sfmtmp,beta_g] = FamaMacbeth(Dboot,X,Yboot,Z,group); 
        aboot_fm(rr,kk) = bfmtmp(1);
        sboot_fm(rr,kk) = sfmtmp(1);
        tboot_fm(rr,kk) = bfmtmp(1)/sfmtmp(1);
        
        % CRS
        pValTHboot(rr,kk) = CRSpVal(beta_g(~isnan(...
        	beta_g(:,1)),1)); 
        % for power simulation
        for aa = 1 : nAltPowerSim
            Yboot_alt = Dboot*altPowerSim(aa)+X*betaY+Uboot; 
            [~,~,beta_g] = FamaMacbeth(Dboot,X,Yboot_alt,Z,group);
            pValTHboot_alt(rr,aa,kk) = CRSpVal(beta_g(~isnan(...
                beta_g(:,1)),1)); 
        end
    end
    
end


%% POWER OPTIMIZATION SUBJECT TO SIZE CONTROL

wap_cce = zeros(l_G+1,1);
wap_im = zeros(l_G,1);

% simulated power for clustering by location
cvalLocation(bb) = max(tinv(1-sigLevel/2,max(group_location)-1),...
    prctile(abs(aboot_cce./sloc),(1-sigLevel)*100));
wap_cce(1) = mean(mean(abs(aboot_cce*ones(1,nAltPowerSim)-ones(Bboot,1)*...
    altPowerSim')./(sloc*ones(1,nAltPowerSim)) > cvalLocation(bb)));

cvalG_cce_temp = zeros(1,l_G);
cvalG_im_temp = zeros(1,l_G);
for kk = 1:l_G
    
    group = group_matrix_km(:,kk);
    
    % simulated power for clustering by k-medoids groupings
    cvalG_cce_temp(kk) = max(tinv(1-sigLevel/2,max(group)-1),...
    	prctile(abs(aboot_cce./sboot_cce(:,kk)),(1-sigLevel)*100));
    wap_cce(1+kk) = mean(mean(abs(aboot_cce*ones(1,nAltPowerSim)-...
        ones(Bboot,1)*altPowerSim')./(sboot_cce(:,kk)*ones(1,...
        nAltPowerSim)) > cvalG_cce_temp(kk)));
    
    % IM
    cvalG_im_temp(kk) = max(tinv(1-sigLevel/2,...
        max(group)-1),prctile(abs(tboot_fm(:,kk)),...
        (1-sigLevel)*100));
    wap_im(kk) = mean(mean(abs(aboot_fm(:,kk)*ones(1,nAltPowerSim)-...
        ones(Bboot,1)*altPowerSim')./(sboot_fm(:,kk)*...
        ones(1,nAltPowerSim)) > cvalG_im_temp(kk)));
    
end
    
cvalG_cce(bb,:) = cvalG_cce_temp;
cvalG_im(bb,:) = cvalG_im_temp;

% optimization including clustering by location
[~,indA_cce] = max(wap_cce);
if indA_cce == 1
    seGstarA(bb) = seLocation(bb);
    cvalGstarA(bb) = cvalLocation(bb);
    GstarA(bb) = max(group_location);        
else
    seGstarA(bb) = seG_cce_temp(indA_cce-1);
    cvalGstarA(bb) = cvalG_cce_temp(indA_cce-1);
    GstarA(bb) = max(group_matrix_km(:,indA_cce-1));
end

% optimization including only k-medoids groupings
[~,indB_cce] = max(wap_cce(2:end));
seGstarB(bb) = seG_cce_temp(indB_cce);
cvalGstarB(bb) = cvalG_cce_temp(indB_cce);
GstarB(bb) = max(group_matrix_km(:,indB_cce));

% IM optimization 
[~,ind_im] = max(wap_im);
aGstar_im(bb,1) = ahatG_fm_temp(ind_im);
seGstar_im(bb,1) = seG_im_temp(ind_im);
cvalGstar_im(bb,1) = cvalG_im_temp(ind_im);
Gstar_im(bb,1) = max(group_matrix_km(:,ind_im));

% CRS optimization
pValAdjG_temp = min(quantile(pValTHboot,sigLevel)-eps,sigLevel);
wap_crs = zeros(5,1);
for kk = 1 : l_G
    wap_crs(kk) = mean(mean(pValTHboot_alt(:,:,kk) <= pValAdjG_temp(kk)));
end
[~,ind_crs] = max(wap_crs);
if wap_crs(ind_crs) == 0
    ind_crs = l_G;
end
pValAdjG_crs(bb,:) = pValAdjG_temp;
pValGstar_crs(bb,1) = pValNullG_crs_temp(ind_crs);
pvalAdjGstar_crs(bb,1) = pValAdjG_temp(ind_crs);
Gstar_crs(bb,1) = max(group_matrix_km(:,ind_crs));
aGstar_crs(bb,1) = ahatG_fm_temp(ind_crs);

% p-values for CRS power curve plot
pValPowerCurveG_crs_temp = zeros(1,nAltPowerCurve,l_G);
pValPowerCurveGstar_crs_temp = zeros(1,nAltPowerCurve);
for aa = 1 : nAltPowerCurve

    Yboot_alt = D*altPowerCurve(aa)+X*betaY+U;  
    
    for kk = 1 : l_G
        group = group_matrix_km(:,kk);
        [~,~,beta_g] = FamaMacbeth(D,X,Yboot_alt,Z,group); % FM results
        pValPowerCurveG_crs_temp(1,aa,kk) = CRSpVal(beta_g(~isnan(...
            beta_g(:,1)),1)); 

    end
    
    
    pValPowerCurveGstar_crs_temp(1,aa) = pValPowerCurveG_crs_temp(...
        1,aa,ind_crs);
    
end
pValPowerCurveG_crs(bb,:,:) = pValPowerCurveG_crs_temp;
pValPowerCurveGstar_crs(bb,:) = pValPowerCurveGstar_crs_temp;


%% Sun and Kim (2015) with optimal bandwidth by Lazarus et al. (2018) 

% parameters
R = 1;
r = 0;
[xloc, yloc] = form_regular_lattice(Xcoord,Ycoord);
loc = [xloc yloc];
q = size(R,1); % number of restrictions
N_loc = size(xloc,1);
K_rule = .4*(N_loc/2)^(2/3);
% K1 = round(sqrt(K_rule+1)-1);
K1 = round(sqrt(K_rule));
K2 = K1;
K = K1 + K2 + K1*K2;
IROR_trans = F_star_IROR_trans_iv(R,r,mY,mD,mZ,loc,K1,K2);

% p-value under the null
theta_hat = ahat(bb);
F_trans = (2*K-q+1)/(2*K)*(sqrt(N_loc)*(R*theta_hat-r))'*IROR_trans*...
    (sqrt(N_loc)*(R*theta_hat-r))/q; 
pValNull_sk(bb) = 1-fcdf(F_trans,q,2*K-q+1);

% p-value under alternatives
pValAlt_sk_temp = zeros(1,nAltPowerCurve);
for aa = 1 : nAltPowerCurve
    theta_hat = ahat(bb)+altPowerCurve(aa);
    F_trans = (2*K-q+1)/(2*K)*(sqrt(N_loc)*(R*theta_hat-r))'*IROR_trans*...
        (sqrt(N_loc)*(R*theta_hat-r))/q; 
    pValAlt_sk_temp(aa) =  1-fcdf(F_trans,q,2*K-q+1);
end
pValAlt_sk(bb,:) = pValAlt_sk_temp;



end
save('../temp/all_results.mat')


%% OUTPUT RESULTS
output


toc
