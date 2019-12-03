%% UPENN, 714, Prof Dirk Krueger, Problem set 02.
% Rodrigo Morales
% based on a code by Anna Cororaton
% November 2019

%Housekeeping
clear; clc; close all;

%% parameters
r               = 0.02;         % interest rate
delta           = 0.8;          % depreciation
rho             = 0.04;         % bbeta = 1 / (1+rho)
bbeta           = 1/(1+rho);    
sigmaepsilon    = 0.2;          % TFP std dev
sigma           = 1;            % Elasticity of subs (utility).
if sigma ==1                    % if sigma ==1, use log utility
    u=@(c) log(c);
else                            % if sigma ~=1, use CRRA
    u=@(c) c.^(1-sigma)./(1-sigma);
end

if delta < 0.8          % if shock is not persistent, do Tauchen
    doRouenhorst = 0; 
else                    % if shock is persistent, do Rouenhorst
    doRouenhorst =1;
end

% CALIBRATE: !!!
T       = 60;     % periods for each simulation (1000 or 5000)
numsim  = 1000;
nk      = 2000;      % number of grid points for K
na      = 11;       % number of grid points for A
N       = na;       % change of name
tol     = 10e-10;   % tolerance level for VFI
maxiter = 1000;     % maximum number of iterations
d       = 100;      % distance metric

%% Get transition matrix (Tauchen or Rouenhorst)

if doRouenhorst == 0        	% do Tauchen
    [a5,ap5]    = tauchen_ram(5,delta,sigmaepsilon,3);
    [a9,ap9]    = tauchen_ram(N,delta,sigmaepsilon,3);
    [a15,ap15]  = tauchen_ram(15,delta,sigmaepsilon,3);
    [vProductivity,mTransition]  = tauchen_ram(N,delta,sigmaepsilon,3);
    vProductivity = vProductivity';
    % previous is equivalent to:  (slight variation because of sensibility to
    % calculation. Bottom line...   works :) ,   but is senstitive :(
    %[Z,Zprob] = tauchen(N,0,delta,(1-delta^2)^(1/2)*sigmaepsilon,3);
    % calculate long run distribution
    a5star      = pstar(ap5);
    a9star      = pstar(ap9);
    a15star     = pstar(ap15);
    prodStatnry = pstar(mTransition);
else   %doRouenhorst == 1;      % do Rouenhorst
    p = (1+delta)/2;
    Psi = sigmaepsilon*sqrt(N-1);
    [a5,ap5]    = rouwenhorst_ram(5,p,p,Psi);
    [a9,ap9]    = rouwenhorst_ram(N,p,p,Psi);
    [a15,ap15]  = rouwenhorst_ram(15,p,p,Psi);
    [vProductivity,mTransition] = rouwenhorst_ram(N,p,p,Psi);
    vProductivity = vProductivity';
    % calculate long run distribution
    a5star      = pstar(ap5);
    a9star      = pstar(ap9);
    a15star     = pstar(ap15);
    prodStatnry = pstar(mTransition);
end

%% create grid for k and a

a       = exp(a9);
% savings 
savmin    = 0.001;
savmax    = vProductivity(N)/(1-bbeta); % -natural borrowing constraint
k = curvspace(savmin,savmax,nk,2)'; % use curved grid to enhance accuracy
%Vinit       = repmat(k,1,na);     % value function init...
%Vinit       = repmat(k*0,1,na);     % value function init...
Vinit = u( (1+r)*repmat(k,1,na)+repmat(a,nk,1) );

%% Value Function Iteration...
%2.2
tic;
[V,IG,S,C] = vfi_01_infty(k,a,ap9,bbeta,r,u,d,tol,maxiter,Vinit);
%01 uses the monotonicity and convexity
%[V,IG,S,C] = vfi_02_infty(k,a,ap9,bbeta,r,u,d,tol,maxiter,Vinit);
% 02 maximizes over whole grid.
toc;
%
% plots
figure;
plot(k,V(:,[1 5 9])); 
title('Value function');
legend('Min shock','Mean shock','Max shock','Location','NorthWest');
xlabel('Assets'); 
ylabel('Value');

%saveas(gcf,'fig01_1_val.eps','epsc2');
%
figure;
plot(k,S(:,[1 5 9])); 
title('Savings Policy function');
xlabel('Assets'); 
ylabel('Policy');
legend('Min shock','Mean shock','Max shock','Location','NorthWest');
%saveas(gcf,'fig01_2_assets.eps','epsc2');

figure;
plot(k,C(:,[1 5 9])); 
title('Consumption Policy function');
xlabel('Assets'); 
ylabel('Policy');
legend('Min shock','Mean shock','Max shock','Location','NorthWest');
%saveas(gcf,'fig01_3_cons.eps','epsc2');


%% Simulation for first 61 periods of an agent's life:

mc = dtmc(mTransition);
shocks_sim1 = zeros(nsim,T);    % contains all the shock for  simulations
shocks_sim1i = zeros(nsim,T);   % contains all the loc shocks 4 simltns
Assets_sim1 = zeros(nsim,T+1);  % contains all the assets simulations
Assets_sim1i = zeros(nsim,T+1); % contains all the loc assets simltns
CS_sim1 = zeros(nsim,T);        % contains all the consumptions simulation

for iisim = 1:numsim
    
    k0i = randi([1 nk],1,1);         % take the initial asset level random;
    shocksT = simulate(mc,T);        % productivity shocks
    shocks_sim1(iisim,:) = exp(vProductivity(shocksT')); % recrd shocks sim
    shocks_sim1i(iisim,:) = shocksT';       % record shocks loc from simln
    Assets_sim1(iisim,1) = k(k0i);          % record asset holding...
    Assets_sim1i(iisim,tt+1) = k0i;         % record asset loc sim
    
    for tt = 1:T
        prodShock = exp(vProductivity(shocksT(tt)));
        prodShocki = shocksT(tt);
        kprime = S(k0i,prodShocki);
        consumpt = C(k0i,prodShocki);
        CS_sim1(iisim,tt) = consumpt;           % record consumption sim
        k0i = IG(k0i,prodShocki);               % get new loc for savings
        Assets_sim1(iisim,tt+1) = kprime;       % record asset sim
        Assets_sim1i(iisim,tt+1) = k0i;         % record asset loc sim
        CS_sim1(iisim,tt) = consumpt;
        
        
        invest = investmentFunction(k,kprime);
        Q = value0(k0i,b0i,prodShocki,zminusi)/k;
        profit_i = profitFunction(prodShock,k);
        Xs(ii,1) = Q; %Q
        Xs(ii,2) = profit_i/k; %profit_i/k
        ys(ii) = invest/k; % i/k
        Xs2(ii,1) = Q; %Q
        Xs2(ii,2) = profit_i/k; %profit_i/k
        Xs2(ii,3) = log(k); %log(k)
        ys2(ii) = b/k; % b/k
        %find next period values:
        k0i = kPolicyIndex(k0i,b0i,prodShocki,zminusi);
        b0i = bPolicyIndex(k0i,b0i,prodShocki,zminusi);
        zminusi = prodShocki;
        
        k0i = IG(k0i,prodShocki);
        
        
        k = kprime;
        b = bprime;
    end
    
    
end


for ii = 1:numsteps
    prodShock = exp(grid_a_log(Zetasii(ii)));
    prodShocki = Zetasii(ii);
    kprime = kPolicy(k0i,b0i,prodShocki,zminusi);
    bprime = bPolicy(k0i,b0i,prodShocki,zminusi);
    invest = investmentFunction(k,kprime);
    Q = value0(k0i,b0i,prodShocki,zminusi)/k;
    profit_i = profitFunction(prodShock,k);
    Xs(ii,1) = Q; %Q
    Xs(ii,2) = profit_i/k; %profit_i/k
    ys(ii) = invest/k; % i/k
    Xs2(ii,1) = Q; %Q
    Xs2(ii,2) = profit_i/k; %profit_i/k
    Xs2(ii,3) = log(k); %log(k)
    ys2(ii) = b/k; % b/k
    %find next period values:
    k0i = kPolicyIndex(k0i,b0i,prodShocki,zminusi);
    b0i = bPolicyIndex(k0i,b0i,prodShocki,zminusi);
    zminusi = prodShocki;
    k = kprime;
    b = bprime;
end
%Burn some values (first 20%):
Xs = Xs((0.2*numsteps):end,:);
ys = ys((0.2*numsteps):end,:);


