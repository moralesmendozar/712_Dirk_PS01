%% UPENN, 714, Prof Dirk Krueger, Problem set 01.
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
sigmaepsilon    = 0.2;          % TFP std dev in PSet, called sigma_y
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

% CALIBRATEION...
T       = 60;       % periods for each simulation
numsim  = 1000;     % number of simulations
nk      = 2000;     % number of grid points for K
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
    % previous is equivalent to: (slight variation becauseof sensibility to
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
% start in the utility of cash at hand:
Vinit = u( (1+r)*repmat(k,1,na)+repmat(a,nk,1) );

%% (Infinite) Value Function Iteration...
%2.2
tic;
[V,IG,S,C] = vfi_01_infty(k,a,ap9,bbeta,r,u,d,tol,maxiter,Vinit);
%01 uses the monotonicity and convexity
%[V,IG,S,C] = vfi_02_infty(k,a,ap9,bbeta,r,u,d,tol,maxiter,Vinit);
% 02 maximizes over whole grid.
toc;

% plots
figure;
plot(k,V(:,[1 5 9])); 
title('2.2) Value function for infinitely lived agent');
legend('Min shock','Mean shock','Max shock','Location','NorthWest');
xlabel('Assets'); 
ylabel('Value');

%saveas(gcf,'fig01_1_val.eps','epsc2');
%
figure;
plot(k,S(:,[1 5 9])); 
title('2.2) Savings Policy function for infinitely lived agent');
xlabel('Assets'); 
ylabel('Policy');
legend('Min shock','Mean shock','Max shock','Location','NorthWest');
%saveas(gcf,'fig01_2_assets.eps','epsc2');

figure;
plot(k,C(:,[1 5 9])); 
title('2.2) Consumption Policy function for infinitely lived agent');
xlabel('Assets'); 
ylabel('Policy');
legend('Min shock','Mean shock','Max shock','Location','NorthWest');
%saveas(gcf,'fig01_3_cons.eps','epsc2');


%% Simulation for first 61 periods of an infinitely lived agent's life:
% https://www.mathworks.com/help/econ/dtmc.simulate.html#d117e311371
numsim = 1000;
mc = dtmc(mTransition);
shocks_sim1 = zeros(numsim,T+1);  % contains all the shock for  simulations
shocks_sim1i = zeros(numsim,T+1); % contains all the loc shocks 4 simltns
Assets_sim1 = zeros(numsim,T+2);  % contains all the assets simulations
Assets_sim1i = zeros(numsim,T+2); % contains all the loc assets simltns
CS_sim1 = zeros(numsim,T);        % contains all consumptions simulation

for iisim = 1:numsim
    
    %k0i = randi([1 nk],1,1);            % take random initial asset level
    [unusedtemp, k0i] = min(abs(k));     % start simulation with no wealth
    x0 = zeros(1,mc.NumStates);  
    x0(1) = 1;                           % start wimulations w/lowest shock
    shocksT = simulate(mc,T,'X0',x0);    % productivity shocks
    shocks_sim1(iisim,:) = exp(vProductivity(shocksT)'); % recrd shocks sim
    shocks_sim1i(iisim,:) = shocksT';    % record shocks loc from simln
    Assets_sim1(iisim,1) = k(k0i);       % record asset holding...
    Assets_sim1i(iisim,1) = k0i;         % record asset loc sim
    
    for tt = 1:T+1
        prodShocki = shocksT(tt);
        prodShock = exp(vProductivity(prodShocki));
        kprime = S(k0i,prodShocki);
        consumpt = C(k0i,prodShocki);
        CS_sim1(iisim,tt) = consumpt;           % record consumption sim
        k0i = IG(k0i,prodShocki);               % get new loc for savings
        Assets_sim1(iisim,tt+1) = kprime;       % record asset sim
        Assets_sim1i(iisim,tt+1) = k0i;         % record asset loc sim
        
    end
    
end
%%
% Plot of mean of simulation:
figure;
plot(1:T+1,mean( shocks_sim1  ) ); 
title('2.2) Simulation starting with zero wealth and worst income shock');
%title('2.2) Average shocks (income) from simulation');
xlabel('time'); 
%ylabel('Policy (savings)');
%saveas(gcf,'fig01_2_assets.eps','epsc2');

%figure;
hold on;
plot(1:T+1,mean( Assets_sim1(:,2:end)  ) ); 
%title('2.2) Average Policy function (savings) from simulation');
xlabel('time'); 
%ylabel('Policy (savings)');
%saveas(gcf,'fig01_2_assets.eps','epsc2');

%figure;
plot(1:T+1,mean( CS_sim1  ) ); 
%title('2.2) Average Consumption from simulation');
xlabel('time'); 
%ylabel('Policy (savings)');
hold off;
legend('Income','Savings','Consumption','Location','East');
%saveas(gcf,'fig01_3_cons.eps','epsc2');


%% (Finite) Value Function Iteration...
%2.2
disp('Running finitely lived agents VFI')
tic;
[VT,IGT,ST,CT] = vfi_01_finite(k,a,ap9,bbeta,r,u,d,tol,maxiter,Vinit,T);
toc;

% plots for finite stuff
figure;
tplot = 55;
plot(k,VT(:,[1 5 9],tplot));
title('2.3) Value function for dying agent in t = 55');
legend('Min shock','Mean shock','Max shock','Location','NorthWest');
xlabel('Assets'); 
ylabel('Value');

%saveas(gcf,'fig01_1_val.eps','epsc2');
%
figure;
plot(k,ST(:,[1 5 9],tplot)); 
title('2.3) Savings Policy function for dying agent, in t = 55');
xlabel('Assets'); 
ylabel('Policy (savings)');
legend('Min shock','Mean shock','Max shock','Location','NorthWest');
%saveas(gcf,'fig01_2_assets.eps','epsc2');

%% Compare young/old
figure;
plot(k,CT(:,[1 5 9],1)); 
title('2.3) Consmptn Policy function for dying agent, in t = 1 (young)');
xlabel('Assets'); 
ylabel('Policy (consumption)');
legend('Min shock','Mean shock','Max shock','Location','NorthWest');

figure;
plot(k,CT(:,[1 5 9],tplot)); 
title('2.3) Consumption Policy function for dying agent, in t = 55');
xlabel('Assets'); 
ylabel('Policy (consumption)');
legend('Min shock','Mean shock','Max shock','Location','NorthWest');

%saveas(gcf,'fig01_3_cons.eps','epsc2');

%% Simulation for a finitely lived agent
%   Similar to previous, but with Value Function changing...
numsim = 1000;
mc = dtmc(mTransition);
shocks_sim1 = zeros(numsim,T);    % contains all the shock for  simulations
shocks_sim1i = zeros(numsim,T);   % contains all the loc shocks 4 simltns
Assets_sim1 = zeros(numsim,T+1);  % contains all the assets simulations
Assets_sim1i = zeros(numsim,T+1); % contains all the loc assets simltns
CS_sim1 = zeros(numsim,T);        % contains all consumptions simulation

for iisim = 1:numsim
    
    [unusedtemp, k0i] = min(abs(k)); % start simulation with no wealth
    x0 = zeros(1,mc.NumStates);  
    x0(1) = 1;                       % start wimulations with lowest shock
    shocksT = simulate(mc,T-1,'X0',x0);  % productivity shocks
    shocks_sim1(iisim,:) = exp(vProductivity(shocksT)'); % recrd shocks sim
    shocks_sim1i(iisim,:) = shocksT';    % record shocks loc from simln
    Assets_sim1(iisim,1) = k(k0i);       % record asset holding...
    Assets_sim1i(iisim,1) = k0i;         % record asset loc sim
    
    for tt = 1:T
        prodShocki = shocksT(tt);
        prodShock = exp(vProductivity(prodShocki));
        kprime = ST(k0i,prodShocki,tt);
        consumpt = CT(k0i,prodShocki,tt);
        CS_sim1(iisim,tt) = consumpt;           % record consumption sim
        k0i = IGT(k0i,prodShocki,tt);           % get new loc for savings
        Assets_sim1(iisim,tt+1) = kprime;       % record asset sim
        Assets_sim1i(iisim,tt+1) = k0i;         % record asset loc sim
        
    end
    
end

%%
% Plot of mean of simulation:
figure;
plot(1:T,mean( shocks_sim1  ) ); 
title('2.2) Simulation starting with zero wealth and worst income shock');
%title('2.2) Average shocks (income) from simulation');
xlabel('time'); 
%ylabel('Policy (savings)');
%saveas(gcf,'fig01_2_assets.eps','epsc2');

%figure;
hold on;
plot(1:T,mean( Assets_sim1(:,2:end)  ) ); 
%title('2.2) Average Policy function (savings) from simulation');
xlabel('time'); 
%ylabel('Policy (savings)');
%saveas(gcf,'fig01_2_assets.eps','epsc2');

%figure;
plot(1:T,mean( CS_sim1  ) ); 
%title('2.2) Average Consumption from simulation');
xlabel('time'); 
%ylabel('Policy (savings)');
hold off;
legend('Income','Savings','Consumption','Location','East');
%saveas(gcf,'fig01_3_cons.eps','epsc2');
