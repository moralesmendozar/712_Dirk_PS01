%% UPENN, 714, Prof Dirk Krueger, Problem set 01.
% Rodrigo Morales
% Nov/dec 2019

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
savmax    = 18;
k = curvspace(savmin,savmax,nk,2)'; % use curved grid to enhance accuracy
kconstrained = k; %save for later
% start in the utility of cash at hand:
Vinit = u( (1+r)*repmat(k,1,na)+repmat(a,nk,1) );

%% 2.2 (Infinite) Value Function Iteration...
quietdown = 0; %makes the function not print iterations...
tic;
%01 uses the monotonicity and convexity
[V,IG,S,C] = vfi_01_infty(k,a,ap9,bbeta,r,u,d,tol,maxiter,Vinit,quietdown);
%opt is faster, but loses precision
%[V,IG,S,C] = vfi_01_infty_opt(k,a,ap9,bbeta,r,u,d,tol,maxiter,Vinit,quietdown);
% 02 maximizes over whole grid.
%[V,IG,S,C] = vfi_02_infty(k,a,ap9,bbeta,r,u,d,tol,maxiter,Vinit);
toc;
%now compute the environment with sigmaepsilon = 0.4
[asigma4,apsigma4]    = rouwenhorst_ram(N,p,p,0.4*sqrt(N-1));
a2sigma4       = exp(asigma4);
tic;
[V2,IG2,S2,C2] = vfi_01_infty(k,a2sigma4,apsigma4,bbeta,r,u,d,tol,...
    maxiter,Vinit,quietdown);
%[V2,IG2,S2,C2] = vfi_01_infty_opt(k,a2sigma4,apsigma4,bbeta,r,u,d,tol,...
    maxiter,Vinit,quietdown);
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
disp('Starting infinite simulation')
numsim = 1000;
mc = dtmc(mTransition);
shocks_sim1 = zeros(numsim,T+1);  % contains all the shock for  simulations
shocks_sim1i = zeros(numsim,T+1); % contains all the loc shocks 4 simltns
Assets_sim1 = zeros(numsim,T+2);  % contains all the assets simulations
Assets_sim1i = zeros(numsim,T+2); % contains all the loc assets simltns
CS_sim1 = zeros(numsim,T);        % contains all consumptions simulation
%simulation matrices for the sigma 0.4
Assets_sim2 = zeros(numsim,T+2);  % contains all the assets simulations
Assets_sim2i = zeros(numsim,T+2); % contains all the loc assets simltns
CS_sim2 = zeros(numsim,T);        % contains all consumptions simulation


for iisim = 1:numsim
    
    %k0i = randi([1 nk],1,1);            % take random initial asset level
    %%%% CHECK THIS STEP WHEN CONSTRAINED IS RELAXED TO START AT REAL ZERO!
    [~, k0i] = min(abs(k));     % start simulation with no wealth
    k0i2 = k0i;
    [unusedtemp, k0i] = min(abs(k));     % start simulation with no wealth
    k0i2 = k0i;
    x0 = zeros(1,mc.NumStates);  
    x0(1) = 1;                           % start wimulations w/lowest shock
    shocksT = simulate(mc,T,'X0',x0);    % productivity shocks
    shocks_sim1(iisim,:) = exp(vProductivity(shocksT)'); % recrd shocks sim
    shocks_sim1i(iisim,:) = shocksT';    % record shocks loc from simln
    Assets_sim1(iisim,1) = k(k0i);       % record asset holding...
    Assets_sim1i(iisim,1) = k0i;         % record asset loc sim
    Assets_sim2(iisim,1) = k(k0i2);       % record asset holding...
    Assets_sim2i(iisim,1) = k0i2;         % record asset loc sim
    
    for tt = 1:T+1
        prodShocki = shocksT(tt);
        prodShock = exp(vProductivity(prodShocki));
        kprime = S(k0i,prodShocki);
        CS_sim1(iisim,tt) = C(k0i,prodShocki);  % record consumption sim
        k0i = IG(k0i,prodShocki);               % get new loc for savings
        Assets_sim1(iisim,tt+1) = kprime;       % record asset sim
        Assets_sim1i(iisim,tt+1) = k0i;         % record asset loc sim
        kprime2 = S2(k0i2,prodShocki);
        CS_sim2(iisim,tt) = C2(k0i2,prodShocki);  % record consumption sim
        k0i2 = IG2(k0i2,prodShocki);               % get new loc for savings
        Assets_sim2(iisim,tt+1) = kprime2;       % record asset sim
        Assets_sim2i(iisim,tt+1) = k0i2;         % record asset loc sim
        
    end
    
end

disp('done infinite simulation')

%% 2.4) Plot of simulation (infty agnts) consumption with given parameters:
% Plot of mean of simulation:
figure;
plot(1:T+1,nanmean( shocks_sim1  ) ); 
title('2.2) Simulation starting with zero wealth and worst income shock');
xlabel('time'); 

%figure;
hold on;
plot(1:T+1,nanmean( Assets_sim1(:,2:end)  ) ); 
plot(1:T+1,nanmean( CS_sim1  ) ); 

plot(1:T+1,nanmean( Assets_sim2(:,2:end)  ) ); 
plot(1:T+1,nanmean( CS_sim2  ) ); 

hold off;
legend('Income shock','Savings \delta = 0.2','Consumption \delta = 0.2',...
    'Savings \delta = 0.4','Consumption \delta = 0.4','Location','East');
%saveas(gcf,'fig01_4_simInftyAgents.eps','epsc2');




%% 2.3) (Finite) Value Function Iteration...
quietdown = 0;
disp('Running finitely lived agents VFI')
tic;
[VTf,IGTf,STf,CTf] = vfi_01_finite(k,a,ap9,bbeta,r,u,d,tol,maxiter,Vinit,T);
toc;

%now compute the environment with sigmaepsilon = 0.4
[asigma4,apsigma4]    = rouwenhorst_ram(N,p,p,0.4*sqrt(N-1));
a2sigma4       = exp(asigma4);
tic;
[VTf2,IGTf2,STf2,CTf2] = vfi_01_finite(k,a2sigma4,apsigma4,bbeta,r,u,d,...
    tol,maxiter,Vinit,T);
toc;

% plots for finite stuff
figure;
tplot = 55;
plot(k,VTf(:,[1 5 9],tplot));
title('2.3) Value function for dying agent in t = 55');
legend('Min shock','Mean shock','Max shock','Location','NorthWest');
xlabel('Assets'); 
ylabel('Value');

%saveas(gcf,'fig03_1_val.eps','epsc2');
%
figure;
plot(k,STf(:,[1 5 9],tplot)); 
title('2.3) Savings Policy function for dying agent, in t = 55');
xlabel('Assets'); 
ylabel('Policy (savings)');
legend('Min shock','Mean shock','Max shock','Location','NorthWest');
%saveas(gcf,'fig03_2_assets.eps','epsc2');

%% Compare young/old
figure;
plot(k,CTf(:,[1 5 9],1)); 
title('2.3) Consmptn Policy function for dying agent, in t = 1 (young)');
xlabel('Assets'); 
ylabel('Policy (consumption)');
legend('Min shock','Mean shock','Max shock','Location','NorthWest');

figure;
plot(k,CTf(:,[1 5 9],tplot)); 
title('2.3) Consumption Policy function for dying agent, in t = 55');
xlabel('Assets'); 
ylabel('Policy (consumption)');
legend('Min shock','Mean shock','Max shock','Location','NorthWest');

%saveas(gcf,'fig03_3_cons.eps','epsc2');

%% Simulation for a finitely lived agent
%   Similar to previous, but with Value Function changing...
disp('Starting finite simulation')
numsim = 1000;
mc = dtmc(mTransition);
shocks_sim1 = NaN(numsim,T);    % contains all the shock for  simulations
shocks_sim1i = NaN(numsim,T);   % contains all the loc shocks 4 simltns
Assets_sim1 = NaN(numsim,T+1);  % contains all the assets simulations
Assets_sim1i = NaN(numsim,T+1); % contains all the loc assets simltns
CS_sim1 = NaN(numsim,T);        % contains all consumptions simulation
Assets_sim2 = NaN(numsim,T+1);  % contains all the assets simulations
Assets_sim2i = NaN(numsim,T+1); % contains all the loc assets simltns
CS_sim2 = NaN(numsim,T);        % contains all consumptions simulation

for iisim = 1:numsim
    
    [~, k0i] = min(abs(k)); % start simulation with no wealth
    k0i2 = k0i;
    x0 = zeros(1,mc.NumStates);  
    x0(1) = 1;                       % start wimulations with lowest shock
    shocksT = simulate(mc,T-1,'X0',x0);  % productivity shocks
    shocks_sim1(iisim,:) = exp(vProductivity(shocksT)'); % recrd shocks sim
    shocks_sim1i(iisim,:) = shocksT';    % record shocks loc from simln
    Assets_sim1(iisim,1) = k(k0i);       % record asset holding...
    Assets_sim1i(iisim,1) = k0i;         % record asset loc sim
    Assets_sim2(iisim,1) = k(k0i2);       % record asset holding...
    Assets_sim2i(iisim,1) = k0i2;         % record asset loc sim
    
    for tt = 1:T
        prodShocki = shocksT(tt);
        prodShock = exp(vProductivity(prodShocki));
        kprime = STf(k0i,prodShocki,tt);
        CS_sim1(iisim,tt) = CTf(k0i,prodShocki,tt);% record consumption sim
        k0i = IGTf(k0i,prodShocki,tt);           % get new loc for savings
        Assets_sim1(iisim,tt+1) = kprime;       % record asset sim
        Assets_sim1i(iisim,tt+1) = k0i;         % record asset loc sim
        kprime2 = STf2(k0i2,prodShocki,tt);
        CS_sim2(iisim,tt) = CTf2(k0i2,prodShocki,tt);% record consmptn sim
        k0i2 = IGTf2(k0i2,prodShocki,tt);         % get new loc for savings
        Assets_sim2(iisim,tt+1) = kprime2;       % record asset sim
        Assets_sim2i(iisim,tt+1) = k0i2;         % record asset loc sim
        
    end
    
end
disp('done finite simulation')
%% 2.5) Plot of simulation (finit agnts) consumption with given parameters:
% Plot of mean of simulation:
figure;
plot(1:T,nanmean( shocks_sim1  ) ); 
title('2.2) Simulation starting with zero wealth and worst income shock');
xlabel('time'); 

%figure;
hold on;
plot(1:T,nanmean( Assets_sim1(:,2:end)  ) ); 
plot(1:T,nanmean( CS_sim1  ) ); 

plot(1:T,nanmean( Assets_sim2(:,2:end)  ) ); 
plot(1:T,nanmean( CS_sim2  ) ); 

hold off;
legend('Income shock','Savings \delta = 0.2','Consumption \delta = 0.2',...
    'Savings \delta = 0.4','Consumption \delta = 0.4','Location','East');
%saveas(gcf,'fig01_5_simFiniteAgents.eps','epsc2');



%% 2.6) Trying to get a hump shape maybe, relax the constraint of assets
%now compute the environment with sigmaepsilon = 0.4

%invert rho and r...
r               = 0.02;         % interest rate
rho             = 0.04;         % bbeta = 1 / (1+rho)

quietdown = 0;
knewBound = -0.4;
disp('Running finitely lived agents VFI with relaxed Constraint')
kHump = [linspace(knewBound,min(k)-0.00001,100)'; k];
tic;
[VTf,IGTf,STf,CTf] = vfi_01_finite(kHump,a,ap9,bbeta,r,u,d,tol,maxiter,Vinit,T);
toc;

%now compute the environment with sigmaepsilon = 0.4
[asigma4,apsigma4]    = rouwenhorst_ram(N,p,p,0.4*sqrt(N-1));
a2sigma4       = exp(asigma4);
tic;
[VTf2,IGTf2,STf2,CTf2] = vfi_01_finite(kHump,a2sigma4,apsigma4,bbeta,r,u,d,...
    tol,maxiter,Vinit,T);
toc;

% plots for finite stuff
figure;
tplot = 55;
plot(kHump,VTf(:,[1 5 9],tplot));
title('2.3) Value function for dying agent in t = 55');
legend('Min shock','Mean shock','Max shock','Location','NorthWest');
xlabel('Assets'); 
ylabel('Value');

%saveas(gcf,'fig03_1_val.eps','epsc2');
%
figure;
plot(kHump,STf(:,[1 5 9],tplot)); 
title('2.3) Savings Policy function for dying agent, in t = 55');
xlabel('Assets'); 
ylabel('Policy (savings)');
legend('Min shock','Mean shock','Max shock','Location','NorthWest');
%saveas(gcf,'fig03_2_assets.eps','epsc2');

% Compare young/old
figure;
plot(kHump,CTf(:,[1 5 9],1)); 
title('2.3) Consmptn Policy function for dying agent, in t = 1 (young)');
xlabel('Assets'); 
ylabel('Policy (consumption)');
legend('Min shock','Mean shock','Max shock','Location','NorthWest');

figure;
plot(kHump,CTf(:,[1 5 9],tplot)); 
title('2.3) Consumption Policy function for dying agent, in t = 55');
xlabel('Assets'); 
ylabel('Policy (consumption)');
legend('Min shock','Mean shock','Max shock','Location','NorthWest');

%saveas(gcf,'fig03_3_cons.eps','epsc2');

% Simulation for a finitely lived agent
%   Similar to previous, but with Value Function changing...
disp('Starting finite simulation')
numsim = 1000;
mc = dtmc(mTransition);
shocks_sim1 = NaN(numsim,T);    % contains all the shock for  simulations
shocks_sim1i = NaN(numsim,T);   % contains all the loc shocks 4 simltns
Assets_sim1 = NaN(numsim,T+1);  % contains all the assets simulations
Assets_sim1i = NaN(numsim,T+1); % contains all the loc assets simltns
CS_sim1 = NaN(numsim,T);        % contains all consumptions simulation
Assets_sim2 = NaN(numsim,T+1);  % contains all the assets simulations
Assets_sim2i = NaN(numsim,T+1); % contains all the loc assets simltns
CS_sim2 = NaN(numsim,T);        % contains all consumptions simulation

for iisim = 1:numsim
    
    [~, k0i] = min(abs(kHump)); % start simulation with no wealth
    k0i2 = k0i;
    x0 = zeros(1,mc.NumStates);  
    x0(1) = 1;                       % start wimulations with lowest shock
    shocksT = simulate(mc,T-1,'X0',x0);  % productivity shocks
    shocks_sim1(iisim,:) = exp(vProductivity(shocksT)'); % recrd shocks sim
    shocks_sim1i(iisim,:) = shocksT';    % record shocks loc from simln
    Assets_sim1(iisim,1) = kHump(k0i);       % record asset holding...
    Assets_sim1i(iisim,1) = k0i;         % record asset loc sim
    Assets_sim2(iisim,1) = kHump(k0i2);       % record asset holding...
    Assets_sim2i(iisim,1) = k0i2;         % record asset loc sim
    
    for tt = 1:T
        prodShocki = shocksT(tt);
        prodShock = exp(vProductivity(prodShocki));
        kprime = STf(k0i,prodShocki,tt);
        CS_sim1(iisim,tt) = CTf(k0i,prodShocki,tt);% record consumption sim
        k0i = IGTf(k0i,prodShocki,tt);           % get new loc for savings
        Assets_sim1(iisim,tt+1) = kprime;       % record asset sim
        Assets_sim1i(iisim,tt+1) = k0i;         % record asset loc sim
        kprime2 = STf2(k0i2,prodShocki,tt);
        CS_sim2(iisim,tt) = CTf2(k0i2,prodShocki,tt);% record consmptn sim
        k0i2 = IGTf2(k0i2,prodShocki,tt);         % get new loc for savings
        Assets_sim2(iisim,tt+1) = kprime2;       % record asset sim
        Assets_sim2i(iisim,tt+1) = k0i2;         % record asset loc sim
        
    end
    
end
disp('done finite simulation')

% Plot of simulation (finit agnts) consumption with given parameters:
% Plot of mean of simulation:
figure;
plot(1:T,nanmean( shocks_sim1  ) ); 
title('2.2) Simulation starting with zero wealth and worst income shock');
xlabel('time'); 

%figure;
hold on;
plot(1:T,nanmean( Assets_sim1(:,2:end)  ) ); 
plot(1:T,nanmean( CS_sim1  ) ); 

plot(1:T,nanmean( Assets_sim2(:,2:end)  ) ); 
plot(1:T,nanmean( CS_sim2  ) ); 

hold off;
legend('Income shock','Savings \delta = 0.2','Consumption \delta = 0.2',...
    'Savings \delta = 0.4','Consumption \delta = 0.4','Location','East');
%saveas(gcf,'fig01_6_simFiniteAgentsDying.eps','epsc2');


%% 2.7) Death probability and specific income structure...  X-)
theta = 0.7;
xtext = textread('incprofile.txt','%f');
%xtext = textscan('incprofile.txt');
ybar = xtext(:,1);
ybar = [ybar; theta*ones(16,1)*ybar(45)];

knewBound = -0.4;
kHump = [linspace(knewBound,min(k)-0.00001,100)'; kconstrained];
k = kHump;
%probas of living next period:
phis = textread('survs.txt','%f');
T = 60;
disp('Running finitely lived agents VFI, with proba of dying')
tic;
[VTd,IGTd,STd,CTd] = vfi_01_finite_dying(k,a,ap9,bbeta,phis,r,...
    u,d,tol,maxiter,Vinit,T,ybar);
toc;
%now compute the environment with sigmaepsilon = 0.4
[asigma4,apsigma4]    = rouwenhorst_ram(N,p,p,0.4*sqrt(N-1));
a2sigma4       = exp(asigma4);
tic;
[VTd2,IGTd2,STd2,CTd2] = vfi_01_finite_dying(k,a2sigma4,apsigma4,bbeta,...
    phis,r,u,d,tol,maxiter,Vinit,T,ybar);
toc;

%% plots for finite stuff, shocks given by data:
figure;
tplot = 55;
plot(k,VTd(:,[1 5 9],tplot));
title('2.7) Value function for dying agent in t = 55');
legend('Min shock','Mean shock','Max shock','Location','NorthWest');
xlabel('Assets'); 
ylabel('Value');

%saveas(gcf,'fig01_1_val.eps','epsc2');
%
figure;
plot(k,STd(:,[1 5 9],tplot)); 
title('2.7) Savings Policy function for dying agent, in t = 55');
xlabel('Assets'); 
ylabel('Policy (savings)');
legend('Min shock','Mean shock','Max shock','Location','NorthWest');
%saveas(gcf,'fig01_2_assets.eps','epsc2');

%% Compare young/old
figure;
plot(k,CTd(:,[1 5 9],1)); 
title('2.7) Consmptn Policy function for dying agent, in t = 1 (young)');
xlabel('Assets'); 
ylabel('Policy (consumption)');
legend('Min shock','Mean shock','Max shock','Location','NorthWest');

figure;
plot(k,CTd(:,[1 5 9],tplot)); 
title('2.7) Consumption Policy function for dying agent, in t = 55');
xlabel('Assets'); 
ylabel('Policy (consumption)');
legend('Min shock','Mean shock','Max shock','Location','NorthWest');

%saveas(gcf,'fig01_3_cons.eps','epsc2');

%% Simulation for a finitely lived agent with proba of dying
%   Similar to previous, but with Value Function changing...
%           and NaN because we might kill him...
numsim = 1000;
T = 59;
mc = dtmc(mTransition);
shocks_sim1 = NaN(numsim,T+1);    % contains all the shock for  simulations
shocks_sim1i = NaN(numsim,T+1);   % contains all the loc shocks 4 simltns
Assets_sim1 = NaN(numsim,T+2);  % contains all the assets simulations
Assets_sim1i = NaN(numsim,T+2); % contains all the loc assets simltns
CS_sim1 = NaN(numsim,T+1);        % contains all consumptions simulation
Assets_sim2 = NaN(numsim,T+2);  % contains all the assets simulations
Assets_sim2i = NaN(numsim,T+2); % contains all the loc assets simltns
CS_sim2 = NaN(numsim,T+1);        % contains all consumptions simulation

for iisim = 1:numsim
    
    [~, k0i] = min(abs(k)); % start simulation with no wealth
    k0i2 = k0i;
    x0 = zeros(1,mc.NumStates);      % for simulation, get zeros size na
    x0(1) = 1;                       % start wimulations with lowest shock
    shocksT = simulate(mc,T,'X0',x0);  % productivity shocks
    shocks_sim1i(iisim,:) = shocksT';    % record shocks loc from simln
    Assets_sim1(iisim,1) = k(k0i);       % record 1st asset holding...
    Assets_sim1i(iisim,1) = k0i;         % record 1st asset loc sim
    Assets_sim2(iisim,1) = k(k0i2);       % record 1st asset holding...
    Assets_sim2i(iisim,1) = k0i2;         % record 1st asset loc sim
    
    for tt = 1:T+1
        prodShocki = shocksT(tt);
        prodShock = exp(vProductivity(prodShocki))*ybar(tt);
        shocks_sim1(iisim,tt) = prodShock;      % record shocks sim
        kprime = STd(k0i,prodShocki,tt);
        consumpt = CTd(k0i,prodShocki,tt);
        CS_sim1(iisim,tt) = consumpt;           % record consumption sim
        k0i = IGTd(k0i,prodShocki,tt);           % get new loc for savings
        Assets_sim1(iisim,tt+1) = kprime;       % record asset sim
        Assets_sim1i(iisim,tt+1) = k0i;         % record asset loc sim
        %once the consumption and savings are set, kill person maybe:
        if rand() > phis(tt)  %kill the person with certain proba...
            break
        end
    end
    
    for tt = 1:T+1
        prodShocki = shocksT(tt);
        prodShock = exp(vProductivity(prodShocki))*ybar(tt);
        shocks_sim1(iisim,tt) = prodShock;      % record shocks sim
        kprime2 = STd2(k0i2,prodShocki,tt);
        consumpt2 = CTd2(k0i2,prodShocki,tt);
        CS_sim2(iisim,tt) = consumpt2;           % record consumption sim
        k0i2 = IGTd2(k0i,prodShocki,tt);           % get new loc for savings
        Assets_sim2(iisim,tt+1) = kprime2;       % record asset sim
        Assets_sim2i(iisim,tt+1) = k0i2;         % record asset loc sim
        %once the consumption and savings are set, kill person maybe:
        if rand() > phis(tt)  %kill the person with certain proba...
            break
        end
    end
    
end
disp('Computed Simulation of dying agent with shocks given by data. ok')

% 2.7end) Plot of simulation (finit agnts) consumption with shocks by data:
% Plot of mean of simulation:
figure;
plot(1:T+1,nanmean( shocks_sim1  ) ); 
title('2.2) Simulation starting with zero wealth and worst income shock');
xlabel('time'); 

%figure;
hold on;
plot(1:T+1,nanmean( Assets_sim1(:,2:end)  ) ); 
plot(1:T+1,nanmean( CS_sim1  ) ); 
plot(1:T+1,nanmean( Assets_sim2(:,2:end)  ) ); 
plot(1:T+1,nanmean( CS_sim2  ) ); 

hold off;
legend('Income shock','Savings \delta = 0.2','Consumption \delta = 0.2',...
    'Savings \delta = 0.4','Consumption \delta = 0.4','Location','Northwest');
%saveas(gcf,'fig01_7_simFiniteAgents_incomeData.eps','epsc2');


%% 2.8)
%compare the consumption generated by fdz-Vllvrd&Krueger(2007):
consTarget = dlmread('consprofile.txt');
consTarget_yearly = ones(size(consTarget,1)/4,1);
for ii = 1:(size(consTarget,1)/4)
    consTarget_yearly(ii) = prod(consTarget( ((ii-1)*4+1):ii*4,2 ) );
end

figure
plot(1:T+1,nanmean( CS_sim2  ) ); 
title('2.8) Average Consumption from simulation vs. Target');
xlabel('time'); 
hold on
plot(3:68,consTarget_yearly)
ylabel('Policy (consumption)');
hold off;
legend('Sim Cons','Cons JFV&DK','Location','Northwest');
%saveas(gcf,'fig01_3_cons.eps','epsc2');

%% 2.9) Blundell,Pistaferri and Preston,
%   Consumption insurance coefficient
%%%% COMPUTE THIS NUMBER FOR THE INFINITELY LIVED AGENT 
% when delta = 0 and 0.99
k = kconstrained;
quietdown = 1; %makes the function not print iterations...
[a9,ap9]    = tauchen_ram(N,0,sigmaepsilon,3);
a = exp(a9);
[V9a,IG9a,S9a,C9a] = vfi_01_infty(k,a,ap9,bbeta,r,u,d,tol,maxiter,Vinit,quietdown);
disp('Starting infinite simulation for delta = 0')
numsim = 1000;
mc = dtmc(mTransition);
shocks_sim1 = zeros(numsim,T+1);  % contains all the shock for  simulations
shocks_sim1i = zeros(numsim,T+1); % contains all the loc shocks 4 simltns
Assets_sim1 = zeros(numsim,T+2);  % contains all the assets simulations
Assets_sim1i = zeros(numsim,T+2); % contains all the loc assets simltns
CS_sim1 = zeros(numsim,T);        % contains all consumptions simulation

for iisim = 1:numsim
    
    %k0i = randi([1 nk],1,1);            % take random initial asset level
    %%%% CHECK THIS STEP WHEN CONSTRAINED IS RELAXED TO START AT REAL ZERO!
    [~, k0i] = min(abs(k));     % start simulation with no wealth
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
        kprime = S9a(k0i,prodShocki);
        CS_sim1(iisim,tt) = C9a(k0i,prodShocki);  % record consumption sim
        k0i = IG9a(k0i,prodShocki);               % get new loc for savings
        Assets_sim1(iisim,tt+1) = kprime;       % record asset sim
        Assets_sim1i(iisim,tt+1) = k0i;         % record asset loc sim
        
    end
    
end

disp('done infinite simulation for delta = 0')

covsCY = zeros(numsim,1);
varYY = zeros(numsim,1);
%for tt = 1:T+1
for ii = 1:numsim
    covConsYtt = nancov(CS_sim1(ii,:),shocks_sim1(ii,:));
    covsCY(ii) = covConsYtt(1,2);
    varYY(ii) = covConsYtt(2,2);
end
% covConsY = mean(covsCY)
% varY = mean(varYY)
% vLogY = mean(nanvar( shocks_sim1 ))
consInsCoeff = 1 - covsCY./varYY;
disp('Results for \delta = 0, for every age')
disp(num2str(consInsCoeff'));
m1 = min(consInsCoeff);
m2 = mean(consInsCoeff);
m3 = max(consInsCoeff);
m4 = median(consInsCoeff);
disp(['min(coef) = ', num2str(m1),'  mean(coef) = ', num2str(m2),...
     '  max(coef) = ', num2str(m3),'  median(coef) = ', num2str(m4)]);

%% repeat simulation for delta = 0.99
p = (1+0.99)/2;
Psi = sigmaepsilon*sqrt(N-1);
[a9,ap9]    = rouwenhorst_ram(N,p,p,Psi);
a = exp(a9);
[~,IG9b,S9b,C9b] = vfi_01_infty(k,a,ap9,bbeta,r,u,d,tol,maxiter,Vinit,quietdown);
disp('Starting infinite simulation for delta = 0')
numsim = 1000;
mc = dtmc(mTransition);
shocks_sim1 = zeros(numsim,T+1);  % contains all the shock for  simulations
shocks_sim1i = zeros(numsim,T+1); % contains all the loc shocks 4 simltns
Assets_sim1 = zeros(numsim,T+2);  % contains all the assets simulations
Assets_sim1i = zeros(numsim,T+2); % contains all the loc assets simltns
CS_sim1 = zeros(numsim,T);        % contains all consumptions simulation

for iisim = 1:numsim
    
    %k0i = randi([1 nk],1,1);            % take random initial asset level
    %%%% CHECK THIS STEP WHEN CONSTRAINED IS RELAXED TO START AT REAL ZERO!
    [~, k0i] = min(abs(k));     % start simulation with no wealth
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
        kprime = S9b(k0i,prodShocki);
        CS_sim1(iisim,tt) = C9b(k0i,prodShocki);  % record consumption sim
        k0i = IG9b(k0i,prodShocki);               % get new loc for savings
        Assets_sim1(iisim,tt+1) = kprime;       % record asset sim
        Assets_sim1i(iisim,tt+1) = k0i;         % record asset loc sim
        
    end
    
end

disp('done infinite simulation for delta = 0.99')

covsCY = zeros(numsim,1);
varYY = zeros(numsim,1);
%for tt = 1:T+1
for ii = 1:numsim
    covConsYtt = nancov(CS_sim1(ii,:),shocks_sim1(ii,:));
    covsCY(ii) = covConsYtt(1,2);
    varYY(ii) = covConsYtt(2,2);
end
% covConsY = mean(covsCY)
% varY = mean(varYY)
% vLogY = mean(nanvar( shocks_sim1 ))
consInsCoeff = 1 - covsCY./varYY;
disp('Results for \delta = 0.99, for every age')
disp(num2str(consInsCoeff'));
m1 = min(consInsCoeff);
m2 = mean(consInsCoeff);
m3 = max(consInsCoeff);
m4 = median(consInsCoeff);
disp(['min(coef) = ', num2str(m1),'  mean(coef) = ', num2str(m2),...
     '  max(coef) = ', num2str(m3),'  median(coef) = ', num2str(m4)]);
