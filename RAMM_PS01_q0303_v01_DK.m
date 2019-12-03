%% UPENN, 714, Prof Dirk Krueger, Problem set 01. Question 3
% Rodrigo Morales
% based on a code by Anna Cororaton
% November 2019
% Aiyagari model, general Equilibrium

%Housekeeping
clear; clc; close all;

% parameters
%r       = 0.02;     % r will be given by guess of K.
ddelta          = 0.08;          % depreciation
delta           = 0.8;         % persistence
aalpha          = 0.36;         % capital elasticity (cobb douglas)
rho             = 0.04;         % bbeta = 1 / (1+rho)
bbeta           = 1/(1+rho);    
sigmaepsilon    = 0.4;          % TFP std dev in PSet, called sigma_y
sigma           = 1;            % Elasticity of subs (utility).
csi             = 0.997;         % Convergence parameter for capital
if sigma ==1                    % if sigma ==1, use log utility
    u=@(c) log(c);
else                            % if sigma ~=1, use CRRA
    u=@(c) c.^(1-sigma)./(1-sigma);
end

if delta >= 0.8        % if shock is persistent, do Rouenhorst
    doRouenhorst = 1; 
else                   % if shock is not persistent, do Tauchen
    doRouenhorst = 0;
end

% CALIBRATEION...
T       = 60;       % periods for each simulation
numsim  = 1000;     % number of simulations
nk      = 200;     % number of grid points for K
na      = 11;       % number of grid points for A
N       = na;       % change of name
tol     = 10e-10;   % tolerance level for VFI
maxiter = 1000;     % maximum number of iterations

% Get transition matrix (Tauchen or Rouenhorst)

if doRouenhorst == 1        	% do Rouenhorst
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
else   %doRouenhorst ==01;      % do Tauchen
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
end

% create grid for k and a, now depends on r...

% savings 
savmin    = 0.001;
savmax    = 15;%vProductivity(N)/(1-bbeta); % -natural borrowing constraint
k = curvspace(savmin,savmax,nk,2)'; % use curved grid to enhance accuracy

%% (Infinite) Value Function Iteration...
% Find the equilibrium kapital and interest rate
iteration   = 0;
rguess      = rho - 0.0001;
Kguessmin      = ((rguess + ddelta)/aalpha)^(1/(aalpha - 1));
Kguess = Kguessmin;
r           = aalpha*Kguess^(aalpha-1)-ddelta;
w           = (1-aalpha)*Kguess^aalpha;
a           = w.*exp(a9);  %now w*a, wage, depends on r...
% start in the utility of cash at hand:
Vinit = u( (1+r)*repmat(k,1,na)+repmat(a,nk,1) );
lambdainit = (1/(na*nk))*ones(nk,na);  %initial stationary distribution
d = 100;
maxiterLoc = 100;
while d > 1e-4 && iteration < maxiterLoc
    iteration = iteration + 1;
    a       = w.*exp(a9);  %now w*a, wage, depends on r...
    %Kguess
    qtdn =1;
    [V,IG,S,C,L_policy] = vfi_01_infty_labor(k,kappa,a,ap9,bbeta,r,u,100,tol,maxiter,Vinit, qtdn);
    
    %lambdastationary = findLambdaStationary(lambdainit,ap9,IG,k,tol,maxiter);
    lambdastationary = findLambdaStationary_MATRIX(lambdainit,ap9,IG,tol,maxiter);
    
    Knew = sum(sum(lambdastationary.*S));
    d = norm(Kguess-Knew,2);
    % update values
    Vinit = V;

%     if (d<1e-5)
%         csi = 3/4*csi;
%     end
    Kguess2 = csi*Kguess + (1-csi)*Knew;
    if Kguess2 == Kguess
        iteration = maxiterLoc+1;
    else
        Kguess = Kguess2;
    end
    %csi = 99*csi/100;
    r = aalpha*Kguess^(aalpha-1)-ddelta;
    lambdainit = lambdastationary;
    if(mod(iteration-1,20)==0)
        fprintf('KGuess Iteration = %d, Sup Diff = %2.8f\n', iteration , d);
    end
    %pause
end
r
Knew
%
figure
[X, Y] = meshgrid(k,a);
surf(X,Y,lambdastationary')
%
figure
plot(k,lambdastationary(:,1))
hold on
plot(k,lambdastationary(:,na))

%%
% Find E[a(r)] for multiple values of r  (or K)
r           = aalpha*Kguessmin^(aalpha-1)-ddelta;
Vinit       = u( (1+r)*repmat(k,1,na)+repmat(a,nk,1) );
lambdainit = (1/(na*nk))*ones(nk,na);  %initial stationary distribution
iteration = 0;
%Ks = curvspace(Kguessmin,savmax,nk,2); % use curved grid to enhance accuracy
Ks = linspace(Kguessmin,savmax,nk); % use curved grid to enhance accuracy
krs = ones(nk,1);
rrs = ones(nk,1);
for Kguess = Ks
    iteration = iteration + 1;
    
    r = aalpha*Kguess^(aalpha-1)-ddelta;
    
    a       = w.*exp(a9);  %now w*a, wage, depends on r...
    qtdn =1;
    [V,IG,S,C,L_policy] = vfi_01_infty_labor(k,kappa,a,ap9,bbeta,r,u,100,tol,maxiter,Vinit, qtdn);
    
    %lambdastationary = findLambdaStationary(lambdainit,ap9,IG,k,tol,maxiter);
    lambdastationary = findLambdaStationary_MATRIX(lambdainit,ap9,IG,tol,maxiter);
    
    Knew = sum(sum(lambdastationary.*S));
    % record values
    rnew = aalpha*Knew^(aalpha-1)-ddelta;
    rrs(iteration) = r;
    krs(iteration) = Knew;
    %Update
    Vinit = V;
    lambdainit = lambdastationary;
    
    if(mod(iteration-1,20)==0)
        fprintf('KGuess Iteration = %d \n', iteration); 
    end
end

% plot of fpk and E[a(r)]
figure
Ks = linspace(Kguessmin,savmax,nk); % use curved grid to enhance accuracy
plot(Ks,aalpha*Ks.^(aalpha-1)-ddelta, 'blue')
hold on
plot(flip(krs),flip(rrs), 'red')
hold off

%%
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
plot(1:T+1,nanmean( shocks_sim1  ) ); 
title('2.2) Simulation starting with zero wealth and worst income shock');
%title('2.2) Average shocks (income) from simulation');
xlabel('time'); 
%ylabel('Policy (savings)');
%saveas(gcf,'fig01_2_assets.eps','epsc2');

%figure;
hold on;
plot(1:T+1,nanmean( Assets_sim1(:,2:end)  ) ); 
%title('2.2) Average Policy function (savings) from simulation');
xlabel('time'); 
%ylabel('Policy (savings)');
%saveas(gcf,'fig01_2_assets.eps','epsc2');

%figure;
plot(1:T+1,nanmean( CS_sim1  ) ); 
%title('2.2) Average Consumption from simulation');
xlabel('time'); 
%ylabel('Policy (savings)');
hold off;
legend('Income','Savings','Consumption','Location','East');
%saveas(gcf,'fig01_3_cons.eps','epsc2');






