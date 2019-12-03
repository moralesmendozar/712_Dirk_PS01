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
    u=@(c) c^(1-sigma)/(1-sigma);
end

if delta < 0.8          % if shock is not persistent, do Tauchen
    doRouenhorst = 0; 
else                    % if shock is persistent, do Rouenhorst
    doRouenhorst =1;
end

% CALIBRATE: !!!
T       = 1000;     % periods for each simulation (1000 or 5000)
numsim  = 1000;
nk      = 2000;      % number of grid points for K
na      = 11;       % number of grid points for A
N       = na;       % change of name
tol     = 10e-10;   % tolerance level for VFI
maxiter = 1000;     % maximum number of iterations
d       = 100;      % distance metric

%% Get transition matrix (Tauchen or Rouenhorst)
doRouenhorst = 0;
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
    p = (1-delta)/2;
    Psi = sigmaepsilon*sqrt(N-1);
    [a5,ap5]    = rouwenhorst_ram(5,p,p,Psi);
    [a9,ap9]    = rouwenhorst_ram(N,p,p,Psi);
    [a15,ap15]  = rouwenhorst_ram(15,p,p,Psi);
    p = (1+delta)/2;
    Psi = sigmaepsilon*sqrt(N-1);
    [vProductivity,mTransition] = rouwenhorst_ram(N,p,p,Psi);
    vProductivity = vProductivity';
    % calculate long run distribution
    a5star      = pstar(ap5);
    a9star      = pstar(ap9);
    a15star     = pstar(ap15);
    prodStatnry = pstar(mTransition);
end

% following commands WILL BE USED IN SIMULATION:
%mc = dtmc(mTransition);
%X = simulate(mc,numsteps);

%% Simulation stuff, can use calc_markov(),  or simulate()

% simulation
[muhat5, sigmahat5, rhohat5]      = calc_markov(T,1,a5,ap5); 
[muhat9, sigmahat9, rhohat9]      = calc_markov(T,1,a9,ap9); 
[muhat15, sigmahat15, rhohat15]   = calc_markov(T,1,a15,ap15);

% theoretical moments
mu5_th  = a5star*a5;
mu9_th  = a9star*a9;
mu15_th = a15star*a15;

sd5_th  = sqrt(sum(a5star*(a5-mu5_th).^2));
sd9_th  = sqrt(sum(a9star*(a9-mu9_th).^2));
sd15_th = sqrt(sum(a15star*(a15-mu15_th).^2));

%% create grid for k and a

a       = exp(a9);
% savings 
savmin    = 0.001;
savmax    = vProductivity(N)/(1-bbeta); % -natural borrowing constraint
k       = linspace(savmin,savmax,nk)';
Vinit       = repmat(k,1,na);     % value function init...

%% Value Function Iteration...
tic;
[V,IG,S,C] = vfi_01_infty(k,a,ap9,bbeta,r,u,d,tol,maxiter,Vinit);
toc;
%%
% plots
figure;
plot(k,V(:,[1 5 9])); 
title('Value function');
legend('Max a','Mean a','Min a','Location','NorthWest');
xlabel('Capital'); 
ylabel('Value');

saveas(gcf,'fig23b_val.eps','epsc2');

figure;
plot(k,G(:,[1 5 9])); 
title('Policy function');
xlabel('Capital'); 
ylabel('Policy');
legend('Max a','Mean a','Min a','Location','NorthWest');
%%
saveas(gcf,'fig23b_pol.eps','epsc2');




%% part 2.3c
figure;
plot(k,inv1(:,[1 5 9])); 
title('Optimal Investment, b0 = 0, b1 = 0.5');
legend('Max a','Mean a','Min a','Location','SouthWest');
xlabel('Capital'); 
ylabel('Investment');

saveas(gcf,'fig23c_inv.eps','epsc2');

ind     = (fin<0);
figure;
plot(k,-fin(:,[1 5 9]).*ind(:,[1 5 9])); 
title('Financing');
xlabel('Capital'); 
ylabel('Financing');
legend('Max a','Mean a','Min a','Location','NorthEast');

saveas(gcf,'fig23c2_fin.eps','epsc2');



%% part 2.3d
inv5 = [inv1(1:(nk/2),5) inv2(1:(nk/2),5) inv3(1:(nk/2),5) inv4(1:(nk/2),5)];
plot(k(1:(nk/2)),inv5);
legend('b0 = 0, b1 = 0.5','b0 = 0, b1=10','b0 = 0, b1+=0.5, b1-=10','b0 = 0.02, b1=0.05','Location','SouthWest');
title('Optimal Investment');

saveas(gcf,'fig23d.eps','epsc2');



