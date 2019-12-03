%% FNCE 937 - Problem Set 1
% Anna Cororaton
% September 2014

clear all; clc; close all;



%% parameters
T       = 5000;     % periods for each simulation (1000 or 5000)
r       = 0.05;     % interest rate
delta   = 0.1;      % depreciation
theta1  = 0.3;      % capital share
theta2  = 0.6;      % labor share
W       = 2;        % wage
sigma   = 0.1;      % TFP std dev
rho     = 0.8;      % TFP persistence (0.8 or 0.95 or 0.99)
b0      = 0;        % adjustment cost parameter
b1p     = 0.5;      % adjustment cost parameter when i>delta*k
b1m     = 0.5;      % adjustment cost parameter when i<=delta*k
nk      = 200;      % number of grid points for K
na      = 9;        % number of grid points for A
tol     = 10e-10;   % tolerance level for VFI
maxiter = 1000;     % maximum number of iterations
d       = 100;      % distance metric



%% part 2.2a
% set long run mean of shocks

gamma   = theta1/(1-theta2);
abar    = ((r+delta)/(gamma*(1-theta2)))^(1-theta2)*W^theta2*theta2^(-theta2)*exp(-sigma^2/(2*(1-theta2)));



%% part 2.2b 

% use tauchen hussey method
[a5,ap5]    = tauchenhussey(5,log(abar),rho,sigma,sigma);
[a9,ap9]    = tauchenhussey(na,log(abar),rho,sigma,sigma);
[a15,ap15]  = tauchenhussey(15,log(abar),rho,sigma,sigma);

% calculate long run distribution
a5star      = pstar(ap5);
a9star      = pstar(ap9);
a15star     = pstar(ap15);



%% part 2.2c, 2.2d

% simulation
[muhat5 sigmahat5 rhohat5]      = calc_markov(T,1,a5,ap5); 
[muhat9 sigmahat9 rhohat9]      = calc_markov(T,1,a9,ap9); 
[muhat15 sigmahat15 rhohat15]   = calc_markov(T,1,a15,ap15);

% theoretical moments
mu5_th  = a5star*a5;
mu9_th  = a9star*a9;
mu15_th = a15star*a15;

sd5_th  = sqrt(sum(a5star*(a5-mu5_th).^2));
sd9_th  = sqrt(sum(a9star*(a9-mu9_th).^2));
sd15_th = sqrt(sum(a15star*(a15-mu15_th).^2));




%% create grid for k and a

% profit calculation written as pi = Ak^gamma
a       = exp(a9);
amax    = max(a);
amin    = min(a);
x       = W^(theta2/(theta2-1))*(theta2^(theta2/(1-theta2))-theta2^(1/(1-theta2)));
Amax    = x*amax^(1/(1-theta2));
Amin    = x*amin^(1/(1-theta2));
Abar    = x*abar^(1/(1-theta2));

% calculate A for all a shocks which is useful later
A       = x*a.^(1/(1-theta2)); 

% capital 
kmin    = (gamma*Amin/(r+delta))^(1/(1-gamma));
kmax    = 4; %((r+delta)/(gamma*Amax))^(1/(gamma-1));
k       = linspace(kmin,kmax,nk)';



%% part 2.3a

pi_Amax = Amax*k.^(gamma);
pi_Amin = Amin*k.^(gamma);
pi_Abar = Abar*k.^(gamma);

figure;
plot(k,[pi_Amax,pi_Abar,pi_Amin]); 
axis tight; 
title('Profit function'); 
legend('Max a','Mean a','Min a','Location','NorthWest');
xlabel('Capital'); 
ylabel('Profit');
saveas(gcf,'fig23a.eps','epsc2');



%% part 2.3b
tic;
[V,G,inv1,phi,fin] = vfi_hw1(k,A,ap9,r,theta1,theta2,delta,b0,b1p,b1m,d,tol,maxiter);
toc;

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




%% part 2.3d - 2
b0      = 0;
b1p     = 10;
b1m     = 10;
[~,~,inv2] = vfi_hw1(k,A,ap9,r,theta1,theta2,delta,b0,b1p,b1m,d,tol,maxiter);

figure;
plot(k,inv2(:,[1 5 9])); 
title('Optimal Investment, b0 = 0, b1=10');
legend('Max a','Mean a','Min a','Location','NorthWest');
xlabel('Capital'); 
ylabel('Investment');

saveas(gcf,'fig23d2.eps','epsc2');



%% part 2.3d - 3
b0      = 0;
b1p     = 0.5;
b1m     = 10;
[~,~,inv3] = vfi_hw1(k,A,ap9,r,theta1,theta2,delta,b0,b1p,b1m,d,tol,maxiter);

figure;
plot(k,inv3(:,[1 5 9])); 
title('Optimal Investment, b0 = 0, b1+=0.5, b1-=10');
legend('Max a','Mean a','Min a','Location','NorthWest');
xlabel('Capital'); 
ylabel('Investment');

saveas(gcf,'fig23d3.eps','epsc2');



%% part 2.3d - 4
b0      = 0.02;
b1p     = 0.05;
b1m     = 0.05;
[~,~,inv4] = vfi_hw1(k,A,ap9,r,theta1,theta2,delta,b0,b1p,b1m,d,tol,maxiter);

figure;
plot(k,inv4(:,[1 5 9])); 
title('Optimal Investment, b0 = 0.02, b1=0.05');
legend('Max a','Mean a','Min a','Location','SouthWest');
xlabel('Capital'); 
ylabel('Investment');

saveas(gcf,'fig23d4.eps','epsc2');



%% part 2.3d
inv5 = [inv1(1:(nk/2),5) inv2(1:(nk/2),5) inv3(1:(nk/2),5) inv4(1:(nk/2),5)];
plot(k(1:(nk/2)),inv5);
legend('b0 = 0, b1 = 0.5','b0 = 0, b1=10','b0 = 0, b1+=0.5, b1-=10','b0 = 0.02, b1=0.05','Location','SouthWest');
title('Optimal Investment');

saveas(gcf,'fig23d.eps','epsc2');



