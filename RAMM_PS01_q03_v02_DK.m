%% UPENN, 714, Prof Dirk Krueger, Problem set 01. Question 3
% Rodrigo Morales
% November 2019
% Aiyagari model, general Equilibrium

%Housekeeping
clear; clc; close all;

% parameters
deltas = [0 0.3 0.6 0.9];%[0 0.3 0.6 0.9];
sigmas = [1 3 5]; %[1 3 5];
sigmaepsilons = [0.2 0.4];%[0.2 0.4];
quietdown = 1;


rguesses =zeros(length(deltas),length(sigmas),length(sigmaepsilons));
rguesses(:,:,1) = [4.16 4.14 4.08; 4.14 4.04 3.909; 4.09 3.87 3.58;
    3.93 3.29 2.52];
rguesses(:,:,2) = [4.06 3.78 3.42; 3.96 3.41 2.80; 3.76 2.78 1.81;
    3.31 1.29 -0.35];
rguesses = rguesses -0.1;
RsModel = zeros(length(deltas),length(sigmas),length(sigmaepsilons));
KsModel = zeros(length(deltas),length(sigmas),length(sigmaepsilons));

ii = 0;
jj = 0;
kk = 0;
countersols = 0;
tic

for delta = deltas
    %fprintf('Solving delta = %2.8f\n', delta);
    ii = ii + 1;
    jj = 0;
    for sigma = sigmas
        %fprintf('Solving sigma = %2.8f\n', sigma);
        jj = jj+1;
        kk = 0;
        for sigmaepsilon = sigmaepsilons
            countersols = countersols + 1;
            fprintf('Iteration of solution: = %d of 24 \n', countersols);
            %fprintf('Solving sigmaepsilon = %2.8f\n', sigmaepsilon);
            kk = kk+1;
            [rmodeli Knewi ki ai nai lambdastationaryi] = function_q03(sigma, sigmaepsilon, delta, quietdown, rguesses(ii,jj,kk) );
            RsModel(ii,jj,kk) = rmodeli;
            KsModel(ii,jj,kk) = Knewi;
        end
    end
end
toc;

RsModel

k = ki;
a = ai;
na = nai;
lambdastationary = lambdastationaryi;

%
figure
[X, Y] = meshgrid(k,a);
surf(X,Y,lambdastationary')
title('Stationary Distribution of kapital and shocks)')
%
figure
plot(k,lambdastationary(:,1))
hold on
plot(k,lambdastationary(:,na))
hold off;
title('Stationary distribution of kapital for worst and highest shock')










%%
% Find E[a(r)] for multiple values of r  (or K)
ddelta          = 0.08;          % depreciation
%delta           = 0.8;         % persistence
aalpha          = 0.36;         % capital elasticity (cobb douglas)
rho             = 0.04;         % bbeta = 1 / (1+rho)
bbeta           = 1/(1+rho);    
%sigmaepsilon    = 0.4;          % TFP std dev in PSet, called sigma_y
%sigma           = 1;            % Elasticity of subs (utility).
csi             = 0.75;   

Kguessmin      = ((rguess + ddelta)/aalpha)^(1/(aalpha - 1));
r           = aalpha*Kguessmin^(aalpha-1)-ddelta;
Vinit       = u( (1+r)*repmat(k,1,na)+repmat(a,nk,1) );
lambdainit = (1/(na*nk))*ones(nk,na);  %initial stationary distribution
iteration = 0;
%Ks = curvspace(Kguessmin,savmax,nk,2); %use curved grid 2 enhance accuracy
Ks = linspace(Kguessmin,savmax,nk); % use curved grid to enhance accuracy
krs = ones(nk,1);
rrs = ones(nk,1);
for Kguess = Ks
    iteration = iteration + 1;
    
    r = aalpha*Kguess^(aalpha-1)-ddelta;
    w           = (1-aalpha)*Kguess^aalpha;
    a       = w.*exp(a9);  %now w*a, wage, depends on r...
    qtdn =1;
    [V,IG,S,C] = vfi_01_infty(k,a,ap9,bbeta,r,u,100,tol,maxiter,Vinit, qtdn);
    
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


