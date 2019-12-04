%% UPENN, 714, Prof Dirk Krueger, Problem set 01. Question 3
% Rodrigo Morales
% November 2019
% Aiyagari model, general Equilibrium

function [r Knew k a na lambdastationary] = function_q03(sigma, ...
    sigmaepsilon, delta, quietdown, rguess )

    % parameters
    %r       = 0.02;     % r will be given by guess of K.
    ddelta          = 0.08;          % depreciation
    %delta           = 0.8;         % persistence
    aalpha          = 0.36;         % capital elasticity (cobb douglas)
    rho             = 0.04;         % bbeta = 1 / (1+rho)
    bbeta           = 1/(1+rho);    
    %sigmaepsilon    = 0.4;          % TFP std dev in PSet, called sigma_y
    %sigma           = 1;            % Elasticity of subs (utility).
    csi             = 0.75;         % Convergence parameter for capital
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
        [a9,ap9]    = rouwenhorst_ram(N,p,p,Psi);
        [vProductivity,mTransition] = rouwenhorst_ram(N,p,p,Psi);
        % calculate long run distribution
        a9star      = pstar(ap9);
        prodStatnry = pstar(mTransition);
    else   %doRouenhorst ==01;      % do Tauchen
        [a9,ap9]    = tauchen_ram(N,delta,sigmaepsilon,3);
        [vProductivity,mTransition]  = tauchen_ram(N,delta,sigmaepsilon,3);
        vProductivity = vProductivity';
        % previous is equivalent to: (slight variation becauseof sensibility to
        % calculation. Bottom line...   works :) ,   but is senstitive :(
        %[Z,Zprob] = tauchen(N,0,delta,(1-delta^2)^(1/2)*sigmaepsilon,3);
        % calculate long run distribution
        a9star      = pstar(ap9);
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
     if ~exist('rguess','var')
         % third parameter does not exist, so default it to something
          rguess      = rho - 0.0001;
     end
    Kguessmin      = ((rguess + ddelta)/aalpha)^(1/(aalpha - 1));
    Kguess = Kguessmin;
    r           = aalpha*Kguess^(aalpha-1)-ddelta;
    w           = (1-aalpha)*Kguess^aalpha;
    a           = w.*exp(a9);  %now w*a, wage, depends on r...
    % start in the utility of cash at hand:
    Vinit = u( (1+r)*repmat(k,1,na)+repmat(a,nk,1) );
    lambdainit = (1/(na*nk))*ones(nk,na);  %initial stationary distribution
    d = 100;
    maxiterLoc = 500;
    while d > 1e-4 && iteration < maxiterLoc
        iteration = iteration + 1;
        w           = (1-aalpha)*Kguess^aalpha;
        a       = w.*exp(a9);  %now w*a, wage, depends on r...
        %Kguess
        qtdn =1;
        [V,IG,S,C] = vfi_01_infty_opt(k,a,ap9,bbeta,r,u,100,tol,maxiter,Vinit,qtdn);

        %lambdastationary = findLambdaStationary(lambdainit,ap9,IG,k,tol,maxiter);
        lambdastationary = findLambdaStationary_MATRIX(lambdainit,ap9,IG,tol,maxiter);

        Knew = sum(sum(lambdastationary.*S));
        %check that the new K won't send r to somewhere crazy, above rho
        if Knew < Kguessmin
            Knew = Kguessmin + 0.01;
            Kguessmin =Kguessmin + 0.000001;
        end
        d = norm(Kguess-Knew,2);
        % update values
        Vinit = V;

    %     if (d<1e-5)
    %         csi = 3/4*csi;
    %     end
        Kguess2 = csi*Kguess + (1-csi)*Knew;
    %     if Kguess2 == Kguess
    %         iteration = maxiterLoc+1;
    %     else
    %         Kguess = Kguess2;
    %     end
        KguessOld   = Kguess;
        Kguess      = Kguess2;

        %csi = 99*csi/100;
        r = aalpha*Kguess^(aalpha-1)-ddelta;
        lambdainit = lambdastationary;

        if quietdown ~=1
            fprintf('KGuess Iteration = %d, Sup Diff = %2.8f\n', iteration , d);
            fprintf('KGuess0 = %2.8f, E[a(r)] = %2.8f, Kbisect = %2.8f\n', KguessOld , Knew, Kguess);
        end
    %     if(mod(iteration,20)==0 || iteration == 1)
    %         fprintf('KGuess Iteration = %d, Sup Diff = %2.8f\n', iteration , d);
    %         fprintf('KGuess0 = %2.8f, E[a(r)] = %2.8f, Kbisect = %2.8f\n', KguessOld , Knew, Kguess);
    %     end
        %pause
    end

    if quietdown ~= 1
        fprintf('r_solved = %2.8f, Knew = %2.8f\n', r , Knew);
    end
    fprintf('Finished solving for sigma= %2.8f, sigmaepsilon= %2.8f, delta = %2.8f\n', sigma, sigmaepsilon , delta);
end


