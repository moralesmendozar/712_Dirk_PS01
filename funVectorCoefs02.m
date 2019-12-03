% RAMM & Shasha Wang PS2 Q2 FNCE 937

% If one wants to do MULTIGRID to speed up the process, just change
% kGridLength to a vector

function [coefsRhoLambda Vrenewed] = funVectorCoefs02(rrho,llambda,...
    Vinit,otherCoefs,kGridLength,saveDoc,numsteps)
%this function finds the coefficients of a linear regression
%   for a simulation of shocks, given the solution of a Bellman 
%   function of a firm with an exogenous rule for default...
% llambda proportional cost of issuing equity
% Parametization
M = otherCoefs(1);
Rf = 1/M;
r = 1/M-1; % interest rate for notation convenience
ddelta = otherCoefs(2);
aalphaK = otherCoefs(3);
aalphaL = otherCoefs(4);
W = otherCoefs(5); % wage
ttaoC=otherCoefs(6); % corporate tax rate
% k_grid & b_grid & a_grid & m_a_prob, a's transition probability matrix
Nk = otherCoefs(7);
kMin = otherCoefs(8);
kSteadyState = otherCoefs(9);
tempK = kSteadyState^((aalphaK+aalphaL-1)/(1-aalphaL));
% set a such that k_SteadyState equals how much you set it to be;
aMean = ((r+ddelta)/aalphaK/tempK)^(1-aalphaL)*(W/aalphaL)^aalphaL; 

% a_grid and a's transition matrix
m = otherCoefs(10); % parameter for tauchen method
Na = otherCoefs(11);
ssigma = otherCoefs(12);
[grid_a_log,m_a_prob] = tauchen(Na,log(aMean),rrho,ssigma,m);
grid_a = exp(grid_a_log)';
grid_a_minus = grid_a;
aMax = max(grid_a);

% k_grid
kMax = (aMax/ddelta)^(1/(1-aalphaK));
kMax = min(2*kSteadyState, kMax); % Tighten the grid
grid_k = curvspace(kMin,kMax,Nk,2)'; % I use curved grid to enhance accuracy

% b_grid grid for bond
% b_grid should be finer to see the difference in default probability under different productivity shocks at steady state
Nb = otherCoefs(13);
grid_b = curvspace(0,kMax,Nb,2)'; % To cover up as wide leverage level as possible

% invariant distribution of a_grid - vDistribution_a0
vDistribution_a0=( 1/Na )*ones(Na,1); % initial guess
vDistribution_a=vDistribution_a0;
distance=100; 
tolerance= otherCoefs(14);
iteration = 0;
while distance>tolerance
    distribution = vDistribution_a'*m_a_prob;
    distance=sum(abs(distribution-vDistribution_a'));
    vDistribution_a = distribution';
    iteration = iteration+1;
end

clear vDistribution_a0 distribution;

% Functions
% create functions for convenience
profitFunction = @(a,k)a.* k.^aalphaK.* ((k.^aalphaK *...
    aalphaL).* a/W).^(aalphaL/(1-aalphaL)) - W * ((k.^aalphaK * ...
    aalphaL).* a/W).^(1/(1-aalphaL));
nonDefaultFunction = @ (profit,k,bond)((1-ttaoC)*profit + ...
    (1-ddelta)*k > bond);
investmentFunction = @(k,kPrime)kPrime - (1-ddelta)*k;
%k_prime usually is k_grid
taxPaymentsFunction = @(k,bond,profit,RbMinus)ttaoC *...
    (profit - ddelta*k - bond.* (RbMinus-1).* ((1-ttaoC)*profit + ...
    (1-ddelta)*k > bond)); % note the non-default indicator
dividentFunction = @(profit,investment,bond,bondPrime,RbMinus,...
    taxPayments,mIsDefaultNextPeriod)(profit - investment  +...
    bondPrime.*(1-mIsDefaultNextPeriod) - RbMinus.* bond - ...
    taxPayments).*(1 + llambda * ((profit - investment + ...
    bondPrime.*(1-mIsDefaultNextPeriod) - RbMinus.* bond - ...
    taxPayments) < 0)); % note the indicator function for issuance cost


%% d) Stationary Distribution of Firms
% Consider now a world with many such firms and no entry or exit. 
% Specifically,
%   suppose that upon hitting the default threshold debt claims
%    are settled so b = 0.
% The restructured firm continues to operate but with capital,
%   k = 0 and the previous
%   productivity shock, z.

% Note now the setting has changed a little bit from question c). 
% All we have to change is the value function continuation in the bellman
% equation, where we don't set value of default to zero. We set if to value
% when k and b are 0.

% Notice this is the only difference from c). So I basically copied the
% codes for c) and made the change accordingly. 
% 1) This period, if the firm % default, instead of setting the value to 0,
% I set it to value at k=0.0000001 and b=0, and policy function entry is 
% set to equal to the entry when k=0.0000001 and b=0.
% 2) Next Period, if the firm % default, instead of setting the 
%    valueTomorrow
% to 0, I set it to value at k=0.0000001 and b=0.
% I use two kinds of endogenous grids to solve the problem - (a-,a,k,b) and
% (a,n). For the first method, we can easily make both the first and second
% adjustment, but for the second method where we only have today's networth
% information, we can only make the second adjustment and assume that
% today's firms are nondefaulters and the 0 networth denotes defaulters.

% Use multigrid method to speed up iteration
Nk = max(kGridLength);
kMax            = 10 * kSteadyState;
grid_b = curvspace(0,kMax,Nb,2)';

% Required matrices and vectors
% Dimensionality is k,b,a,aMinus
    
kPolicyIndex = zeros(kGridLength(1),Nb,Na,Na);
kPolicy = zeros(kGridLength(1),Nb,Na,Na);
bPolicyIndex = zeros(kGridLength(1),Nb,Na,Na);
bPolicy = zeros(kGridLength(1),Nb,Na,Na);

value   = zeros(kGridLength(1),Nb,Na,Na);
value0 = Vinit;%ones(kGridLength(1),Nb,Na,Na); % initial guess

tic
for i=1:length(kGridLength)
    
    grid_k           = curvspace(kMin,kMax,kGridLength(i),2)';
    % Since the profitFunction takes so much time, let's calculate it all
    % at once to retrieve later
    mANkByNa = repmat(grid_a',kGridLength(i),1);
    mKNkByNa = repmat(grid_k,1,Na);
    profitNkByNa = profitFunction(mANkByNa,mKNkByNa); % Nk by Na
    
    % Calculate Default Probability
    mCutOffValue = zeros(kGridLength(i),Nb);
    mDefaultProbability = zeros(kGridLength(i),Nb,Na);

    vDenominator = grid_k.^aalphaK .* (aalphaL * grid_k.^aalphaK / W).^(aalphaL/(1-aalphaL)) - W * (aalphaL * grid_k.^aalphaK/W).^(1/(1-aalphaL));

    for ib = 1:Nb
    %     vNumerator = max(0.000000001, (vBond(ib) - (1-ddelta)*k_grid)/(1-ttaoC));
        vNumerator =  (grid_b(ib) - (1-ddelta)*grid_k)/(1-ttaoC);

    %     vCutOffValue = (vNumerator ./v_Denominator).^(1-aalphaL); % I droped
    %     this line because it generates complex numbers

        vCutOffValue = vNumerator ./vDenominator;
        mCutOffValue(:,ib) = vCutOffValue; % cutoff value^(1-aalphaL) is cutoff productivity, but I don't use productivity per se in order to avoid complex numbers
        for ia = 1:Na
            for ik = 1:kGridLength(1)
                 temp= (vCutOffValue(ik) > grid_a.^(1/(1-aalphaL)));
                mDefaultProbability(ik,ib,ia) = sum(temp.* m_a_prob(ia,:)');
            end
        end
    end

    mRf = Rf*ones(kGridLength(1),Nb,Na);
    mCarryOnProbability = 1-mDefaultProbability;
    mIsDefaultNextPeriod = (mDefaultProbability==1);
    mRb = min(Rf/mCarryOnProbability,1000000);
        
    tolerance = 0.00001;
    iteration = 0;
    distance = 100;

    kPrime = repmat(grid_k,1,Nb); % Nk*Nb matrix
    bondPrime = repmat(grid_b',kGridLength(i),1); % Nk*Nb matrix
    
    tic
    while distance > tolerance
%         tic
        for ia=1:Na
            a = grid_a(ia);
            isDefaultNextPeriod = mIsDefaultNextPeriod(:,:,ia); % Nk*Nb matrix
            
            for ik=1:kGridLength(i)
                k = grid_k(ik);        
%                 labor = laborFunction(a,k);
%                 profit = profitFunction(a,k,labor);        
%                 profit = profitFunction(a,k); 
                profit = profitNkByNa(ik,ia);
                investment = investmentFunction(k,kPrime); % Nk*Nb matrix

                for ib = 1:Nb
                    bond = grid_b(ib);
                    
                    for iaMinus = 1:Na
                        RbMinus = mRb(ik,ib,iaMinus); % scalar
                        aMinus = grid_a(iaMinus);



                        if (1-ttaoC)*profit + (1-ddelta)*k <= bond % if default this period
                            value (ik,ib,ia,iaMinus)=value0(1,1,ia,iaMinus); % after restructuring, you continue to run the firm at square 1 with 0 capital and 0 bond
                            
                            kPolicyIndex(ik,ib,ia,iaMinus) = kPolicyIndex(1,1,ia,iaMinus) ; % policy function is the same as when you are at square 1
                            bPolicyIndex(ik,ib,ia,iaMinus) = bPolicyIndex(1,1,ia,iaMinus) ;

                            kPolicy(ik,ib,ia,iaMinus) = grid_k(kPolicyIndex(1,1,ia,iaMinus));
                            bPolicy(ik,ib,ia,iaMinus) = grid_b(bPolicyIndex(1,1,ia,iaMinus));

                        else % if not default this period

                            taxPayments = taxPaymentsFunction(k,bond,profit,RbMinus); % scalar
                            divident = dividentFunction(profit,investment,bond,bondPrime,RbMinus,taxPayments,isDefaultNextPeriod); % Nk*Nb matrix

                            valueTomorrow = zeros(kGridLength(i),Nb,Na);% k',b',a'
                            
                            for iaPrime = 1:Na % iterate over all possible states for tomorrow
                                aPrime = grid_a(iaPrime);
                                profitPrime = repmat(profitNkByNa(:,iaPrime),1,Nb);% Nk*Nb
                               
%                                 valueTomorrow(:,:,iaPrime) =  m_a_prob(ia,iaPrime) * ... % ��Ҫ����default֮��value����Ϊ0������Ϊset bond and capital to 0��value
%                                     (value0(:,:,iaPrime,ia).* ((1-ttaoC)*profitPrime + (1-ddelta)*kPrime > bondPrime)... % if not default next period
%                                     + value0(1,1,iaPrime,ia).* ((1-ttaoC)*profitPrime + (1-ddelta)*kPrime <= bondPrime));% if default next period
%                                 
                                valueTomorrow(:,:,iaPrime) =  m_a_prob(ia,iaPrime) *( ...% ��Ҫ����default֮��value����Ϊ0������Ϊset bond and capital to 0��value
                                    (1-isDefaultNextPeriod).*(...% if known today not for sure that default will happen next period
                                    value0(:,:,iaPrime,ia).* ((1-ttaoC)*profitPrime + (1-ddelta)*kPrime > bondPrime)... 
                                    + value0(1,1,iaPrime,ia).* ((1-ttaoC)*profitPrime + (1-ddelta)*kPrime <= bondPrime))...
                                    +isDefaultNextPeriod.* ... % if known today for sure WILL default next period
                                    (repmat(value0(:,1,iaPrime,ia),1,Nb).* ((1-ttaoC)*profitPrime + (1-ddelta)*kPrime > bondPrime)... % next period if not default 
                                    + value0(1,1,iaPrime,ia).* ((1-ttaoC)*profitPrime + (1-ddelta)*kPrime <= bondPrime)));% next period if default 
                            end
                            mExpectedValueTomorrow = sum(valueTomorrow,3); % sum by the third dimension to get a Nk*Nb matrix

                            x = divident + M * mExpectedValueTomorrow;

                            [rows,cols]=find(x==max(max(x)));

    %                 if (1-ttaoC)*profit + (1-ddelta)*k <= bond
    %                     value(ik,ib,ia) = 0;
    %                 else
    %                     value(ik,ib,ia) = x(max(rows),max(cols));
    %                 end

                            kPolicyIndex(ik,ib,ia,iaMinus) = min(rows);
                            bPolicyIndex(ik,ib,ia,iaMinus) = min(cols);

                            kPolicy(ik,ib,ia,iaMinus) = grid_k(min(rows));
                            bPolicy(ik,ib,ia,iaMinus) = grid_b(min(cols));
                            value(ik,ib,ia,iaMinus) = max(max(x));
                        end
                    end
                end
            end
        end

%         toc
        %hold on
        %mesh(aa, kk, value')
        distance = sum(sum(sum(sum(abs(value(:,:,:,:)-value0(:,:,:,:))))));
        value0 = value;
        iteration = iteration + 1;

        if mod(iteration,5) == 0
            display("iteration =    " + iteration + "   difference =   " +...
                distance )
        end
    end

    display("iteration =    " + iteration + "   difference =   " + distance + ". Converged")
    if i ~= length(kGridLength)
        value0 = interp1(grid_k,value,linspace(kMin, kMax, kGridLength(i+1)));
        value  = value0;
%         kPolicy = interp1(grid_k,kPolicy,linspace(kMin, kMax, kGridLength(i+1)));
        kPolicy         = zeros(kGridLength(i+1),Nb,Na,Na);
        kPolicyIndex    = zeros(kGridLength(i+1),Nb,Na,Na);
        bPolicy         = zeros(kGridLength(i+1),Nb,Na,Na);
        bPolicyIndex    = zeros(kGridLength(i+1),Nb,Na,Na);
    end
    
end

% I left the laptop in the office to let it run.
%    "iteration =    1636   difference =   9.9406e-06. Converged"

toc
save('resultD','value','kPolicy','bPolicy','kPolicyIndex','bPolicyIndex')

figure(6);
[bb,kk]=meshgrid(grid_b, grid_k);
mesh(bb, kk, value(:,:,1,round((Na+1)/2)));% yesterday's productivity is mean

for ia = 2:Na
    hold on;
    mesh(bb, kk, value(:,:,ia,round((Na+1)/2)));    
end

title('Value Under Different Shocks given mean $z^{-}$','interpreter','latex')
ylabel('Capital Stock $k$','interpreter','latex')
xlabel('Debt $b$','interpreter','latex')
zlabel('Value','interpreter','latex')
% zlim([-0.5*max(max(max(max(value)))),max(max(max(max(value))))])
zlim([0,max(max(max(max(value))))])
savefig('q1d_value_3D')

figure(7)
mesh(bb, kk, kPolicy(:,:,1,round((Na+1)/2)));% yesterday's productivity is mean

for ia = 2:Na
    hold on;
    mesh(bb, kk, kPolicy(:,:,ia,round((Na+1)/2)));    
end

title('Policy $k^\prime$ Under Different Shocks given mean $z^{-}$',...
    'interpreter','latex')
ylabel('Capital Stock $k$','interpreter','latex')
xlabel('Debt $b$','interpreter','latex')
zlabel('$k^\prime$','interpreter','latex')
savefig('q1d_kPolicy_3D')


figure(8)
mesh(bb, kk, bPolicy(:,:,1,round((Na+1)/2)));% yesterday's productivity is mean

for ia = 2:Na
    hold on;
    mesh(bb, kk, bPolicy(:,:,ia,round((Na+1)/2)));    
end

title('Policy $b^\prime$ Under Different Shocks given mean $z^{-}$',...
    'interpreter','latex')
ylabel('Capital Stock $k$','interpreter','latex')
xlabel('Debt $b$','interpreter','latex')
zlabel('$bond^\prime$','interpreter','latex')
savefig('q1d_bPolicy_3D')


%% Stationary Distribution
% Compute the stationary distribution of firms. 

distributionStationary0 = (1/(Nk*Nb*Na*Na))*ones(Nk,Nb,Na,Na);
distance=100;
tolerance=0.0000000001;
iteration=0;

while distance>tolerance
    distributionStationary1 = zeros(Nk,Nb,Na,Na);
    for ia=1:Na
        for iaMinus=1:Na
            for ib=1:Nb
                for ik=1:Nk
                    ikPrime = kPolicyIndex(ik,ib,ia,iaMinus);          
                    ibPrime = bPolicyIndex(ik,ib,ia,iaMinus);

                    prob = distributionStationary0(ik,ib,ia,iaMinus);
                    for iaPrime=1:Na
                        prob_aPrime = prob*m_a_prob(ia,iaPrime);
                        distributionStationary1(ikPrime,ibPrime,...
                            iaPrime,ia) = distributionStationary1(ikPrime,ibPrime,iaPrime,ia) + prob_aPrime;
                    end
                    
                end
            end
        end
    end
    
    distance=sum(sum(sum(sum(abs(distributionStationary0-distributionStationary1)))));
    distributionStationary0 = distributionStationary1;
    iteration = iteration + 1;
end

% Plot the distribution
figure(9);
[bb,kk]=meshgrid(grid_b, grid_k);
aMinusDescription = ["low","","medium","","high"];

for iaMinus = 1:Na
    subplot(1,Na,iaMinus);
    mesh(bb, kk, distributionStationary0(:,:,1,iaMinus));
    
    for ia = 2:Na
        hold on;
        mesh(bb,kk,distributionStationary0(:,:,ia,iaMinus));
    end
    title(['Distribution $z^{-}$ ',aMinusDescription(iaMinus)],...
        'interpreter','latex');
    ylabel('Capital Stock $k^\prime$','interpreter','latex')
    xlabel('Debt $b^\prime$','interpreter','latex')
    zlabel('Probability Mass','interpreter','latex')
end
    
savefig('q1d_stationary_distribution_3D')

%% Table
% Use the invariant distribution to construct
% a table reporting the cross-sectional average values of:
% (1) probability of default, p(��);
% (2) required return on risky bonds, Rb(��);
% (3) leverage ratio, b/k;
% (4) investment to capital ratio, i/k.
% (5) fraction of firms issuing equity;

% To do this question, I created 4D array to reduce the layer of loops.

mK4D=repmat(grid_k,1,Nb,Na,Na);% Nk*Nb*Na*Na matrix
mBond4D=repmat(grid_b',Nk,1,Na,Na); % Nk*Nb*Na*Na matrix

mA4D=repmat(grid_a',Na,1,Nb,Nk);
mA4D=permute(mA4D,[4,3,2,1]);% transform/reallocate the dimension to get
%       a Nk*Nb*Na*Na matrix

mAMinus4D=repmat(grid_a,1,Na,Nb,Nk);
mAMinus4D=permute(mAMinus4D,[4,3,2,1]);% transform/reallocate the dimension to get a Nk*Nb*Na*Na matrix

%% (1) probability of default, p(??);

% mLabor4D = laborFunction(mA4D,mK4D); % Nk*Nb matrix
% mProfit4D = profitFunction(mA4D,mK4D,mLabor4D);% Nk*Nb matrix
mProfit4D = profitFunction(mA4D,mK4D);% Nk*Nb matrix

mDefault4D = 1-nonDefaultFunction(mProfit4D,mK4D,mBond4D);% Nk*Nb 0-1 matrix

defaultProbability = sum(sum(sum(sum(mDefault4D.*distributionStationary0))));
fprintf('Average Default Probability is %2.10f\n', defaultProbability);

% mK = repmat(grid_k,1,Nb); % Nk*Nb matrix
% mBond = repmat(grid_b',Nk,1); % Nk*Nb matrix

% defaultAll = zeros(Nk,Nb,Na,Na);
% for iaMinus = 1:Na
%     for ia = 1:Na
%         a = grid_a(ia);
%         mLabor = laborFunction(a,mK); % Nk*Nb matrix
%         mProfit = profitFunction(a,mK,mLabor);% Nk*Nb matrix
%         mDefault = 1-nonDefaultFunction(mProfit,mK,mBond);% Nk*Nb 0-1 matrix
%         defaultAll(:,:,ia,iaMinus) = mDefault;
%     end
% end
% 
% defaultProbability = sum(sum(sum(sum(defaultAll.*distributionStationary0))));


%% (2) required return on risky bonds, Rb(��);
% riskyBondReturn = sum(sum(sum(sum((min(1000000000,Rf/(1-mDefault4D))).*distributionStationary0))));
riskyBondReturn = sum(sum(sum(sum((min(10000000000000,Rf/(1-mDefault4D))).*((mDefault4D~=1).*distributionStationary0)))));
fprintf('Required return on risky bonds on average is %2.8f\n', riskyBondReturn);

% riskyBondReturn = Rf/(1-defaultProbability);
% fprintf('Required return on risky bonds is %2.8f\n', riskyBondReturn);

%% (3) leverage ratio, b/k;
leverageRatio = sum(sum(sum(sum(mBond4D ./ mK4D   .*distributionStationary0))));
fprintf('Average leverage ratio is %2.8f\n', leverageRatio);

% leverageRatioAll = zeros(Nk,Nb,Na,Na);
% for iaMinus = 1:Na
%     for ia = 1:Na
%         a = grid_a(ia);
%         leverageRatioAll(:,:,ia,iaMinus) = mBond./mK;
%     end
% end
% 
% leverageRatio = sum(sum(sum(sum(leverageRatioAll.*distributionStationary0))));
% fprintf('Average leverage ratio is %2.8f\n', leverageRatio);

%% (4) investment to capital ratio, i/k.
investment2Capital=sum(sum(sum(sum((kPolicy./mK4D + 1 - ddelta).*...
    distributionStationary0))));
fprintf('Average investment to capital ratio is %2.8f\n', investment2Capital);

%% (5) fraction of firms issuing equity;
    % mIsDefaultNextPeriod4D = repmat(mIsDefaultNextPeriod,Na,1);�в�ͨ
mIsDefaultNextPeriod4D = zeros(kGridLength(1),Nb,Na,Na);
for iaMinus=1:Na
    aMinus = grid_a(iaMinus);
    mIsDefaultNextPeriod4D(:,:,:,iaMinus)=mIsDefaultNextPeriod;
end

mRbMinus4D = zeros(kGridLength(1),Nb,Na,Na);
for ia = 1:Na
    mRbMinus4D(:,:,ia,:)=mRb;
end
%     a = grid_a(ia);
%     for iaMinus = 1:Na
%         aMinus = grid_a(iaMinus);
%         mRbMinus4D(:,:,ia,iaMinus)=mRb(:,:,iaMinus);
%     end
% end
    
mInvestment4D = investmentFunction(mK4D,mK4D);
mTaxPayments4D = taxPaymentsFunction(mK4D,mBond4D,mProfit4D,mRbMinus4D);
mDivident4D = dividentFunction(mProfit4D,mInvestment4D,mBond4D,mBond4D,...
    mRbMinus4D,mTaxPayments4D,mIsDefaultNextPeriod4D);

mIsIssuingEquity = (mDivident4D<0);
fractionOfFirmsIssuingEquity = sum(sum(sum(sum(  mIsIssuingEquity .*distributionStationary0))));
fprintf('The fraction of firms issuing equity is %2.8f\n', fractionOfFirmsIssuingEquity);

T = table(defaultProbability,riskyBondReturn,leverageRatio,investment2Capital,fractionOfFirmsIssuingEquity)
save('table','T')
%     defaultProbability    riskyBondReturn       leverageRatio      investment2Capital    fractionOfFirmsIssuingEquity
%     __________________    ________________    _________________    __________________    ____________________________
% 
%   3.56994771855542e-09    1.01010100649501    0.539685179800697     2.00309413226414          0.030296293479109      

% Nk=20, Nb=20, Na=5
%     defaultProbability    riskyBondReturn    leverageRatio    investment2Capital    fractionOfFirmsIssuingEquity
%     __________________    _______________    _____________    __________________    ____________________________
% 
%         3.5699e-09            1.0101            0.65464             1.989                     0.030296      


%%% Run the linear regression to get beta and gamma
% first find Xs and ys with a simulation

mc = dtmc(m_a_prob);
Zetasii = simulate(mc,numsteps); %are the shocks
zminusi = 1; %take the lowst one as the initial;
k0i = 1;
b0i = 1;
%Xs = ones(numsteps,3); %contains the elements for regression
Xs = ones(numsteps,2); %contains the elements for regression
ys = ones(numsteps,1); %contains the dependent variable for regression
%Xs2 = ones(numsteps,4);
Xs2 = ones(numsteps,3);
ys2 = ones(numsteps,1);
k = kMin;
b = 0;
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
Xs2 = Xs2((0.2*numsteps):end,:);
ys2 = ys2((0.2*numsteps):end,:);
mdl = fitlm(Xs,ys);
betas = mdl.Coefficients.Estimate;
mdl2 = fitlm(Xs2,ys2);
gammas = mdl2.Coefficients.Estimate;

coefsRhoLambda = [betas' gammas'];
Vrenewed = value0;

%erase variables which can give trouble when recovered
rrhoUsed = rrho;
llambdaUsed = llambda;
clear rrho llambda;
%save data
if saveDoc ==1
    filename = '00_data_01';
    save(filename)
end

end
