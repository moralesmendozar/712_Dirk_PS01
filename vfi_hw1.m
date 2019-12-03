function [V,G,inv,phi,fin] = vfi_hw1(k,A,QQ,r,theta1,theta2,delta,b0,b1p,b1m,d,tol,maxiter)
% VFI
% Anna Cororaton
% September 2014

nk      = length(k);
na      = length(A);
gamma   = theta1/(1-theta2);


% initialize
K       = repmat(k,1,nk);  
KP      = repmat(k',nk,1);  % k prime 
Vold    = zeros(nk,na);     % initial V
V       = zeros(nk,na);     % value function 
IG      = zeros(nk,na);     % index of policy


% dividends for all possible K
ind     = (KP > K); % check i>delta*k
F       = zeros(nk,nk,na);
for j = 1:na
    F(:,:,j) = A(j)*K.^gamma - KP + (1-delta-b0)*K ...
                - b1p.*(ind).*(((KP./K)-1).^2).*K ...
                - b1m.*(1-ind).*(((KP./K)-1).^2).*K;
end


% VFI
iter    = 0;
while d > tol && iter < maxiter 
    iter = iter + 1
    for i = 1:nk
        for j = 1:na 
            [V(i,j),IG(i,j)] = max(F(i,:,j)' + Vold*QQ(j,:)'./(1+r));
        end
    end
    d = abs(max(max(V-Vold)));
    Vold = V;
end


% policy
G       = k(IG); 

% investment
inv     = G - (1-delta)*repmat(k,1,na);

% adjustment costs
ind2    =  (G > repmat(k,1,na));
phi     = (b0)*repmat(k,1,na) ...
            + b1p.*ind2.*(((G./repmat(k,1,na))-1).^2).*repmat(k,1,na) ...
            + b1m.*(1-ind2).*(((G./repmat(k,1,na))-1).^2).*repmat(k,1,na);

% financing
fin     = repmat(A',nk,1).*repmat(k,1,na).^gamma - inv - phi;
