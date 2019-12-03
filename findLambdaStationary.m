function [lambda] = findLambdaStationary(lambdainit,matTransttn,IG,k,tol, maxiter)

d = 100;
nk = size(IG,1);
lambda = lambdainit;
iter = 0;
while d > tol && iter<maxiter
    iter = iter+1;
    for ii = 1:nk
        %ki = k(ii);
        IGi = IG == ii;
        lambda(ii,:) = sum(lambdainit.*IGi*matTransttn);
    end
    
    d = norm(lambdainit-lambda,2);
    lambdainit = lambda;
    if(mod(iter,20)==0)
        fprintf('lambdaStat Iteration = %d, Sup Diff = %2.8f\n', iter, d); 
    end
end

end