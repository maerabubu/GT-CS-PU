function [Xwh, mu, invMat, whMat] = whiten(X,epsilon)

% INPUT
% X: rows are the instances, columns are the features
% epsilon: small number to compensate for nearly 0 eigenvalue [DEFAULT =
% 0.0001]
%
% OUTPUT
% Xwh: whitened data, rows are instances, columns are features
% mu: mean of each feature of the orginal data
% invMat: the inverse data whitening matrix
% whMat: the whitening matrix

if ~exist('epsilon','var')
    epsilon = 0.0001;
end

mu = mean(X); 
X = bsxfun(@minus, X, mu);
A = X'*X;
[V,D,notused] = svd(A);
whMat = sqrt(size(X,1)-1)*V*sqrtm(inv(D + eye(size(D))*epsilon))*V';
Xwh = X*whMat;  
invMat = pinv(whMat);

end

