function [R,piv,y] = qrmcp(B,y)
%
% [R,piv,y] = qrmcp(B,y) computes the QR factorization of B with 
%             minimum-column pivoting: 
%                  Q'BP = R (underdetermined B), 
%                  Q'BP = [R; 0] (underdetermined B)
%             and computes Q'*y. The orthogonal matrix Q is not produced. 
%
% Inputs:
%    B - m-by-n real matrix to be factorized
%    y - m-dimensional real vector to be transformed to Q'y
%
% Outputs:
%    R - m-by-n real upper trapezoidal matrix (m < n)
%        n-by-n real upper triangular matrix (m >= n)
%    piv - n-dimensional permutation vector representing P
%    y - m-vector transformed from the input y by Q, i.e., y := Q'*y 

% Authors: Xiao-Wen Chang, www.cs.mcgill.ca/~chang
%          Xiaohu Xie, Tianyang Zhou
% Copyright (c) 2006-2016. Scientific Computing Lab, McGill University.
% October 2006; Last revision: June 2016


[m,n] = size(B);

% Initialization
colNormB = zeros(2,n);
piv = 1:n;

% Compute the 2-norm squared of each column of B 
for j = 1:n
    colNormB(1,j) = (norm(B(:,j)))^2;
end

n_dim = min(m-1,n);

for k = 1 : n_dim
    % Find the column with minimum 2-norm in B(k:m,k:n)
    [~, i] = min(colNormB(1,k:n) - colNormB(2,k:n));
    q = i + k - 1;
    
    % Column interchange
    if q > k
        piv([k,q]) = piv([q,k]);
        colNormB(:,[k,q]) = colNormB(:,[q,k]);
        B(:,[k,q]) = B(:,[q,k]);
    end

    % Compute and apply the Householder transformation  I-tau*v*v'
    if norm(B(k+1:m,k)) > 0 % A Householder transformation is needed
	    v = B(k:m,k);
        rho = norm(v);     
	    if v(1) >= 0
            rho = -rho;
        end
        v(1) = v(1) - rho; % B(k,k)+sgn(B(k,k))*norm(B(k:n,k))
        tao = -1 / (rho * v(1));
        B(k,k) = rho;
        if m < n
           B(k+1:m,k) = 0;
        end
        B(k:m,k+1:n) = B(k:m,k+1:n) - tao * v * (v' * B(k:m,k+1:n));    
        % Update y by the Householder transformation
        y(k:m) = y(k:m,:) - tao * v * (v' * y(k:m));
    end
  
    % Update colnormB(2,k+1:n)
    colNormB(2,k+1:n) = colNormB(2,k+1:n) + B(k,k+1:n) .* B(k,k+1:n);
end

if m < n
   R = B;
else
   R = triu(B(1:n,1:n));
end