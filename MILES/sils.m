function X = sils(B,y,p)
%
% X = sils(B,y,p) produces p optimal solutions to the standard integer 
% least squares problem min_{x}||y-Bx||
%
% Inputs:
%    B - m-by-n real matrix with full column rank
%    y - m-dimensional real vector
%    p - number of optimal solutions and its default value is 1
%
% Output:
%    X - n-by-p integer matrix (in double precision), whose j-th column 
%        is the j-th optimal solution, i.e., its residual is the j-th 
%        smallest, ||y-B*X(:,1)|| <= ...<= ||y-B*X(:,p)||
%

% Subfunctions: sils_reduction, sils_search

% Main References:
% [1] X. Xie, X.-W. Chang, and M. Al Borno. Partial LLL reduction, 
%     Proceedings of IEEE GLOBECOM 2011, 5 pages.
% [2] X.-W. Chang, X. Yang, and T. Zhou, MLAMBDA: A Modified LAMBDA Method 
%     for Integer Least-squares Estimation, Journal of Geodesy, 79 (2005), 
%     pp. 552-565. 
% [3] A. Ghasemmehdi and E. Agrell, Faster Recursions in Sphere Decoding,
%     IEEE Transactions on Information Theory, 57 (2011), pp. 3530-3536. 

% Authors: Xiao-Wen Chang, www.cs.mcgill.ca/~chang 
%          Tianyang Zhou
% Copyright (c) 2006-2016. Scientific Computing Lab, McGill University.
% October 2006. Last revision: June 2016
 

% Check input arguments
if nargin < 2 % input error
    error('Not enough input arguments!')
end

if nargin < 3
    p = 1;
end

if p <= 0 % input error
    error('Third input argument must be an integer bigger than 0!')
end

[m,n] = size(B);

if rank(B) < n
	error('Matrix does not have full column rank!')
end

if m ~= size(y,1) || size(y,2) ~= 1  % Input error
    error('Input arguments have a matrix dimension error!')
end


% Reduction - reduce the problem to the triangular form
[R,Z,y] = sils_reduction(B,y);

% Search - find the p optimal solustions to the reduced problem
Zhat = sils_search(R,y(1:n),p);

% Perform the unimodual transformation to obtain the solutions to
%   the original problem
X = Z * Zhat;
