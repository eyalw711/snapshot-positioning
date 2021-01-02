function [X,Z] = smils(A,B,y,p)
%
% [X,Z] = smils(A,B,y,p) produces p pairs of optimal solutions to the 
% standard mixed integer least squares problem min_{x,z}||y-Ax-Bz||,
% where x and z are real and integer vectors, respectively.
%
% Inputs:
%    A - m-by-k real matrix
%    B - m-by-n real matrix, [A,B] has full column rank
%    y - m-dimensional real vector
%    p - The number of optimal solutions and its default value is 1
%
% Outputs:
%    X - k-by-p real matrix
%    Z - n-by-p integer matrix (in double precision). 
%        The pair {Xhat(:,j),Zhat(:,j)} is the j-th optimal solution
%        i.e., its residual is the j-th smallest, so
%        ||y-A*X(:,1)-B*Z(:,1)||<=...<=||y-A*X(:,p)-B*Z(:,p)||
%

% Subfunction: sils 

% Main reference
% X.-W. Chang and T. Zhou, MILES: MATLAB package for solving Mixed Integer 
% LEast Squares problems, GPS Solutions, 11 (2007), pp. 289-294. 

% Authors: Xiao-Wen Chang, www.cs.mcgill.ca/~chang 
%          Tianyang Zhou
% Copyright (c) 2006-2015. Scientific Computing Lab, McGill University.
% October 2006; Last revision: December 2015


% Check the input arguments
if nargin < 3 % input error
    error('Not enough input arguments!')
end

if nargin <= 3
    p = 1;
end

if p <= 0 % Input error
    error('Fourth input argument must be an integer bigger than 0!')
end

[m,k] = size(A);
n = size(B,2);
if m ~= size(B,1) || m ~= size(y,1) || size(y,2) ~= 1 % Input error
    error('Input arguments have a matrix dimension error!')
end

if rank([A,B]) < k + n
    error('Augumented matrix is rank defficient!')
end

% Compute the QR factorization of A
[Q,R] = qr(A);

% Compute the p optimal integer least squares solutions
Z = sils(Q(:,k+1:m)'*B, Q(:,k+1:m)'*y, p);

% Compute the corresponding real least squares solutions
X = R(1:k,:)\(Q(:,1:k)'*(y*ones(1,p) - B*Z));

