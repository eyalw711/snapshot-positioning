function Zhat = sils_search(R,y,p)
%
% Zhat = sils_search(R,y,p) produces p optimal solutions to 
% the upper triangular integer least squares problem min_{z}||y-Rz|| 
% by a depth-first search algorithm.
%
% Inputs:
%    R - n-by-n real nonsingular upper triangular matrix
%    y - n-dimensional real vector
%    p - the number of optimal solutions with a default value of 1
%
% Output:
%    Zhat - n-by-p integer matrix (in double precision), whose j-th column 
%           is the j-th optimal solution, i.e., its residual is the j-th 
%           smallest, so ||y-R*Zhat(:,1)|| <= ...<= ||y-R*Zhat(:,p)||

% Main References:
% [1] X.-W. Chang, X. Yang, and T. Zhou, MLAMBDA: A Modified LAMBDA Method 
%     for Integer Least-squares Estimation, Journal of Geodesy, 79 (2005), 
%     pp. 552-565. 
% [2] A. Ghasemmehdi and E. Agrell, Faster Recursions in Sphere Decoding,
%     IEEE Transactions on Information Theory, 57 (2011), pp. 3530-3536.
%  

% Authors: Xiao-Wen Chang, www.cs.mcgill.ca/~chang
%          Tianyang Zhou, Xiangyu Ren
% Copyright (c) 2006-2016. Scientific Computing Lab, McGill University.
% October 2006. Last revision: June 2016.


% ------------------------------------------------------------------
% --------  Initialization  ----------------------------------------
% ------------------------------------------------------------------

n = size(R,1);

% Current point
z = zeros(n,1);

% c(k)=(y(k)-R(k,k+1:n)*z(k+1:n))/R(k,k)
c = zeros(n,1);

% d(k): left or right search direction at level k   
d = zeros(n,1); 

% Partial squared residual norm for z
% prsd(k) = (norm(y(k+1:n) - R(k+1:n,k+1:n)*z(k+1:n)))^2
prsd = zeros(n,1); 

% Store some quantities for efficiently calculating c
% S(k,n) = y(k),
% S(k,j-1) = y(k) - R(k,j:n)*z(j:n) = S(k,j) - R(k,j)*z(j), j=k+1:n
S = zeros(n,n);
S(:,n) = y;

% path(k): record information for updating S(k,k:path(k)-1) 
path = n*ones(n,1); 

% The level at which search starts to move up to a higher level
ulevel = 0; 

% The p candidate solutions (or points) 
Zhat = zeros(n,p); 

% Squared residual norms of the p candidate solutions
rsd = zeros(p,1); 

% Initial squared search radius
beta = inf; 

% The initial number of candidate solutions
ncand = 0;   

% ------------------------------------------------------------------
% --------  Search process  ----------------------------------------
% ------------------------------------------------------------------

c(n) = y(n) / R(n,n);
z(n) = round(c(n));
gamma = R(n,n) * (c(n) - z(n));
% Determine enumeration direction at level n
if c(n) > z(n)
    d(n) = 1;
else
    d(n) = -1;
end

k = n;

while 1
    % Temporary partial squared residual norm at level k
    newprsd = prsd(k) + gamma * gamma;

    if newprsd < beta
        if k ~= 1 % move to level k-1
            % Update path  
            if ulevel ~= 0 
                path(ulevel:k-1) = k;
                for j = ulevel-1 : -1 : 1
                     if path(j) < k
                           path(j) = k;
                     else
                         break;  % Note path(1:j-1) >= path(j)
                     end
                end
            end
            
            % Update S
            k = k - 1;
            for j = path(k) : -1 : k+1
                S(k,j-1) = S(k,j) - R(k,j) * z(j);
            end
            
            % Update the partial squared residual norm
            prsd(k) = newprsd;
            
            % Find the initial integer
            c(k) = S(k,k) / R(k,k);
            z(k) = round(c(k));
            gamma = R(k,k) * (c(k) - z(k));
            if c(k) > z(k) 
                d(k) = 1;
            else
                d(k) = -1;
            end
            
            ulevel = 0; 
            
        else % A new point is found, update the set of candidate solutions
            if ncand < p % Add the new point
                ncand = ncand + 1;
                Zhat(:,ncand) = z;
                rsd(ncand) = newprsd;
                if ncand == p
                    beta = rsd(p);
                end
            else % Insert the new point and remove the worst one
                i = 1;
                while i < p && rsd(i) <= newprsd
                    i = i + 1;
                end
                Zhat(:,i:p) = [z, Zhat(:,i:p-1)];
                rsd(i:p) = [newprsd; rsd(i:p-1)];
                beta = rsd(p);
            end

            z(1) = z(1) + d(1);
            gamma = R(1,1)*(c(1)-z(1));
            if d(1) > 0
                d(1) = -d(1) - 1;
            else
                d(1) = -d(1) + 1;
            end
        end
    else  
        if k == n % The p optimal solutions have been found
            break
        else  % Move back to level k+1
            if ulevel == 0
               ulevel = k;
            end
            k = k + 1; 
            % Find a new integer at level k  
            z(k) = z(k) + d(k);
            gamma = R(k,k) * (c(k) - z(k));
            if d(k) > 0
                d(k) = -d(k) - 1;
            else
                d(k) = -d(k) + 1;
            end
        end 
    end
end


