function [xhat,dual_poly_coeffs,Tu] =  ast_sdpt3(y,tau)
% AST_SDPT3 Solves AST for line spectral signals using SDPT3
%
% [xhat,dual_poly_coeffs,Tu] =  AST_SDPT3(y,tau)
%
% Denoises uniform samples of a mixture of complex sinusoids and
% returns the solution x of the atomic soft thresholding problem
%            1            2
% minimize  ---|| y - x || + tau ||x|| (AST)
%    x       2                        A
% where the atoms are complex sinusoids. AST is recast as an SDP
% solved using SDPT3.
%
% The function returns the optimal value x, and the dual polynomial coefficients
% q = (y - x)/mu. The dual polynomial coefficients specify a trigonometric
% polynomial whose absolute value reaches unity at the supporting
% frequencies of the solution.
% 
% This function requires the installation of SDPT3 package which is
% available at http://www.math.nus.edu.sg/~mattohkc/sdpt3.html
%
% See also AST_ADMM, AST_CVX

n = length(y);

% Write y' = [s; t; u1; uR; uI; xR; xI]
% Now, SDPT3 solves
%
% minimize b'*y
% subject to smat(At{1,1} y' + C{1,1}) is PSD
%            At{2,1} y'+ C{2,1} is SOC

b = 0.5*speye(4*n+1,1); % minimize s/2

En = speye(n);
E2 = speye(2);
PX = sparse([0,-1; 1 0]);

% Define various indices into y
i_s  = 1;
i_t  = 2;
i_u1 = 3;
i_uR_s = 3; % Add 1 to n-1
i_uR   = (1:n-1)+i_uR_s;
i_uI_s = 2+n; % Add 1 to n-1
i_uI   = (1:n-1)+i_uI_s;
i_xR_s = 2*n+1;
i_xR   = (1:n) + i_xR_s;
i_xI_s = 3*n+1;
i_xI   = (1:n) + i_xI_s;

% Semidefinite Constraint
blk{1,1} = 's'; blk{1,2} = 2*n+2;
% SOC Constraint
blk{2,1} = 'q'; blk{2,2} = 2*n+2;

% Semidefinite
At{1,1} = sparse(nchoosek(2*n+3,2),4*n+1);
At{1,1}(:,i_t) = svec(blk(1,:),kron(E2,sparse(diag([zeros(n,1);1])))); %t
At{1,1}(:,i_u1) = svec(blk(1,:),kron(E2,blkdiag(En,0))); %u1
for j=1:n-1
    At{1,1}(:,i_uR_s+j)  =svec(blk(1,:),kron(E2,blkdiag(toeplitz(En(:,j+1)),0))); %uR
    At{1,1}(:,i_uI_s+j)=svec(blk(1,:),kron(PX,blkdiag(toeplitz(-En(:,j+1),En(:,j+1)),0))); %uI
end
% For x
for j=1:n
    At{1,1}(:,i_xR_s+j)=svec(blk(1,:),kron(E2,[sparse(n,n),En(:,j);En(j,:),0]));
    At{1,1}(:,i_xI_s+j)=svec(blk(1,:),kron(PX,[sparse(n,n),En(:,j);-En(j,:),0]));
end

C{1,1} = sparse(2*n+2,2*n+2);
% SOC
At{2,1} = sparse(2*n+2,4*n+1);
At{2,1}(1,:) = [1,-tau,-tau*n,sparse(1,2*n-2),2*reshape(real(y),1,n),-2*reshape(imag(y),1,n)];
At{2,1}(2:n+1,i_xR)     = 2*speye(n);
At{2,1}(n+2:2*n+1,i_xI) = 2*speye(n);
At{2,1}(2*n+2,:) = -[1,-tau,-tau*n,sparse(1,2*n-2),2*reshape(real(y),1,n),-2*reshape(imag(y),1,n)];

C{2,1}      = sparse(2*n+2,1);
C{2,1}(1)   = 1;
C{2,1}(end) = 1;

opts = sqlparameters;
opts.printlevel = 0;
[~,X,yy]=sqlp(blk,At,C,b,opts);

% Primal Solution x
xhat = -yy(i_xR)+1i*yy(i_xI);

% Dual Solution
dual_poly_coeffs  = -X{1}(1:n,n+1)/tau-1i*X{1}(n+2:2*n+1,n+1)/tau;
dual_poly_coeffs = (y-xhat)/tau;

% Extract u, where Tu is the Toeplitz matrix of AST
u = -yy(3)*eye(n,1)-[0;yy(i_uR)]+1i*[0;yy(i_uI)];
Tu = toeplitz(u);
end