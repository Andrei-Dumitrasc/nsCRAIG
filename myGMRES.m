function [x,relresnorms,V,H,beta0]=myGMRES(matvecdata,b,x0,tol,maxit,ortho_type)

% GMRES (Generalized Minimum RESidual) based on the Arnoldi process.
% See e.g. Y. Saad, Iterative methods for sparse linear systems,
%   SIAM, 2003, https://doi.org/10.1137/1.6609780898718003
%
%  [x,relresnorms,V,H,beta0]=myGMRES(matvecdata,b,x0,tol,maxit,ortho_type)
%
% INPUT:
% matvecdata - struct with the matrix/matrices used in the multiplication;
% b - right hand side;
% x0 - initial guess;
% tol - relative value for the stopping criterion;
% maxit - maximum allowed number of iterations;
% ortho_type - orthogonalization method:
%           0 - Gram-Schmidt;
%           >=1 - modified Gram-Schmidt, repeated this many times;
%
% OUTPUT:
% x - final approximate solution x
% relresnorms - relative residual 2-norms
% V - orthogonal basis vectors
% H - upper Hessenberg matrix
% beta0 - initial residual norm;

% allocation, preparation, etc
n=length(x0);
V=zeros(n,maxit);
H=zeros(maxit+1,maxit);
HR=zeros(maxit+1,maxit);
G=zeros(2,maxit);
relresnorms=zeros(maxit+1,1);

% first step
matvecdata.b=x0;
r=b-custom_matvec(matvecdata);
beta0=norm(r); relresnorms(1)=1;
v=r/beta0; V(:,1)=v;

% Arnoldi iteration loop
for i=1:maxit
    matvecdata.b=V(:,i);
    w=custom_matvec(matvecdata); % matrix vector product
    [h,w]=ortho(w,V(:,1:i),ortho_type); % orthogonalize
    beta=norm(w);
    H(i+1,i)=beta;
    H(1:i,i)=h;
    V(:,i+1)= w/beta;

    % residual norm
    HR(1:i,i)=h;
    % apply previous rotations
    for ri=1:i-1
        htemp=HR(ri,i);
        HR(ri,i)=G(1,ri)*htemp+G(2,ri)*HR(ri+1,i);
        HR(ri+1,i)=-G(2,ri)*htemp+G(1,ri)*HR(ri+1,i);
    end
    hb=[HR(i,i);beta];
    rho=norm(hb);
    G(1:2,i)=hb/rho;  % new Givens rotation c and s
    HR(i,i)=rho;
    relresnorms(i+1)=abs(relresnorms(i)*G(2,i));

    if relresnorms(i+1) < tol
        break
    end
end

% trim
V(:,i+2:end)=[];
H(:,i+1:end)=[]; H(i+2:end,:)=[];
relresnorms(i+2:end)=[];

% solver
bet=beta0;
Hb=zeros(i+1,1);
% apply all rotations to the RHS
for ri=1:i
    Hb(ri)=G(1,ri)*bet;
    Hb(ri+1)=-G(2,ri)*bet;
    bet=Hb(ri+1);
end
Hb(i+1)=[];

x=x0+V(:,1:i)*(HR(1:i,1:i)\Hb);
