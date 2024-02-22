function [V,Q,B,H,residual,relresnorms,vel,pres] = nsCRAIG(A,b,maxit,tol,M,decoM,ortho_type)

% Generalized Golub-Kahan bidiagonalization for saddle-point matrices with
% a nonsymmetric leading block. 
% See https://doi.org/10.48550/arXiv.2310.06952
% ========
% [V,Q,B,H,residual,relresnorms,vel,pres] = nsGKB(A,b,m,tol,M,decoM,ortho_type)
% ========
% INPUT
% A - rectangular matrix  ( (1,2)-block of saddle-point matrix )
% b - start vector / lower RHS / initial residual
% maxit - maximum # of steps
% tol - tolerance
% M - Positive definite matrix, redefines orthogonality for the left
%     vectors (elliptic sing. vectors)
% decoM - M's decomposition object. Optional.
% ortho_type - 0: Gram-Schmidt;
%         >0: iterative modified Gram-Schmidt, run this many times.
%
% OUTPUT
% V - left orthogonal matrix
% Q - right orthogonal matrix
% B - upper bidiagonal matrix
% H - upper Hessenberg matrix
% residual - last residual vector
% relresnorms - relative residual vector norms (convergence history). 
%               rho in the paper
% vel - velocity (primal solution variable)
% pres - pressure (dual solution variable)

% allocation, preparation, etc
[m,n]=size(A);
V=zeros(m,maxit);
Q=zeros(n,maxit);
H=zeros(maxit+1,maxit);
B=zeros(maxit,maxit);
relresnorms=zeros(maxit+1,1);
x=zeros(maxit+1,1);

% normalize starting vector
be1=norm(b);
if be1~=1
    q=b/be1;
end
Q(:,1)=q;

% find first left vector
if isempty(decoM)
    decoM=M;
end
v=decoM\(A*q);

al=sqrt(v'*M*v);
B(1,1)=al;
v=v/al;
V(:,1)=v;

relresnorms(1)=1;
x(1)=be1/al;
% main loop
for j=1:maxit
    % new right vector
    q=A'*v;
    
    [h,q]=ortho(q,Q(:,1:j),ortho_type); % orthogonalization
    H(1:j,j)=h;    
    
    % compute the quantities for the stopping criterion
    be=norm(q);
    solve_res_norm=be*abs(x(j));
    rrn=solve_res_norm/be1;
    relresnorms(j+1)=rrn;
    
    if relresnorms(j+1) <= tol || j==maxit
        % stop iterating or ...
        break
    else
        % ... continue iterating
        q=q/be;
        Q(:,j+1)=q;
        H(j+1,j)=be;
        
        % new left vector
        w=decoM\(A*q);        
        w=w-V(:,j)*be;
        al=sqrt(w'*M*w);
        v=w/al;
        
        V(:,j+1)=v;
        
        B(j+1,j+1)=al;
        B(j,j+1)=be;
        x(j+1)=-be*x(j)/al;
        
    end
end
residual=q;

% trim
V(:,j+1:end)=[];
Q(:,j+1:end)=[];
relresnorms(j+2:end)=[];
B(:,j+1:end)=[]; B(j+1:end,:)=[];
H(:,j+1:end)=[]; H(j+1:end,:)=[];

% compute solutions (final iterates)
HS=H*B;
pres=-Q*(HS\(be1*[1; zeros(j-1 ,1)]));
vel=-decoM\(A*pres);
