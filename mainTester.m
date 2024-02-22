% Tests the nonsymmetric CRAIG algorithm (nsCRAIG) on two linearized Navier-Stokes
% problems, and optionally performs comparisons with other approaches 
% and algorithms (FOM and GMRES).

clear
clc

cla
crt_fig=gcf();
crt_fig.WindowStyle='docked';

subplotidx=1;

% load problem
% Choices used in the paper: 'NSStepGrid5Visc1Ov100', 'NSCavityGrid5Visc1Ov200'
for pbtype={ 'NSStepGrid5Visc1Ov100', 'NSCavityGrid5Visc1Ov200'}
pb.type=pbtype{1}
if ~isempty(regexpi(pb.type,'Visc','once'))
    load(['test_matrices\',pb.type]);
    M=W; clear W;
    % rescale the problem
    nMexU=sqrt(exU'*M*exU);    
    g=M*exU+A*exP; % recompute RHS
    r=A'*exU;
end

% parameters
[m,n]=size(A);
for solve_tol=[1e-3]        solve_tol
for reortho=[1]   reortho    % 0: Gram-Schmidt. >0: Modified GS, repeated.
nsCRAIG_max_iter= min(n,200)

% optional parts on/off (1/0)
compareWithFOM=1; % FOM on the Schur complement equation

compareWithGMRES=1; % GMRES on the entire system
GMRES_max_iter=min(m+n,1000);
restart_GMRES=[1 0] ; % restart 1, no restart 0. 
% Restart parameter chosen automatically in order to limit GMRES 
% to using at most as much memory as nsCRAIG has used.

% ========== main section ==========

decoM=decomposition(M,'lu');
% initial transformation
u0=decoM\g;
% start vector / initial residual / first right vector
start_vec=r-A'*u0;

% nonsymmetric nsCRAIG
tic
[V,Q,B,H,residual,relresnorms,u,p]=nsCRAIG(A,start_vec,nsCRAIG_max_iter,...
    solve_tol,M,decoM,reortho);
'nsCRAIG time'
toc

u=u+u0;
d=u-exU;
'final solver error norms (relative)'
errU=sqrt( (d'*M*d) / (exU'*M*exU) )
errP=norm(p-exP)/norm(exP)

% plotting
subplot(1,2,subplotidx); subplotidx=subplotidx+1;
hold on
ca=gca;
ca.YScale='log';
ylabel('residual norm')
xlabel('iterations')

semilogy(relresnorms,'-*','DisplayName','nsCRAIG')
title(pb.type)
legend show

%%
% ======== optional parts ========

if compareWithFOM
    % reference problem (FOM on Schur complement)    
    b=-start_vec;
    fom_maxiter= nsCRAIG_max_iter
    matvecdata.type="schur_comp";
    matvecdata.A=A; matvecdata.M=decoM;
    
    tic
    [P_fom,relresnorms_fom,V_fom,H_fom,Y_fom] = myFOM(matvecdata,b,...
        zeros(n,1), solve_tol, fom_maxiter, reortho);
    'FOM time'
    toc    
    
    semilogy(relresnorms_fom,'-d','DisplayName','FOM Schur Comp.')
    'relative difference between FOM and GK final iterates (pressure)'
    GKvsFOM=norm(p-P_fom(:,end))/norm(exP)
end

if compareWithGMRES
    compareWithGMRES_script
end

end
end
end