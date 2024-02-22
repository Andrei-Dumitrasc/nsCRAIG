% right preconditioning with blkdiag(M, I)
myG_matvecdata.type="MI_right_prec_gmres";
myG_matvecdata.A=A;
myG_matvecdata.M=decoM;

% RHS
Gb=[g;r];


for restart_as_nsCRAIG=restart_GMRES
    % restart
    if restart_as_nsCRAIG
        % allow GMRES only as much memory as GKB used, then restart;
        nsCRAIG_iter=length(relresnorms);
        mem_nsCRAIG=m+n*(nsCRAIG_iter+1);
        restart_m=floor(mem_nsCRAIG/(m+n));
    else
        restart_m=GMRES_max_iter; % full / unrestarted
    end
    % choose solver
    tic
    elapsed_its=0;
    relresnorms_G=1;
    start_x=zeros(m+n,1);
    beta0=norm(Gb);
    while relresnorms_G(end)>solve_tol
        [X_GMRES,relresnorms_G_part,V_myG,H_myG,beta0_myG] = ...
            myGMRES(myG_matvecdata,Gb,start_x,solve_tol,...
            min([restart_m,GMRES_max_iter-elapsed_its]),reortho);
        elapsed_its=elapsed_its+length(relresnorms_G_part)-1;
        % renormalize residual norm w.r.t. original problem            
        relresnorms_G=[relresnorms_G; relresnorms_G_part(2:end)*beta0_myG/beta0];
        start_x=X_GMRES;
    end
    'myGMRES time'
    if restart_as_nsCRAIG, 'with restarts', end 
    toc    

    u_GMRES=X_GMRES(1:m,end);
    p_GMRES=X_GMRES(m+1:end,end);
    % "reverse" preconditioning
    u_GMRES=decoM\u_GMRES;
    
    % convergence tracking
    gmres_name='GMRES';
    if restart_as_nsCRAIG
        gmres_name=[gmres_name,'(',num2str(restart_m),')'];
    else
        % 'memory GMRES / memory nsCRAIG'
        mem_GMRES=(m+n)*(length(relresnorms_G)+1);
        ax=gca;
        ax.Title.String=[ax.Title.String, '. memory ratio=', ...
            num2str(mem_GMRES/mem_nsCRAIG) ];
        
    end
    semilogy(relresnorms_G,'-p','DisplayName',gmres_name)
    du=u-u_GMRES;
    'difference between GMRES and GK final iterates'
    GKvs_GMRES_U=sqrt( (du'*M*du) / (exU'*M*exU) )
    GKvs_GMRES_P=norm(p-p_GMRES)/norm(exP)
end