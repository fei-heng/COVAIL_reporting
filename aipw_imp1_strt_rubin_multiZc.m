function [eff_beta_hat,sig_eff_beta,p_values_t1,p_values_t2, F10, F50, F90]=aipw_imp1_strt_rubin_multiZc...
    (beta_ini,X,Z,Zc,Zm,Vpred,v,delta,sel,beta_n,wei,h,n0,n,K,Z_c,rep_n,a,b,aa,strata,tt,zpd2)

% add by Fei - 07/13/2023
% expand the function aipw_imp1_strt_rubin() to handle multiple Zc covariates

nstrt=length(unique([strata]));

eff_beta_hat = zeros(1, n0*beta_n);
% cov26 = zeros(1, n0); % extract cov(beta2, beta6) for the CRR analysis

Z_1c=[ones(n,1),Z_c];

% replace NaN by 0 in Z, renamed as ZZ
ZZ = Z;
ZZ(isnan(ZZ)) = 0;

% XX
XX = (repmat(X,1,n)>=repmat(X',n,1));

j=1;
for j=1:n0
    %j
    if (j==1)
        beta00=beta_ini(:,j);
    else
        beta00=beta11;
    end
    
    iter=1;
    Error=1;
    beta11=beta00;
    epsilon=0.001;
    
    while( Error>epsilon && iter < 20)
        %% predict Z
        
        predZ=zeros(n,beta_n);
        e_b1z1=zeros(n,1);
        e_b1z1_z1=zeros(n,1);
        e_b1z1_z1_2=zeros(n,1);
        for jj=1:nstrt
            Z_1c_sub=Z_1c;
            Z_1c_sub(strata~=jj-1,:)=NaN;
            
            
            eb=zeros(n,1);
            predZ_jj=ZZ;
            for k=1:beta_n
                eb=beta00(k)*ZZ(:,k)+eb; % Calculate beta*Z_ki(t)
                
                if k<=size(Zm,2)
                    [theta2,~,~] = regress(Z(:,k),Z_1c_sub);
                    predZ_jj(:,k)=sum(repmat(theta2',n,1).*Z_1c,2); %calculate E(Z|O)
                end
            end
            predZ_jj(strata~=jj-1,:)=0;
            eb(isnan(eb))=0;
            
            % e(b1z1)
            [theta2,~,~] = regress(exp(Zm*beta00(1)),Z_1c_sub);
            e_b1z1_jj=sum(repmat(theta2',n,1).*Z_1c,2);
            e_b1z1_jj(strata~=jj-1,:)=0;
            % e(b1z1)z1
            [theta2,~,~] = regress(exp(Zm*beta00(1)).*Zm,Z_1c_sub);
            e_b1z1_z1_jj=sum(repmat(theta2',n,1).*Z_1c,2);
            e_b1z1_z1_jj(strata~=jj-1,:)=0;
            % e(b1z1)z1^2
            [theta2,~,~] = regress(exp(Zm*beta00(1)).*Zm.*Zm,Z_1c_sub);
            e_b1z1_z1_2_jj=sum(repmat(theta2',n,1).*Z_1c,2);
            e_b1z1_z1_2_jj(strata~=jj-1,:)=0;
            
            predZ=predZ_jj+predZ;
            e_b1z1=e_b1z1_jj+e_b1z1;
            e_b1z1_z1=e_b1z1_z1_jj+e_b1z1_z1;
            e_b1z1_z1_2=e_b1z1_z1_2_jj+e_b1z1_z1_2;
        end
        
        e_b2z2 = zeros(n, size(Zc,2));
        for k = 1:size(Zc,2)
            e_b2z2(:,k)=exp(Z(:,k+size(Zm,2))*beta00(k+size(Zm,2))); % each component of e^(b2'z2)
        end
        % prod(..., dim=2): product of each row
        % note that z1 is for Zm, z2 is for Zc
        e_b2z2=prod([e_b2z2],2);  % e^(b2'z2)
        pred_eb0=prod([e_b1z1 , e_b2z2],2); %e^(b1z1+b2z2)
        pred_eb0(isnan(pred_eb0))=0;
        e_b1z1_b2z2_z1=prod([e_b1z1_z1 , e_b2z2],2); %e^(b1z1+b2z2)z1
        e_b1z1_b2z2_z1_2=prod([e_b1z1_z1_2 , e_b2z2],2); %e^(b1z1+b2z2)z1^2
        
        % e^(bz)z
        pred_eb1=[e_b1z1_b2z2_z1, repmat(pred_eb0,1,size(Zc,2)).*Zc];
        pred_eb1(isnan(pred_eb1))=0;
        % e^(bz)z^2
        pred_eb2=zeros(n*beta_n, beta_n);
        pred_eb2(1:n, 1)=e_b1z1_b2z2_z1_2;
        for icol = 2:beta_n
            pred_eb2(1:n, icol)=e_b1z1_b2z2_z1.*Zc(:,icol-1);
        end
        for irow = 2:beta_n
            pred_eb2(((irow-1)*n+1):(irow*n), 1)=e_b1z1_b2z2_z1.*Zc(:,irow-1);
        end
        for irow = 2:beta_n
            for icol = 2:beta_n
                pred_eb2(((irow-1)*n+1):(irow*n), icol)=pred_eb0.*Zc(:,irow-1).*Zc(:,icol-1);
            end
        end
        
        %         for k=1:beta_n
        %             pred_eb1_wd(:,(k-1)*n+1:k*n)=repmat(pred_eb1(:,k),1,n);
        %             pred_eb2_wd(:,(k-1)*n+1:k*n)=repmat(pred_eb2(:,k),1,n);
        %         end
        
        %%
        
        AS0=zeros(beta_n, n);
        AS1=zeros(beta_n, n);
        
        S0kk1=zeros(nstrt, n);
        S1kk1=zeros(nstrt, n);
        S2kk1=zeros(nstrt, n);
        
        for k=1:beta_n
            for kk=1:nstrt
                % sum/nansum: sum by column
                S0kk=(nansum(XX.*repmat((strata==kk-1).*exp(eb).*wei+...
                    (strata==kk-1).*pred_eb0.*(1-wei),1,n)))./...
                    sum(repmat(strata==kk-1,1,n));
                S1kk=(nansum(XX.*repmat((strata==kk-1).*exp(eb).*wei.*ZZ(:,k)+...
                    (strata==kk-1).*pred_eb1(:,k).*(1-wei),1,n)))...
                    ./sum(repmat(strata==kk-1,1,n));
                S0kk1(kk,:)=S0kk.*(strata==kk-1)';
                S1kk1(kk,:)=S1kk.*(strata==kk-1)';
            end
            % sum(...,1): sum by column
            AS0(k,:)=sum(S0kk1,1); % repeat S0kk1 beta_n times
            AS1(k,:)=sum(S1kk1,1);
        end
        
        
        Ubeta=0;
        Jacb=zeros(beta_n,beta_n);
        Jacb_rep=zeros(beta_n,beta_n);
        for rep=1:rep_n
            V=Vpred(:,rep);
            
            Ubeta_rep=nansum(repmat(K((V-v(j))/h)/h.*delta.*wei,1,beta_n)...
                .*(Z-(AS1./AS0)'))+nansum(repmat(K((V-v(j))/h)/h.*delta...
                .*(1-wei),1,beta_n).*(predZ-(AS1./AS0)'));
            
            for k=1:beta_n
                for l=1:beta_n
                    for kk=1:nstrt
                        S2kk=(nansum(XX.*repmat((strata==kk-1).*wei.*exp(eb).*ZZ(:,k).*ZZ(:,l)+...
                            (strata==kk-1).*(1-wei).*pred_eb2((l-1)*n+1:l*n,k),1,n)))./...
                            sum(repmat(strata==kk-1,1,n));
                        S2kk1(kk,:)=S2kk.*(strata==kk-1)';
                    end
                    AS2=sum(S2kk1,1);
                    %                     Jacb_rep(k,l)=-nansum(K((V-v(j))/h)/h.*delta.*wei.*...
                    %                         (AS2./AS0(k,:)-(AS1(k,:)./AS0(k,:).*AS1(l,:)./AS0(l,:)))')-...
                    %                         nansum(K((V-v(j))/h)/h.*delta.*(1-wei).*...
                    %                         (AS2./AS0(k,:)-(AS1(k,:)./AS0(k,:).*AS1(l,:)./AS0(l,:)))');
                    Jacb_rep(k,l)=-nansum(K((V-v(j))/h)/h.*delta.*...
                        (AS2./AS0(k,:)-(AS1(k,:)./AS0(k,:).*AS1(l,:)./AS0(l,:)))');
                end
            end
            
            Ubeta=Ubeta+Ubeta_rep;
            Jacb=Jacb+Jacb_rep;
        end
        if (rcond(Jacb)<0.00000001)
            break;
        end
        try  beta11=beta00-pinv(Jacb)*Ubeta';
            Error=nansum(abs(beta11-beta00));
            iter=iter+1;
            beta00=beta11;
        catch
            break;
        end
    end
    eff_beta_hat(1,j:n0:beta_n*n0)=beta11';
    
    %% variance correction using Rubin's idea
    
    Ubeta_m=zeros(rep_n, beta_n);
    
    for rep=1:rep_n
        V=Vpred(:,rep);
        
        Ubeta_rep=nansum(repmat(K((V-v(j))/h)/h.*delta.*wei,1,beta_n)...
            .*(Z-(AS1./AS0)'))+nansum(repmat(K((V-v(j))/h)/h.*delta...
            .*(1-wei),1,beta_n).*(predZ-(AS1./AS0)'));
        
        Ubeta_m(rep,:)=Ubeta_rep;
    end
    
    Var_rubin = (1+1/rep_n)/(rep_n-1) * (Ubeta_m') * Ubeta_m;
    
    
    %% variance estimation
    A=Jacb;
    
    % new
    temp=zeros(beta_n,beta_n);
    for rep=1:rep_n
        V=Vpred(:,rep);
        
        Var1a_rep=repmat(K((V-v(j))/h)/h.*delta,1,beta_n).*(Z-(AS1./AS0)') ; % N part for wei
        Var1b_rep=repmat(K((V-v(j))/h)/h.*delta,1,beta_n).*(predZ-(AS1./AS0)'); % N part for 1-wei
        Var1a_rep(isnan(Var1a_rep))=0;
        Var1b_rep(isnan(Var1b_rep))=0;
        
        Var2_rep=zeros(n,beta_n); % lambda part for wei
        Var3_rep=zeros(n,beta_n); % lambda part for 1-wei
        
        Var2kk=zeros(nstrt,n);
        Var3kk=zeros(nstrt,n);
        for k=1:beta_n
            for kk=1:nstrt
                S0kk=(nansum(XX.*repmat((strata==kk-1).*exp(eb).*wei+...
                    (strata==kk-1).*pred_eb0.*(1-wei),1,n)))./...
                    sum(repmat(strata==kk-1,1,n));
                %                 dlambdakk=1./S0kk.*delta'.*...%wei'.*... % should we delete wei here? - Fei 07/14/2023
                %                     nansum((repmat(X.*(strata==kk-1).*delta,1,n)==repmat(X',n,1)).*...
                %                     repmat(strata==kk-1,1,n)')./nansum(repmat(strata==kk-1,1,n));
                dlambdakk=1./S0kk.*delta'.*(strata==kk-1)'./nansum(repmat(strata==kk-1,1,n));
                dlambdakk(isnan(dlambdakk))=0;
                Var2a= nansum( (repmat(Z(:,k),1,n) - repmat(AS1(k,:)./AS0(k,:),n,1)).*repmat(exp(eb),1,n).*XX.*...
                    repmat(K((V'-v(j))/h)/h.*dlambdakk,n,1),2);
                Var2kk(kk,:)=Var2a'.*(strata==kk-1)' ;
                
                Var3a= nansum((repmat(predZ(:,k),1,n)- repmat(AS1(k,:)./AS0(k,:),n,1)).*pred_eb0.*XX.*...
                    repmat(K((V'-v(j))/h)/h.*dlambdakk,n,1),2);
                Var3kk(kk,:)=Var3a'.*(strata==kk-1)' ;
            end
            Var2_rep(:,k)=nansum(Var2kk,1)';
            Var3_rep(:,k)=nansum(Var3kk,1)';
        end
        
        %         mean((Var1a_rep-Var2_rep))
        %         mean((Var1b_rep-Var3_rep))
        
        % Q^2
        temp = temp + (repmat(wei,1,beta_n).*(Var1a_rep-Var2_rep)+...
            repmat(1-wei,1,beta_n).*(Var1b_rep-Var3_rep))'*...
            (repmat(wei,1,beta_n).*(Var1a_rep-Var2_rep)+...
            repmat(1-wei,1,beta_n).*(Var1b_rep-Var3_rep));
    end
    temp=temp.*rep_n;
    
    Sigma_temp = pinv(A)'* temp *pinv(A)...
        + pinv(A/rep_n)'*Var_rubin*pinv(A/rep_n);
    %cov26(1,j) = Sigma_temp(2,beta_n);
    sig_eff_beta(1,j:n0:beta_n*n0)=sqrt(diag(Sigma_temp));
end


%% new test 12262019
G=3;
if (G==4)
    v_ind=[4,10,16,22];
end
if (G==3)
    v_ind=[1,13,25];
end

q=1;
D=10000;

k=1;
for k=1:beta_n
    sig_Lbeta=sig_eff_beta(1, n0*(k-1)+v_ind);
    Lbeta=eff_beta_hat(1, n0*(k-1)+v_ind);
    
    T1_a_temp=(Lbeta.^2)./(sig_Lbeta.^2);
    T1_a1=max(T1_a_temp);
    T1_a2=sum(T1_a_temp);
    
    T1_a_tuta=chi2rnd(1,D,G*q);
    T1_a1_tuta=max(T1_a_tuta,[],2);
    T1_a2_tuta=sum(T1_a_tuta,2);
    %T1_a2_tuta=chi2rnd(G*q,D,1);
    
    T1_m_temp=(Lbeta)./(sig_Lbeta);
    T1_m1=min(T1_m_temp);
    T1_m2=sum(T1_m_temp);
    
    T1_m_tuta=normrnd(0,sqrt(q),D,G);
    T1_m1_tuta=min(T1_m_tuta,[],2);
    T1_m2_tuta=sum(T1_m_tuta,2);
    %T1_m2_tuta=normrnd(0,sqrt(G*q),D,1);
    
    pv_a1(k)=mean(T1_a1_tuta>T1_a1);
    pv_a2(k)=mean(T1_a2_tuta>T1_a2);
    pv_m1(k)=mean(T1_m1_tuta<T1_m1);
    pv_m2(k)=mean(T1_m2_tuta<T1_m2);
    
    A=zeros(G-1,G);
    for Ai=1:G-1
        for Aj=1:G
            if Ai==Aj
                A(Ai,Aj)=-1;
            elseif Ai==(Aj-1)
                A(Ai,Aj)=1;
            end
        end
    end
    
    Qbeta=A*Lbeta';
    Cov_Lbeta=diag(sig_Lbeta.^2);
    Cov_Qbeta=A*Cov_Lbeta*A';
    
    T2_a2=Qbeta'*pinv(Cov_Qbeta)*Qbeta;
    
    T2_a_tuta=chi2rnd(1,D,(G-1)*q);
    %T2_a2_tuta=sum(T2_a_tuta,2);
    T2_a2_tuta=chi2rnd((G-1)*q,D,1); % Fei 06172024
    
    Cov_Qbeta_sqrt = sqrtm(Cov_Qbeta);
    T2_m2=sum(pinv(Cov_Qbeta_sqrt)*Qbeta);
    
    T2_m_tuta=normrnd(0,sqrt(q),D,G-1);
    T2_m2_tuta=sum(T2_m_tuta,2);
    
    % -------- delete T2_a1, T2_m1 from the results --------
    T2_a_temp=((Lbeta(2:G)-Lbeta(1:G-1)).^2)./(sig_Lbeta(2:G).^2+sig_Lbeta(1:G-1).^2);
    T2_a1=max(T2_a_temp);
    %T2_a2=sum(T2_a_temp);
    
    T2_a_tuta=chi2rnd(1,D,(G-1)*q);
    T2_a1_tuta=max(T2_a_tuta,[],2);
    %T2_a2_tuta=sum(T2_a_tuta,2);
    
    T2_m_temp=(Lbeta(2:G)-Lbeta(1:G-1))./sqrt(sig_Lbeta(2:G).^2+sig_Lbeta(1:G-1).^2);
    T2_m1=max(T2_m_temp);
    %T2_m2=sum(T2_m_temp);
    
    T2_m_tuta=normrnd(0,sqrt(q),D,G-1);
    T2_m1_tuta=max(T2_m_tuta,[],2);
    %T2_m2_tuta=sum(T2_m_tuta,2);
    % ------------------------------------------------------
    
    pv_a1_t2(k)=mean(T2_a1_tuta>T2_a1);
    pv_a2_t2(k)=mean(T2_a2_tuta>T2_a2);
    pv_m1_t2(k)=mean(T2_m1_tuta>T2_m1);
    pv_m2_t2(k)=mean(T2_m2_tuta>T2_m2);
end


p_values_t1=[pv_a1;pv_a2;pv_m1;pv_m2];
p_values_t2=[pv_a1_t2;pv_a2_t2;pv_m1_t2;pv_m2_t2];


%%%%%%%%%%%%%%%% Cumulative incidence function  %%%%%%%%%%%%%%%%%%
kk=1;
F10=zeros(nstrt,n0);
F50=zeros(nstrt,n0);
F90=zeros(nstrt,n0);
for kk=1:nstrt
    
    n2=30;
    v2=linspace(a,b,n2);
    d_v=diff(v2);
    
    eff_beta_u=zeros(n2,size(Z,2));
    for k=1:size(Z,2)
        rr=ksrmv(v',eff_beta_hat((k-1)*length(v)+(1:length(v)))',v(2)-v(1),v2');
        eff_beta_u(:,k)=rr.f;  %smooth beta over v2
    end
    
    A0_tx=0;
    A0_tx2=0;
    
    
    rep=1;
    for rep=1:rep_n
        A0_tx_rep=[];
        A0_tx_rep2=[];
        V=Vpred(:,rep);
        j=1;
        for j=1:n2
            
            e_b1z1=zeros(n,1);
            
            
            
            jj=1;
            for jj=1:nstrt
                Z_1c_sub=Z_1c;
                Z_1c_sub(strata~=jj-1,:)=NaN;
                
                eb=zeros(n,1);
                for k=1:beta_n
                    eb=eff_beta_u(j,k)*ZZ(:,k)+eb; % Calculate beta*Z_ki(t)
                end
                eb(isnan(eb))=0;
                
                [theta2,~,~] = regress(exp(Zm*eff_beta_u(j,1)),Z_1c_sub);
                e_b1z1_jj=sum(repmat(theta2',n,1).*Z_1c,2);
                e_b1z1_jj(strata~=jj-1,:)=0;
                e_b1z1=e_b1z1_jj+e_b1z1;
            end
            for k = 1:size(Zc,2)
                e_b2z2(:,k)=exp(Z(:,k+size(Zm,2))*eff_beta_u(j,k+size(Zm,2))); % each component of e^(b2'z2)
            end
            e_b2z2=prod([e_b2z2],2);  % e^(b2'z2)
            pred_eb0=prod([e_b1z1 , e_b2z2],2); %e^(b1z2+b2z2)
            
            Z_1c=[ones(n,1),Z_c];
            
            S0kk=(nansum(XX.*repmat((strata==kk-1).*exp(eb).*wei,1,n)+...
                XX.*repmat((strata==kk-1).*pred_eb0.*(1-wei),1,n)))./nansum(repmat(strata==kk-1,1,n));
            dlambdakk=1./S0kk./sum(repmat(strata==kk-1,1,n));
            A0_tx_rep(:,j)=nansum(XX.*repmat((strata'==kk-1),n,1).*...
                repmat(delta',n,1).*repmat( K((V'-v2(j))/h)/h.*dlambdakk,n,1),2);
            A0_tx_rep2(:,j)=nansum(XX.*(repmat(v2(j),n,n)>=repmat(V',n,1)).*...
                repmat((strata'==kk-1),n,1).*repmat(delta',n,1).*repmat(dlambdakk,n,1 ) ,2);
            
        end
        A0_tx=A0_tx+A0_tx_rep/rep_n; %in the paper (20)
        A0_tx2=A0_tx2+A0_tx_rep2/rep_n; %in the paper (20)
    end
    
    qtl=[quantile(Z(:,1),0.1),quantile(Z(:,1),0.5), quantile(Z(:,1),0.9) ];  % quantiles to be estimate for Z1
    %   qtl=[0,quantile(Z(:,1),0.5), quantile(Z(:,1),0.9)];  % quantiles to be estimate for Z1
    %   qtl=[-1,0, 1 ];
    
    F_list=zeros(3,n0);
    for iqt=1:3
        Z_pd=[qtl(iqt),zpd2];
        
        j=1;
        eb00a=zeros(1,n2);
        for j=1:n2
            eb00=0;
            for k=1:size(Z,2)
                eb00=eff_beta_u(j,k)*Z_pd(k)+eb00;
            end
            eb00a(j)=eb00;
        end
        e_A=sum(A0_tx.*d_v(1).*exp(eb00a),2); %in the paper (21) in the exp part
        
        rep=1;
        for rep=1:rep_n
            V=Vpred(:,rep);
            F_rep=zeros(1,n0);
            for j=1:n0
                eb0=0;
                for k=1:beta_n
                    eb0=eff_beta_hat(j+n0*(k-1))*Z_pd(k)+eb0;
                end
                
                e_b1z1=zeros(n,1);
                jj=1;
                for jj=1:nstrt
                    Z_1c_sub=Z_1c;
                    Z_1c_sub(strata~=jj-1,:)=NaN;
                    
                    eb=zeros(n,1);
                    for k=1:beta_n
                        eb=eff_beta_hat(j+n0*(k-1))*ZZ(:,k)+eb; % Calculate beta*Z_ki(t)
                    end
                    eb(isnan(eb))=0;
                    
                    [theta2,~,~] = regress(exp(Zm*eff_beta_hat(j)),Z_1c_sub);
                    e_b1z1_jj=sum(repmat(theta2',n,1).*Z_1c,2);
                    e_b1z1_jj(strata~=jj-1,:)=0;
                    e_b1z1=e_b1z1_jj+e_b1z1;
                    
                end
                for k = 1 : size(Zc,2)
                    e_b2z2(:,k)=exp(Z(:,k+size(Zm,2))*eff_beta_hat(size(Zm,2)*n0+j+n0*(k-1))); % each component of e^(b2'z2)
                end
                e_b2z2=prod([e_b2z2],2);  % e^(b2'z2)
                pred_eb0=prod([e_b1z1 , e_b2z2],2); %e^(b1z2+b2z2)
                
                
                S0kk=(nansum(XX.*repmat((strata==kk-1).*exp(eb).*wei,1,n)+...
                    XX.*repmat((strata==kk-1).*pred_eb0.*(1-wei),1,n)))./nansum(repmat(strata==kk-1,1,n));
                dlambdakk=1./S0kk./sum(repmat(strata==kk-1,1,n));
                F_rep(j)= nansum( (X<=tt).*(strata==kk-1).*delta.*exp(-e_A).* K((V-v(j))/h)/h .*dlambdakk')*exp(eb0) ; %in the paper (22)
            end
            
            F_list(iqt,:)=F_list(iqt,:)+F_rep/rep_n;
        end
    end
    
    F10(kk,:)=F_list(1,:);
    F50(kk,:)=F_list(2,:);
    F90(kk,:)=F_list(3,:);
end

