
%% Data preparation
df_marks=readtable('covail_rbd_dms_escape_marks_v6.csv');
% disp(df_marks.Properties.VariableNames);

% remove the second swap info for ptid = 933
% {'COV.00933|SWAB2|04/04/23|XBB.1.5.33|XBB.1.5'}
df_marks(122, :) = [];

% Define the format of dates
opts = detectImportOptions('covail_data_processed_20250524.csv');

% Specify the correct format for the Date column
opts = setvaropts(opts, 'primary_vax1_date', 'InputFormat', 'MM/dd/yy');
opts = setvaropts(opts, 'infect_date1', 'InputFormat', 'MM/dd/yy');
opts = setvaropts(opts, 'studydose1date', 'InputFormat', 'yyyy-MM-dd');
opts = setvaropts(opts, 'studydose2date', 'InputFormat', 'MM/dd/yy');
opts = setvaropts(opts, 'Actual_visit_date_15_2', 'InputFormat', 'yyyy-MM-dd');

df = readtable('covail_data_processed_20250524.csv', opts);
% disp(df.Properties.VariableNames);

% remove sequence information from noncases
noncase_ptid = df.Ptid(df.EventIndPrimaryD15 == 0);
keep_idx = ~ismember(df_marks.ptid, noncase_ptid);
df_marks_filtered = df_marks(keep_idx, :);

df_filtered = df(~ismember(df.arm, [3, 10, 11, 13, 14, 15]) & ...
    df.TwophasesampIndD15 == 1 & ...
    df.Immunemarkerset == 1, :);

% join df and df_marks
df_all=outerjoin(df_filtered, df_marks_filtered, ...
    'LeftKeys', 'Ptid', 'RightKeys', 'ptid', 'Type', 'left', 'MergeKeys', true);

% define immune marker
df_all.treatment_actual = categorical(df_all.treatment_actual);
df_all.immune_marker = df_all.Day15pseudoneutid50_BA_1;
ind_arm17 = df_all.treatment_actual == 'Omicron BA.4/5 + Wildtype/Prototype (Pfizer 2)';
df_all.immune_marker(ind_arm17) = df_all.Day15pseudoneutid50_BA_4_BA_5(ind_arm17);


%% define indicators and parameters for specific analyses

% Analysis 1 - Omicron Group (arms 2, 4-6, 8-9, 12, 16, 17)
% Analysis 1.1 Pooling over naive and non-naive
% Analysis 1.3 Restricting to naive==1

% Analysis 2 - Prototype Group (arms 1 and 7):
% Analysis 2.1 Pooling over naive and non-naive
% Analysis 2.3 restricting to naive==1

analyses = {
    struct('a_name','a1_1', 'subset_ind', df_all.TrtB==0,...
    'use_naive', true,  'V_list', {'zhd_spike', 'zhd_rbd', 'hybrid_vx_rbd_rescaled'})
    struct('a_name','a1_3', 'subset_ind', df_all.TrtB==0 & df_all.naive==1,...
    'use_naive', false, 'V_list', {'zhd_spike', 'zhd_rbd', 'hybrid_vx_rbd_rescaled'})
    struct('a_name','a2_1', 'subset_ind', df_all.TrtB==1,...
    'use_naive', true,  'V_list', {'zhd_spike'})
    struct('a_name','a2_3', 'subset_ind', df_all.TrtB==1 & df_all.naive==1,...
    'use_naive', false, 'V_list', {'zhd_spike'})
    };

%% main code for analysis
for iA = 1:length(analyses)
    % extract analysis settings
    a = analyses{iA};
    a_name = a.a_name;
    subset_ind = a.subset_ind;
    use_naive = a.use_naive;
    V_list = {a.V_list};
    
    % subset
    df1 = df_all(subset_ind, :);
    
    % X: minimum of time to COVID-19 or right-censoring
    X = df1.EventTimePrimaryD15/365.25;
    
    % Delta: indicator of the minimum of time to COVID-19 or right-censoring is COVID-19 (only include cases with viral load):
    delta = df1.EventIndPrimaryD15;
    
    % Phase 1 covariates
    standardized_risk_score = df1.standardized_risk_score;
    FOIstandardized = df1.FOIstandardized;
    
    % Phase 2 covariates - immune marker (IM)
    %Note that the immune markers of interest were measured in 100% of trial participants, such that the two-phase sampling
    %part of Sun et al. (2020) is not needed (i.e., all phase-two weights have value 1.0).
    sel = df1.TwophasesampIndD15;
    wei = df1.TwophasesampIndD15.*df1.wt_D15;
    
    % Parameters
    n=size(X,1);
    strata=[zeros(n,1)];
    K=@(x)0.75*(1-x.^2).*(abs(x)<=1);
    
    % covariates
    if use_naive
        Zc=[standardized_risk_score FOIstandardized df1.naive];
    else
        Zc=[standardized_risk_score FOIstandardized];
    end
    Zm=df1.immune_marker;
    % prctile(Zm,[0 10 50 90])
    mu_Zm=nanmean(Zm);
    sigma_Zm=nanstd(Zm);
    Z=[Zm Zc];
    Z_c=[Zc delta Zc.*delta];
    beta_n=size(Z,2);
    
    %For hotdeck multiple imputation, the only variable to use is the calendar time of the COVID-19 endpoint.
    FirstDate = datetime('3/30/22', 'InputFormat', 'MM/dd/yy');
    df1.Av = days(df1.infect_date1 - FirstDate);
    
    % V: Genetic distance marks of interest
    for iV=1:length(V_list)
        V = df1.(V_list{iV});
        
        fprintf('Running %s, variable %s (%d of %d)\n', ...
            a_name, V_list{iV}, iA, length(analyses));
        
        %implement hotdeck imputation
        rep_n = 10;
        L_n=5;
        V_imp = repmat(V, 1, rep_n);
        
        % observed marks
        temp = [V delta df1.Av];
        temp(isnan(V), :)=[];
        n_temp = size(temp, 1);
        
        
        % index for cases with missing marks
        ind_all = [1:n]';
        ind_mis = ind_all(isnan(V) & delta==1);
        n_mis = length(ind_mis);
        
        % impute mark values from bootstrap samples of observed marks
        for i_mis = 1:n_mis
            Vx_mis = df1.Av(ind_mis(i_mis));
            
            ind_B = randsample(n_temp,n_temp,true);
            temp_B = temp(ind_B, :);
            V_obs = temp_B(:, 1);
            Vx_obs = temp_B(:, 3);
            
            ID=knnsearch(Vx_obs,Vx_mis,'K',L_n);
            V_imp(ind_mis(i_mis),:) = V_obs(randsample(ID,rep_n,true));
        end
        %check = [V delta V_imp];
        Vpred=V_imp;
        
        
        %transformation using normalcdf
        mu_V = nanmean(V);
        sigma_V = nanstd(V);
        V(~isnan(V)) = normcdf(V(~isnan(V)), mu_V, sigma_V);
        for kk=1:rep_n
            temp = Vpred(:,kk);
            temp(~isnan(temp)) = normcdf(temp(~isnan(temp)), mu_V, sigma_V);
            Vpred(:,kk) = temp;
        end
        
        % bandwidth for V
        temp = max(diff(unique(V(~isnan(V)))));
        h = 8*nanstd(V) * sum(delta)^(-1/3);
        
        % grid points for V
        n0=25;
        v=linspace(min(V),max(V),n0);  % v: transformed scale
        
        % estimation
        beta_ini=zeros(beta_n,size(v,2));%initial value
        % define limits for integrals
        a=min(V);
        b=max(V);
        aa=a+(v(2)-v(1));
        % define quantiles for CIF
        qtl=[quantile(Z(:,1),0.1),quantile(Z(:,1),0.5), quantile(Z(:,1),0.9)];  % quantiles to be estimate for Z1
        
        if use_naive
            for inaive=0:1
                
                % 	tt: the value of t used for estimating CIF rate
                % 	zpd2: the value of Z1 used for estimating CIF rate
                tt=0.25; % 3 months
                zpd2=[median(standardized_risk_score), median(FOIstandardized), inaive];
                
                tic
                [eff_beta_hat,sig_eff_beta,p_values_t1,p_values_t2, F_0, F_50, F_90]=...
                    aipw_imp1_strt_rubin_multiZc(beta_ini,X,Z,Zc,Zm,Vpred,v,delta,sel,...
                    beta_n,wei,h,n0,n,K,Z_c,rep_n,a,b,aa,strata,tt,zpd2);
                toc
                
                save(sprintf("res/%s_V%d_L5M10_naive%d.mat", a_name, iV, inaive),...
                    "eff_beta_hat","sig_eff_beta",...
                    "p_values_t1","p_values_t2","F_0","F_50","F_90","h","v","a","b","aa",...
                    "mu_V", "sigma_V", "mu_Zm", "sigma_Zm", "qtl")
            end
        else
            % 	tt: the value of t used for estimating CIF rate
            % 	zpd2: the value of Z1 used for estimating CIF rate
            tt=0.25; % 3 months
            zpd2=[median(standardized_risk_score), median(FOIstandardized)];
            
            tic
            [eff_beta_hat,sig_eff_beta,p_values_t1,p_values_t2, F_0, F_50, F_90]=...
                aipw_imp1_strt_rubin_multiZc(beta_ini,X,Z,Zc,Zm,Vpred,v,delta,sel,...
                beta_n,wei,h,n0,n,K,Z_c,rep_n,a,b,aa,strata,tt,zpd2);
            toc
            
            save(sprintf("res/%s_V%d_L5M10.mat", a_name, iV),...
                "eff_beta_hat","sig_eff_beta",...
                "p_values_t1","p_values_t2","F_0","F_50","F_90","h","v","a","b","aa",...
                "mu_V", "sigma_V", "mu_Zm", "sigma_Zm", "qtl")
        end
    end
end

