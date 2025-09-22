function BSC = computeBSC(FCvec,age,UseEmpiricalRegCoef,diag,networkname)% only use the networkname when computing split BSC

% basic information
[n_subject,~] = size(FCvec);
n_ROI = 94;

% loading regressing and SVM coefficients
if n_subject<=10 && UseEmpiricalRegCoef == 0
    error('Too few subjects, inappropriate for regression.\n')
end
if n_ROI == 94
    load('aal2_SVMcoef.mat')
    if UseEmpiricalRegCoef == 1
        load('aal2_regcoef.mat')
    end
else if n_ROI == 264
        load('power264_SVMcoef.mat')
        if UseEmpiricalRegCoef == 1
            load('power264_regcoef.mat')
        end
else
    error('Current version only support compute BSC based on AAL2 and Power264.\n')
end
end

% load network names for each FC
load("network_indexes_for_each_FC.mat");

% regressing out age terms
if size(age,2)~=1
    age = age';
end
age_term = [age,age.^2,age.^3];
if UseEmpiricalRegCoef == 0
    fprintf("Using data-driven parameters...\n")
    selectedRows = strcmp(diag, 'HC');
    %selectedRows = diag==0;
    for i_ROI = 1:(n_ROI*(n_ROI-1)/2)
        glmstruct = fitglm(age_term(selectedRows,:),FCvec(selectedRows,i_ROI));
        b = glmstruct.Coefficients.Estimate;
        FC_regressed(:,i_ROI) = FCvec(:,i_ROI)-[ones(n_subject,1),age_term]*b;
    end
    %writematrix(FC_regressed, 'D:\BSC and MDD\data\XY_FC_age_adjusted.csv');

else
    fprintf("Using empirical parameters...\n")
    FC_regressed(age>=45,:) = FCvec(age>=45,:) - [ones(sum(age>=45),1),age_term(age>=45,:)]*b_UKB;
    FC_regressed(age<45,:) = FCvec(age<45,:) - [ones(sum(age<45),1),age_term(age<45,:)]*b_HCP;
end

% compute BSC
network_FC_data = [ones(n_subject,1) FC_regressed/12];
clfscore = network_FC_data(:,strcmp(Network_index, networkname))*SVMbeta(strcmp(Network_index, networkname));
BSC = normcdf(clfscore);
