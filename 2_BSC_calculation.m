function BSC = computeBSC(FCvec, age, UseEmpiricalRegCoef, diag, networkname)
% computeBSC computes the Brain Sex Continuum (BSC) score
% from functional connectivity (FC) data with optional age regression.
%
% INPUTS:
%   FCvec               - Functional connectivity matrix in vectorized form 
%                         (subjects × number of unique FCs).
%   age                 - Age of each subject (vector).
%   UseEmpiricalRegCoef - Flag: 0 = estimate regression parameters from data,
%                         1 = use empirical regression coefficients.
%   diag                - Diagnosis labels for subjects (used to select HC for regression).
%   networkname         - Network name used when computing split BSC.
%
% OUTPUT:
%   BSC                 - Brain Signal Complexity score for each subject.
%
% NOTES:
%   - Supports only AAL2 (94 ROIs) and Power264 (264 ROIs).
%   - Uses pre-trained SVM coefficients (aal2_SVMcoef.mat, power264_SVMcoef.mat).
%   - Can regress out age effects either empirically (UKB/HCP) or via data-driven GLM.

%% Basic information
[n_subject, ~] = size(FCvec);  % number of subjects and FC features
n_ROI = 94;                    % number of ROIs (fixed to 94 for this implementation)

%% Load regression and SVM coefficients
if n_subject <= 10 && UseEmpiricalRegCoef == 0
    error('Too few subjects for reliable data-driven regression.');
end

if n_ROI == 94
    load('aal2_SVMcoef.mat')   % load pre-trained SVM coefficients for AAL2 atlas
    if UseEmpiricalRegCoef == 1
        load('aal2_regcoef.mat')  % load empirical regression coefficients for AAL2
    end
elseif n_ROI == 264
    load('power264_SVMcoef.mat') % load pre-trained SVM coefficients for Power264 atlas
    if UseEmpiricalRegCoef == 1
        load('power264_regcoef.mat') % load empirical regression coefficients for Power264
    end
else
    error('Only AAL2 (94 ROIs) and Power264 (264 ROIs) are currently supported.');
end

%% Load network assignment for each FC
load("network_indexes_for_each_FC.mat"); 
% Network_index variable assigns each FC to a functional network

%% Age regression step
% Ensure age is a column vector
if size(age,2) ~= 1
    age = age';
end

% Polynomial expansion: linear, quadratic, cubic terms of age
age_term = [age, age.^2, age.^3];

if UseEmpiricalRegCoef == 0
    % ---------- Data-driven regression ----------
    fprintf("Using data-driven parameters...\n")
    % Select only healthy controls (HC) to estimate regression coefficients
    selectedRows = strcmp(diag, 'HC');
    % selectedRows = diag==0; % alternative if diag is coded numerically
    
    % For each FC edge, fit GLM with age terms and regress out age effects
    for i_ROI = 1:(n_ROI*(n_ROI-1)/2)
        glmstruct = fitglm(age_term(selectedRows,:), FCvec(selectedRows,i_ROI));
        b = glmstruct.Coefficients.Estimate; % regression coefficients
        FC_regressed(:,i_ROI) = FCvec(:,i_ROI) - [ones(n_subject,1), age_term]*b;
    end
    % Optionally: save adjusted FC
    % writematrix(FC_regressed, 'D:\BSC and MDD\data\XY_FC_age_adjusted.csv');
    
else
    % ---------- Empirical regression ----------
    fprintf("Using empirical parameters...\n")
    % Apply empirical regression coefficients separately for older (>=45) and younger (<45) subjects
    FC_regressed(age >= 45,:) = FCvec(age >= 45,:) - [ones(sum(age >= 45),1), age_term(age >= 45,:)]*b_UKB;
    FC_regressed(age < 45,:)  = FCvec(age < 45,:)  - [ones(sum(age < 45),1), age_term(age < 45,:)]*b_HCP;
end

%% Compute BSC
% Add bias term (column of ones) and normalize FC values
network_FC_data = [ones(n_subject,1) FC_regressed/12];

% Compute classifier score using network-specific SVM coefficients
clfscore = network_FC_data(:,strcmp(Network_index, networkname)) * ...
           SVMbeta(strcmp(Network_index, networkname));

% Apply normal cumulative distribution function (probit transform) 
% to convert classifier score to probability-like BSC score
BSC = normcdf(clfscore);

end
