%% This function computes a graph reflective of structural covariance for
% a given measure of interest (e.g., thickness)
%
% Inputs:
% - Data has size n_regions x n_subjects
% - Rho is the density of edges desired (in percent)
%
% Outputs:
% - SCM has size n_regions x n_regions: it represents the pre-treated SCM
% - is_connected defines whether the end matrix is fully connected or not
% at the required density of edges
% - Rho_init is the initial edge density before trimming
% - G is the output matrix at the desired edge density Rho
function [G, Rho_init, is_connected, SCM] = Generate_SCM(Data, Rho, Type, Measure)

    switch Measure
        
        case 'Pearson'
            
            % Computes the SCM using Pearson's correlation coefficient
            SCM = corr(Data');
            
        case 'Partial'
            
            % Computes the SCM using Pearson's correlation coefficient
            SCM = partialcorr(Data');
            
        case 'Covariance'
            
            % Computes the SCM using Pearson's correlation coefficient
            SCM = cov(Data');
            
        otherwise
            
            errordlg('NOOOOOOOOO');
            
    end
            
            
    
    
    % Number of regions
    N = size(SCM,1);
    
    % Maximum number of edges
    N_max = N*(N-1)/2;
    
    % Thresholds to retain only a certain edge density
    SCM_vec = jUpperTriMatToVec(SCM);
    
    % How do we want to construct the SCM: only positive-valued edges, only
    % negative-valued edges, or both (absolute-valued)?
    switch Type
        case 'Positive'
            
            SCM_vec(SCM_vec < 0) = 0;
            
        case 'Negative'
            
            SCM_vec(SCM_vec > 0) = 0;
            SCM_vec = abs(SCM_vec); 
            
        case 'Both'
            
    end
    
    [~,idx] = sort(abs(SCM_vec),'descend');
    
    % Initial edge density before further trimming
    Rho_init = 100*sum(abs(SCM_vec) > 0)/N_max;

    % Required number of remaining edges
    N_final = floor(Rho*N_max/100);
    
    if Rho > Rho_init
        disp('Cannot reach the required density: too few non-null edges!');
    end
    
    SCM_vec(idx(N_final+1:end)) = 0;
    
    G = jVecToUpperTriMat(SCM_vec,N) + jVecToUpperTriMat(SCM_vec,N)';
    %G = G + diag(ones(N,1));
    
    % To assess if a graph is fully connected, we compute the Laplacian and
    % check the second lowest eigenvalue (if not null, then fully
    % connected)
    D = diag(sum(G,2));
    Lambda = eig(D-G);
    Lambda = sort(Lambda,'ascend');

    is_connected = (Lambda(2) > 0);
end