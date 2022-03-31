%% We want to verify, for given connections, whether they are stronger in
% the expected subject subset case
%
% Inputs:
% - X1 is the data for the HC group, X2 for the Baseline group
% - N is the number of total baseline recordings
% - Rho is the edge density at play
% - n_null is the number of null iterations to perform
% - idx_OI contains the indices of the connections to examine
% - idx_rem_TSTH is the index vector with which subjects to remove first
% - n_rem is the maximum number of subjects to remove
function [Delta,Delta_trimmed_good,Delta_trimmed_bad] = ET_SCM_Compute_EdgeStuff2_HB(X1,X2,Rho,idx_OI,idx_rem_TSTH,n_rem)

    % Actual delta
    Delta = jUpperTriMatToVec(Trim_SCM(Generate_SCM(X1, 100, 'Both','Pearson')-Generate_SCM(X2, 100, 'Both','Pearson'),Rho));
    
    % Selecting only the connections of interest
    Delta = Delta(idx_OI);

    % We will gradually remove up to n_rem subjects
    for rem = 1:n_rem

        X2_trim = X2;
        X2_trim(:,idx_rem_TSTH(1:rem)) = [];
        
        X2b_trim = X2;
        X2b_trim(:,idx_rem_TSTH(end:end-rem+1)) = [];

        tmp1 = jUpperTriMatToVec(Trim_SCM(Generate_SCM(X1, 100, 'Both','Pearson')-Generate_SCM(X2_trim, 100, 'Both','Pearson'),Rho));
        Delta_trimmed_good(:,rem) = tmp1(idx_OI);
        
        tmp2 = jUpperTriMatToVec(Trim_SCM(Generate_SCM(X1, 100, 'Both','Pearson')-Generate_SCM(X2b_trim, 100, 'Both','Pearson'),Rho));
        Delta_trimmed_bad(:,rem) = tmp2(idx_OI);
    end
end