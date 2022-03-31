%% Thresholds our actual data using a null distribution
function [p_BON,p_FDR] = Threshold_Data_simpler(t,t_null,alpha,tail)

    n_atoms = size(t,1);
    n_regressors = size(t,2);
    n_null = size(t_null,3);

    % For each region and each regressor, computes the thresholds
    t_threshold = prctile(t_null,[alpha/2,100-(alpha/2)],3);

    % Computes the p-values
    for r = 1:n_atoms
        for reg = 1:n_regressors
            switch tail
                case 'both'
                    pval(r,reg) = min([2/n_null*min([sum(squeeze(t_null(r,reg,:)) > t(r,reg)),sum(squeeze(t_null(r,reg,:)) < t(r,reg))]),1]);
                case 'up'
                    pval(r,reg) = sum(squeeze(t_null(r,reg,:)) > t(r,reg))/n_null;
                case 'down'
                    pval(r,reg) = sum(squeeze(t_null(r,reg,:)) < t(r,reg))/n_null;
                % By default, two-tailed
                otherwise
                    pval(r,reg) = min([2/n_null*min([sum(squeeze(abs(t_null(r,reg,:))) > abs(t(r,reg))),sum(squeeze(t_null(r,reg,:)) < t(r,reg))]),1]);
            end
        end
    end
    
    % At this stage, p-values are uncorrected. I want to consider different
    % correction approaches
    
    % Bonferroni correction
    p_BON = n_atoms*pval;
    
    % FDR correction
    q = alpha;
    FDR_method = 'pdep';
    
    [~, ~, ~, p_FDR] = fdr_bh(pval,q,FDR_method);
    
end