%% Data has size n_regions x n_rho
% Null has size n:regions x n_rho x n_null
function [p_BON,p_FDR] = Assess_AUC(Data,Null,tail,Alpha,Rho_range,CM)

    n_regions = size(Data,1);
    n_null = size(Null,3);
    
    % Computes the AUC for each region
    AUC_Actual = sum(Data,2);
    
    % Computes the AUC for the null data; becomes size n_regions x n_null
    AUC_Null = squeeze(sum(Null,2));
    
    % Computes the p-values
    for r = 1:n_regions
        switch tail
            case 'both'
                pval(r) = min([2/n_null*min([sum(squeeze(AUC_Null(r,:)) > AUC_Actual(r)),sum(squeeze(AUC_Null(r,:)) < AUC_Actual(r))]),1]);
            case 'up'
                pval(r) = sum(squeeze(AUC_Null(r,:)) > AUC_Actual(r))/n_null;
            case 'down'
                pval(r) = sum(squeeze(AUC_Null(r,:)) < AUC_Actual(r))/n_null;
            % By default, two-tailed
            otherwise
                pval(r) = min([2/n_null*min([sum(squeeze(AUC_Null(r,:)) > AUC_Actual(r)),sum(squeeze(AUC_Null(r,:)) < AUC_Actual(r))]),1]);
        end
    end
    
    if n_regions == 1
        pval = min([2/n_null*min([sum(squeeze(AUC_Null) > AUC_Actual),sum(squeeze(AUC_Null) < AUC_Actual)]),1]);
    end
    
    % Here, we compute the p-values for the AUC statistics, separately for
    % each region
    
    % Bonferroni correction
    p_BON = n_regions*pval;
    
    % FDR correction
    q = Alpha;
    FDR_method = 'pdep';
    
    [~, ~, ~, p_FDR] = fdr_bh(pval,q,FDR_method);

    % We define for which regions there is at least one significant outcome
    is_sign = (p_BON < Alpha) | (p_FDR < Alpha);
    % is_sign = logical(ones(1,68));
    
    idx_AUC = find(is_sign == 1);
    
    % Plots the overall data results
    figure;
    set(gca,'Box','off');
    boxplot(AUC_Null','plotstyle','compact','colors',CM,'outliersize',4,'symbol','.');
    hold on
    set(gca,'Box','off');
    for r = 1:n_regions
        plot(r,AUC_Actual(r),'Marker','square','MarkerEdgeColor','k',...
            'MarkerFaceColor','k');
    end
    
    CM = cbrewer('seq','PuBuGn',100);
    
    % Plots individual significant cases more in details
    for s = 1:sum(is_sign)

        % Gathers the necessary data
        
        % Size n_rho x n_null
        tmp_plot_null = squeeze(Null(idx_AUC(s),:,:));
        
        % Size n_rho
        tmp_plot = squeeze(Data(idx_AUC(s),:));
        
        figure;
        set(gca,'Box','off');
        
        hold on
        
        for i = 1:68
            plot(Rho_range(1:5),prctile(tmp_plot_null,100-2.5/n_regions*i,2)','Marker','none',...
                'LineWidth',0.25,'Color',CM(end-i,:));
            plot(Rho_range(1:5),prctile(tmp_plot_null,2.5/n_regions*i,2)','Marker','none',...
                'LineWidth',0.25,'Color',CM(end-i,:));
        end
        
        plot(Rho_range(1:5),tmp_plot,'LineStyle','none','Marker','square',...
            'MarkerEdgeColor','none','MarkerFaceColor','k');
    end
end