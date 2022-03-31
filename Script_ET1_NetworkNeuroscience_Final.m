%% This script performs graph theoretical analyses on drug-resistant 
% essential tremor patients before and 1 year after stereotactic 
% radiosurgical thalamotomy of the ventro-intermediate nucleus of the 
% thalamus. It also includes the comparison of these patients 
% (pre-intervention)to matched healthy controls. 
%
% The following types of analysis are covered:
% 1. Comparisons (HC vs baseline ET, baseline vs 1-year ET) in terms of 
% structural covariance edges (i.e., covariance across subjects between 
% pairs of areas), at the distribution and individual edge levels. This can
% be run for (a) cortical thickness, (b) surface area or (c) mean curvature
% cortical morphometric features
% 2. Analysis of covariance relationships between pairs of modalities for
% each of the regions at hand, for the three possible pairs of morphometric
% featues (cortical thickness/surface area, surface area/mean curvature,
% cortical thickness/mean curvature)
%
% The present analyses are included in the research article "Morphometric 
% features of drug-resistant essential tremor and recovery after 
% stereotactic radiosurgical thalamotomy", currently Under Review in
% Network Neuroscience.
%
% Written by Thomas A. W. Bolton (Department of Radiology, CHUV, Lausanne)
% and last checked on March 22nd 2022


%% Step 1. Loading of the data

% We add the paths to the data and utilities to use for analyses
addpath(genpath('Data'));
addpath(genpath('Utilities'));

% Loads the data from all 3 groups, as well required covariates
load ET_Clinical
load ET_Freesurfer
load ET_Freesurfer_SC
load Age_HC
load Gender_HC
load ET_Covariates
load CodeBook

% Specifies the type of metrics to probe; here, we will consider cortical
% thickness, surface area and mean curvature
Metric = 'Thickness';
Metric2 = 'Area';
Metric3 = 'MeanCurvature';

% The metrics are loaded, with also subcortical and cerebellar volume
% concatenated to the grey matter regional data (for the subcortex,
% Freesurfer only provides estimates of volume, and no further
% decomposition into area, thickness or the likes)

% Cortical thickness
X_HC = [(eval(['HC_',Metric]))';HC_Volume_SC'];
X_BASE = [(eval(['BASE_',Metric]))';BASE_Volume_SC'];
X_YEAR = [(eval(['YEAR_',Metric]))';YEAR_Volume_SC'];

% Surface area
Y_HC = [(eval(['HC_',Metric2]))';HC_Volume_SC'];
Y_BASE = [(eval(['BASE_',Metric2]))';BASE_Volume_SC'];
Y_YEAR = [(eval(['YEAR_',Metric2]))';YEAR_Volume_SC'];

% Mean curvature
Z_HC = [(eval(['HC_',Metric3]))';HC_Volume_SC'];
Z_BASE = [(eval(['BASE_',Metric3]))';BASE_Volume_SC'];
Z_YEAR = [(eval(['YEAR_',Metric3]))';YEAR_Volume_SC'];

% Number of subjects in total
n_HC = size(X_HC,2);
n_BASE = size(X_BASE,2);
n_YEAR = size(X_YEAR,2);

% Number of regions in the parcellation at play
n_regions = size(X_HC,1);
n_regions_cortical = 68;
n_regions_subcortical = 19;

%%%%%%%%%% IN ORDER TO DISPLAY BRAIN PLOTS, ADD FREESURFER TO YOUR PATH



%% Step 2. Definition of analytical parameters and initialization

%%%%%% PARAMETERS

% Will we plot the information or not?
is_plot = 1;

% Will we compute a non-parametric null distribution or not?
is_null = 1;

% Type of graph to generate: choose between Positive (only retains positive
% SCM edges), Negative (only retains and sign flips the negative edges), or
% Both (keeps both edge types with their signs). Here, we want to conduct
% edge-wise assessments, so we want to consider both signs
Type = 'Both';

% Type of covariance to compute; here, simply Pearson's correlation
% coefficient
Corr_type = 'Pearson';

% Number of edges in total with our number of regions, prior to exclusion
n_edges = n_regions*(n_regions-1)/2;

%%%%% Optimal density for the structural covariance matrices (in terms of
%%%%% coefficient of determination R^2). This value was selected based on
%%%%% past literature, so that it removes "weak correlation" edges and
%%%%% keeps "strong correlation" edges
Threshold_density = 0.2;

% Covariates to regress out for "traditional" analyses: note that we demean
% the data for each group in addition to the rest. Age, gender and total
% grey matter volume are included
Cov_BASE = [ones(n_BASE,1),Age,Gender,BASE_TGV];
Cov_YEAR = [ones(n_YEAR,1),Age+1,Gender,YEAR_TGV];
Cov_HC = [ones(n_HC,1),Age_HC,Gender_HC,HC_TGV];

% Number of null realisations to generate for edge-wise analyses. We
% selected this number so that it is larger than 38*n_edges (the number
% of samples required for Bonferroni correction of results at the edge-wise
% level)
n_null_edges = 300000;

% Number of null realisations to generate for cross-modality analyses. We
% selected this number so that it is larger than 38*n_regions*3 (the number
% of samples required for Bonferroni correction, considering that we
% consider three cross-modality relationships)
n_null_regions = 50000;

% Significance level at which we want to test in the
% analyses, and type of tail for the tests ('both' = two-tailed). We have
% no prior assumptions about the direction of potential differences, so we
% will use two-tailed testing
Alpha = 0.05;
tail = 'both';

% Initial density at which to construct structural covariance matrices;
% here, we want to keep all edges initially (and remove noisy ones later
% on)
Rho = 100;

% How many subjects do we want to remove at most in our associations with
% clinical symptoms analyses?
n_rem = 10;

% Mapping to lobes (including subcortical areas and cerebellum, regarded as
% lobes of their own). Region 87 is the brainstem. Mapping follows the
% suggestion from Klein and Tourville, Frontiers in Neuroscience 2012, "101
% labeled brain images and a consistent human cortical labeling protocol"
LobesMapping = [[5,9,1,7,5,5,3,5,9,7,1,7,1,5,5,1,1,1,1,7,3,9,1,3,9,1,...
    1,3,5,3,1,5,5,11],1+[5,9,1,7,5,5,3,5,9,7,1,7,1,5,5,1,1,1,1,7,3,9,1,...
    3,9,1,1,3,5,3,1,5,5,11]]';
LobesMapping(69:77) = 13;
LobesMapping(78:86) = 14;
LobesMapping(87) = 15;


%%%%%%%% COLORMAPS

% The following colormaps are generated thanks to the cbrewer MATLAB
% toolbox, and used for various plotting purposes throughout the script
CM_RB = cbrewer('div','RdBu',1001);
CM_RB(CM_RB<0) = 0;

CM_RYG = cbrewer('div','RdYlGn',1001);
CM_RYG(CM_RYG < 0) = 0;

CM_Paired = cbrewer('qual','Paired',12);
CM_Paired(CM_Paired < 0) = 0;

CM_YG = cbrewer('seq','YlGn',1001);
CM_YG(CM_YG < 0) = 0;

CM_Subjects = colorcube(n_BASE);
CM_Subjects2 = flipud(cbrewer('div','RdYlGn',n_BASE));
CM_Subjects2(CM_Subjects2 < 0) = 0;

CM_RYB_thresh = cbrewer('div','RdYlBu',1000);
CM_RYB_thresh(CM_RYB_thresh<0) = 0;
CM_RYB_thresh = [CM_RYB_thresh(1:500,:);[0.6,0.6,0.6];...
    CM_RYB_thresh(501:end,:)];

CM_RB_thresh = cbrewer('div','RdBu',1000);
CM_RB_thresh(CM_RB_thresh<0) = 0;
CM_RB_thresh = [CM_RB_thresh(1:500,:);[0.6,0.6,0.6];...
    CM_RB_thresh(501:end,:)];



%% Step 3. Regressing out covariates from the data
% For the regression, the "y_regress_ss" function is used

% Cortical thickness
for r = 1:n_regions
    [~,X_BASE_res(r,:)] = y_regress_ss(X_BASE(r,:)',Cov_BASE);
    [~,X_YEAR_res(r,:)] = y_regress_ss(X_YEAR(r,:)',Cov_YEAR);
    [~,X_HC_res(r,:)] = y_regress_ss(X_HC(r,:)',Cov_HC);
end

% Surface area
for r = 1:n_regions
    [~,Y_BASE_res(r,:)] = y_regress_ss(Y_BASE(r,:)',Cov_BASE);
    [~,Y_YEAR_res(r,:)] = y_regress_ss(Y_YEAR(r,:)',Cov_YEAR);
    [~,Y_HC_res(r,:)] = y_regress_ss(Y_HC(r,:)',Cov_HC);
end

% Mean curvature
for r = 1:n_regions
    [~,Z_BASE_res(r,:)] = y_regress_ss(Z_BASE(r,:)',Cov_BASE);
    [~,Z_YEAR_res(r,:)] = y_regress_ss(Z_YEAR(r,:)',Cov_YEAR);
    [~,Z_HC_res(r,:)] = y_regress_ss(Z_HC(r,:)',Cov_HC);
end

% Note that the data for cortical thickness is prefixed with the letter X,
% that for surface area with Y, and for mean curvature with Z. In
% subsequent edge-wise analyses, the user wants to adjust the inputs
% according to the case at hand



%% Step 4. Generation of the Structural Community Matrices
% In what follows, replace "X" by "Y" or "Z" to consider modalities 2 or 3
% instead of 1

% We generate the structural covariance matrices for each of the three
% groups, at 100% density, with Pearson's correlation coefficient
SCM_HC = Generate_SCM(X_HC_res, Rho, Type,Corr_type);
SCM_BASE = Generate_SCM(X_BASE_res, Rho, Type,Corr_type);
SCM_YEAR = Generate_SCM(X_YEAR_res, Rho, Type,Corr_type);

% Sets diagonal elements to zero
for i = 1:n_regions
    SCM_HC(i,i) = 0;
    SCM_BASE(i,i) = 0;
    SCM_YEAR(i,i) = 0;
end

% Coefficients of determination (R^2) vectorised for each group
S_YEAR = jUpperTriMatToVec(SCM_YEAR.^2);
S_BASE = jUpperTriMatToVec(SCM_BASE.^2);
S_HC = jUpperTriMatToVec(SCM_HC.^2);

% We remove the edges for which there is never a strong correlation in any
% of the three considered groups, and store this information in the
% logical vector "idx_toremove"
idx_toremove = (S_HC<Threshold_density & S_BASE<Threshold_density & ...
    S_YEAR<Threshold_density);
disp(['There are ',num2str(sum(idx_toremove)),' out of ',...
    num2str(n_edges),' edges removed!']);

% The matrices for each group are vectorised, and trimmed according to our
% selected cutoff
VCM_HC = jUpperTriMatToVec(SCM_HC);
VCM_HC(idx_toremove) = 0;
VCM_BASE = jUpperTriMatToVec(SCM_BASE);
VCM_BASE(idx_toremove) = 0;
VCM_YEAR = jUpperTriMatToVec(SCM_YEAR);
VCM_YEAR(idx_toremove) = 0;

% Then, they are sent back to the matrix representational space
SCM_HC = jVecToUpperTriMat(VCM_HC,n_regions);
SCM_BASE = jVecToUpperTriMat(VCM_BASE,n_regions);
SCM_YEAR = jVecToUpperTriMat(VCM_YEAR,n_regions);

% Differences between groups, for both contrasts of interest
Delta_SCM_YB = SCM_YEAR - SCM_BASE;
Delta_SCM_HB = SCM_HC - SCM_BASE;


% Plotting of the resulting outputs
if is_plot
    
    % Individual structural covariance matrices for the three groups
    % (displayed in Figure 1A, B and C, top row)
    figure;
    imagesc(SCM_HC);
    colormap(flipud(CM_RYB_thresh));
    caxis([-1,1]);
    axis off
    
    figure;
    imagesc(SCM_BASE);
    colormap(flipud(CM_RYB_thresh));
    caxis([-1,1]);
    axis off
    
    figure;
    imagesc(SCM_YEAR);
    colormap(flipud(CM_RYB_thresh));
    caxis([-1,1]);
    axis off
    
    % Structural covariance contrasts for the two differences of interest
    % (displayed in Figure 1A, B and C, bottom row, left and middle)
    figure;
    imagesc(Delta_SCM_YB);
    colormap(flipud(CM_RB_thresh));
    caxis([-1.6,1.6]);
    
    figure;
    imagesc(Delta_SCM_HB);
    colormap(flipud(CM_RB_thresh));
    caxis([-1.6,1.6]);
       
    % Summarizing histograms of structural covariance values for the three
    % groups (displayed in Figure 1A, B and C, bottom row, right). Note
    % that we have to set the "zero" elements to NaN so that we avoid an
    % ugly artificial peak at 0
    figure;
    hold on
    SCM_HC(SCM_HC==0) = NaN;
    h = histogram(jUpperTriMatToVec(SCM_HC),100);
    set(h,'FaceColor',[70,148,73]/255,'EdgeColor','None','FaceAlpha',0.8);
    
    SCM_BASE(SCM_BASE==0) = NaN;
    h2 = histogram(jUpperTriMatToVec(SCM_BASE),100);
    set(h2,'FaceColor',[56,61,150]/255,'EdgeColor','None','FaceAlpha',0.2);
    
    SCM_YEAR(SCM_YEAR==0) = NaN;
    h3 = histogram(jUpperTriMatToVec(SCM_YEAR),100);
    set(h3,'FaceColor',[175,54,60]/255,'EdgeColor','None','FaceAlpha',0.2);
    
    xlim([-1,1]);
    ylim([0,100]);
end



%% Step 5. Comparison across groups for each edge ("classical" approach)

% Differences of interest in vectorised format
Velta_SCM_YB = jUpperTriMatToVec(Delta_SCM_YB);
Velta_SCM_HB = jUpperTriMatToVec(Delta_SCM_HB);

% Will contain the null differences used for significance testing
Delta_SCM_null_YB = NaN(n_edges,1,n_null_edges);
Delta_SCM_null_HB = NaN(n_edges,1,n_null_edges);

% For loop for null data generation
for n = 1:n_null_edges
    
    n

    % In each of the null realization cases, we want to shuffle the
    % subjects, recompute structural covariance matrices, and store the
    % information

    % First, we shuffle the subjects for the "base - year" contrast; we put
    % all the available data in tmp_X, of size 2*n_BASE x n_regions. We
    % create a vector with random indices, and use it for sampling both
    % "null groups". Each group has size n_BASE
    tmp_X = [X_BASE_res';X_YEAR_res'];
    idx = randperm(2*n_BASE);
    tmp_X = tmp_X(idx,:);
    tmp_BASE = tmp_X(1:n_BASE,:);
    tmp_YEAR = tmp_X(n_BASE+1:end,:);

    % We perform the same process for the HC - ET pre-intervention
    % contrast, where this time we have a total of n_HC + n_BASE subjects
    tmp_X = [X_BASE_res';X_HC_res'];
    idx = randperm(n_HC + n_BASE);
    tmp_X = tmp_X(idx,:);
    tmp_BASE2 = tmp_X(1:n_BASE,:);
    tmp_HC = tmp_X(n_BASE+1:end,:);

    % Generates null data for each group and subsequent difference. We
    % compute it for all edges, but of course we will only analyse the ones
    % retained above
    Delta_SCM_null_YB(:,1,n) = jUpperTriMatToVec(Generate_SCM(tmp_YEAR',...
        Rho,Type,Corr_type)-Generate_SCM(tmp_BASE',Rho,Type,Corr_type));
    Delta_SCM_null_HB(:,1,n) = jUpperTriMatToVec(Generate_SCM(tmp_HC',...
        Rho,Type,Corr_type)-Generate_SCM(tmp_BASE2', Rho, Type,Corr_type));
end

% Only for the edges that we want to retain (~idx_toremove), we compute
% FDR-corrected p-values
[p_BON_edges_YB,p_FDR_edges_YB] = ...
    Threshold_Data_simpler(Velta_SCM_YB(~idx_toremove),...
    Delta_SCM_null_YB(~idx_toremove,:,:),Alpha,tail);

[p_BON_edges_HB,p_FDR_edges_HB] = ...
    Threshold_Data_simpler(Velta_SCM_HB(~idx_toremove),...
    Delta_SCM_null_HB(~idx_toremove,:,:),Alpha,tail);

% Computation of the indices of the significant edges: we start with all,
% then we remove the noisy ones, then we select the significant ones
indices = 1:n_edges;
indices(idx_toremove) = [];
idx_sign_edges_YB = indices(p_FDR_edges_YB<Alpha);
idx_sign_edges_HB = indices(p_FDR_edges_HB<Alpha);



%% Step 6. Associations with the extent of clinical symptoms
% We want to perform additional validations on the extracted significant
% connections. 
% For the HC - Baseline ET comparison, we expect that connections
% associated to the motor symptoms will evolve more markedly in structural
% covariance when strongly impaired patients are removed, as opposed to
% when mildly impaired ones are removed (i.e., the most impaired patients
% largely drive the results).
% For the 1-year - baseline ET comparison, we expect that connections 
% associated to the extent of motor recovery will vary more if the "best
% recoverers" are not included, as compared to not including the "worst"
% recoverers.

% We sort the subjects according to their baseline tremor extent (TSTH, 
% Tremor Score on Treated Hand). Subjects with a
% larger score tremble more, thus the first indices of the vector will
% reflect the least impaired subjects tremor-wise
[~,idx_rem_TSTH] = sort(TSTH_Baseline,'ascend');

% Same process, but sorting in terms of TSTH recovery. First indices
% reflect the subjects who improved the least
[~,idx_rem_TSTHPerc] = sort(TSTH_PercAllev,'ascend');

% This function sees how the connection value evolves when recomputing by
% removing the subjects that were the least impaired as opposed to the
% most.
% The _bad output reflects the case for which the least impaired subjects
% were removed, and the _good output the others
if ~isempty(idx_sign_edges_HB)
    [Delta_HB,Delta_trimmed_good_HB,Delta_trimmed_bad_HB] =...
        ET_SCM_Compute_EdgeStuff2_HB(Y_HC_res,Y_BASE_res,Rho,...
        idx_sign_edges_HB,idx_rem_TSTH,n_rem);
end

% This function sees how the connection value evolves when recomputing by
% removing the subjects that improved the most as opposed to the least
% The _good outputs reflects the evolution of structural covariance when
% the subjects who recover the least are removed, and the _bad one, when
% the subjects who recover the most are removed
if ~isempty(idx_sign_edges_YB)
    [Delta_YB,Delta_trimmed_good_YB,Delta_trimmed_bad_YB] = ...
        ET_SCM_Compute_EdgeStuff2_YB(X_YEAR_res,X_BASE_res,Rho,...
        idx_sign_edges_YB,idx_rem_TSTHPerc,n_rem);
end

% Plots
if is_plot

    % Which connections are significant in a matrix format for both
    % contrasts
    edge_vec_YB = zeros(n_edges,1);
    edge_vec_YB(idx_sign_edges_YB) = 1;
    
    figure;
    set(gca,'Box','off');
    imagesc(jVecToUpperTriMat(edge_vec_YB,n_regions));
    colormap(flipud(cbrewer('div','RdGy',1000)));
    caxis([-2,2]);
    
    edge_vec_HB = zeros(n_edges,1);
    edge_vec_HB(idx_sign_edges_HB) = 1;
    
    figure;
    set(gca,'Box','off');
    imagesc(jVecToUpperTriMat(edge_vec_HB,n_regions));
    colormap(flipud(cbrewer('div','RdGy',1000)));
    caxis([-2,2]);

    % Summary of the additional analyses for each significant connection,
    % for both contrasts
    for conn = 1:length(idx_sign_edges_YB)
        figure;
        title(['Connection ID ',num2str(idx_sign_edges_YB(conn))]);
        hold on;
        set(gca,'Box','off');
        plot(0,Delta_YB(conn),'ks','MarkerSize',10);
        plot(1:n_rem,Delta_trimmed_good_YB(conn,:),'b','LineWidth',2);
        plot(1:n_rem,Delta_trimmed_bad_YB(conn,:),'r','LineWidth',2);
        xlim([-0.5,n_rem]);
    end
    
    for conn = 1:length(idx_sign_edges_HB)
        figure;
        title(['Connection ID ',num2str(idx_sign_edges_HB(conn))]);
        hold on;
        set(gca,'Box','off');
        plot(0,Delta_HB(conn),'bs','MarkerSize',10);
        plot(1:n_rem,Delta_trimmed_good_HB(conn,:),'b','LineWidth',2);
        plot(1:n_rem,Delta_trimmed_bad_HB(conn,:),'r','LineWidth',2);
        xlim([-0.5,n_rem]);
    end

    % Plots linked to each significant connection for both contrasts
    for conn = 1:length(idx_sign_edges_YB)

        % Takes the actual value
        tmp_actual = Velta_SCM_YB(idx_sign_edges_YB(conn));
    
        % Takes the null values
        tmp_null = squeeze(Delta_SCM_null_YB(idx_sign_edges_YB(conn),1,:));

        figure;
        title(['Connection ID ',num2str(idx_sign_edges_YB(conn))]);
        hold on
        h = histogram(tmp_null,1000,'Normalization','pdf');
        set(h,'FaceColor',[94,60,108]/255,'EdgeColor','None',...
            'FaceAlpha',0.8);
        plot([tmp_actual,tmp_actual],[0,3.5],'k--');
    end

    for conn = 1:length(idx_sign_edges_HB)

        % Takes the actual value
        tmp_actual = Velta_SCM_HB(idx_sign_edges_HB(conn));
    
        % Takes the null values
        tmp_null = squeeze(Delta_SCM_null_HB(idx_sign_edges_HB(conn),1,:));

        figure;
        title(['Connection ID ',num2str(idx_sign_edges_HB(conn))]);
        hold on
        h = histogram(tmp_null,1000,'Normalization','pdf');
        set(h,'FaceColor',[94,60,108]/255,'EdgeColor','None',...
            'FaceAlpha',0.8);
        plot([tmp_actual,tmp_actual],[0,3.5],'k--');
    end

    % Also plots the connections in a brain (this will only plot
    % cortico-cortical connections and associated nodal degrees)
%     figure;
%     h = gca;
%     tmp_data = jVecToUpperTriMat(edge_vec_HB,n_regions);
%     tmp_data = tmp_data(1:n_regions_cortical,1:n_regions_cortical);
%     PlotBrainGraph_Kmeans_withtitle(h,tmp_data,sum(tmp_data),CodeBook,0.1,...
%         0.1,0.1,1,1,'hot',...
%         'spring',1,0.8,[],'','');
% 
%     figure;
%     h = gca;
%     tmp_data = jVecToUpperTriMat(edge_vec_YB,n_regions);
%     tmp_data = tmp_data(1:n_regions_cortical,1:n_regions_cortical);
%     PlotBrainGraph_Kmeans_withtitle(h,tmp_data,sum(tmp_data),CodeBook,0.1,...
%         0.1,0.1,1,1,'hot',...
%         'spring',1,0.8,[],'','');
end



%% Step 7. Assessment of cross-property covariance
% For this section, X_ and Y_ should be replaced by whichever prefix
% matches the morphometric properties of interest

% Correlation across modalities is computed: note that this is only
% possible for cortical regions, because we only have volume quantified for
% cerebellum and subcortex. One value is obtained per region
COV_BASE = diag(corr(X_BASE_res(1:n_regions_cortical,:)',...
    Z_BASE_res(1:n_regions_cortical,:)'));
COV_HC = diag(corr(X_HC_res(1:n_regions_cortical,:)',...
    Z_HC_res(1:n_regions_cortical,:)'));
COV_YEAR = diag(corr(X_YEAR_res(1:n_regions_cortical,:)',...
    Z_YEAR_res(1:n_regions_cortical,:)'));

% Differences of interest
COV_YB = COV_YEAR - COV_BASE;
COV_HB = COV_HC - COV_BASE;

% Null data scheme, where we use a similar strategy as for edge-wise
% analyses above
if is_null
    for n = 1:n_null_regions

        n
        
        % The data for both morphometric properties of interest is sampled,
        % and then shuffled across both groups, first for the 1-year -
        % baseline contrast...
        tmp_X = [X_BASE_res(1:n_regions_cortical,:)';...
            X_YEAR_res(1:n_regions_cortical,:)'];
        tmp_Y = [Z_BASE_res(1:n_regions_cortical,:)';...
            Z_YEAR_res(1:n_regions_cortical,:)'];
        idx = randperm(n_BASE+n_YEAR);
        tmp_X = tmp_X(idx,:);
        tmp_Y = tmp_Y(idx,:);
        tmp_BASE_X = tmp_X(1:n_BASE,:);
        tmp_BASE_Y = tmp_Y(1:n_BASE,:);
        tmp_YEAR_X = tmp_X(n_BASE+1:end,:);
        tmp_YEAR_Y = tmp_Y(n_BASE+1:end,:);

        % ... and then for the HC - ET baseline contrast
        tmp_X = [X_BASE_res(1:n_regions_cortical,:)';...
            X_HC_res(1:n_regions_cortical,:)'];
        tmp_Y = [Z_BASE_res(1:n_regions_cortical,:)';...
            Z_HC_res(1:n_regions_cortical,:)'];
        idx = randperm(n_BASE+n_HC);
        tmp_X = tmp_X(idx,:);
        tmp_Y = tmp_Y(idx,:);
        tmp_BASE2_X = tmp_X(1:n_BASE,:);
        tmp_BASE2_Y = tmp_Y(1:n_BASE,:);
        tmp_HC_X = tmp_X(n_BASE+1:end,:);
        tmp_HC_Y = tmp_Y(n_BASE+1:end,:);
        
        % Null covariances are computed
        COV_YB_null(:,n) = diag(corr(tmp_YEAR_X,tmp_YEAR_Y))-...
            diag(corr(tmp_BASE_X,tmp_BASE_Y));
        COV_HB_null(:,n) = diag(corr(tmp_HC_X,tmp_HC_Y))-...
            diag(corr(tmp_BASE2_X,tmp_BASE2_Y));
    end
end

% Computes the p-values
for r = 1:n_regions_cortical
    pval_cov_classical_YB(r) = 2/(n_null_regions)*...
        min([sum(squeeze(COV_YB_null(r,:)) > COV_YB(r)),...
        sum(squeeze(COV_YB_null(r,:)) < COV_YB(r))]);
    pval_cov_classical_HB(r) = 2/(n_null_regions)*...
        min([sum(squeeze(COV_HB_null(r,:)) > COV_HB(r)),...
        sum(squeeze(COV_HB_null(r,:)) < COV_HB(r))]);
end

q = Alpha;
FDR_method = 'pdep';

[~, ~, ~, p_cov_YB_FDR] = fdr_bh(pval_cov_classical_YB,q,FDR_method);
[~, ~, ~, p_cov_HB_FDR] = fdr_bh(pval_cov_classical_HB,q,FDR_method);

if is_plot
    
    % Results are illustrated by sets of boxplots (null distributions) and
    % actual values (dark squares). They are presented in Figure 2
    figure;
    set(gca,'Box','off');
    boxplot(COV_YB_null','plotstyle','compact','colors',...
        CM_Paired(LobesMapping(1:n_regions_cortical),:),...
        'outliersize',4,'symbol','.');
    hold on
    set(gca,'Box','off');
    plot(COV_YB,'Marker','square','MarkerEdgeColor','k',...
        'MarkerFaceColor','k','LineStyle','none');
    
    figure;
    set(gca,'Box','off');
    boxplot(COV_HB_null','plotstyle','compact','colors',...
        CM_Paired(LobesMapping(1:n_regions_cortical),:),...
        'outliersize',4,'symbol','.');
    hold on
    set(gca,'Box','off');
    plot(COV_HB,'Marker','square','MarkerEdgeColor','k',...
        'MarkerFaceColor','k','LineStyle','none');
end