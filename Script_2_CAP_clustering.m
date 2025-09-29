%% CAP clustering script (seed-free, TbCAPs, no GUI)
% Takes in consensus results from Script_1_CAP_Discovery.m 


%% Load consensus results
%load('/Users/npsabet/TbCAPS/SavedData/Consensus_20250927_1647.mat')
%load('/home/nsabet/TbCAPS/SavedData/Consensus_20250929_full.mat')

whos  % shows all variables (sizes, types, memory use)

% You should see at least: Consensus, K_range, Pcc, N
disp('Consensus size:')
disp(size(Consensus))
disp('K_range:')
disp(K_range)
disp('Number of folds (N):')
disp(N)
disp('Subsampling % (Pcc):')
disp(Pcc)

% Compute clustering quality curve
% [~, Qual] = ComputeClusteringQuality(Consensus,[]);
% figure;
% plot(K_range, Qual,'-o','LineWidth',2)
% xlabel('K'); ylabel('Quality metric')
% title('Consensus clustering quality (Controls-pre)')
% Best K should usually be around 5–8 (look for an elbow/peak in Qual)

%% Clustering into CAPs

% Set parameters
K_opt = 4;      % chosen from consensus quality curve
Pp    = 100;    % percentage of positive voxels to retain
Pn    = 100;    % percentage of negative voxels to retain
n_rep = 100;    % number of k-means repetitions for stability

% Run final clustering
[CAP,~,~,idx] = Run_Clustering(cell2mat(Xon), ...
    K_opt, mask, brain_info, Pp, Pn, n_rep, [], []);

% Save CAP results and all relevant parameters
basis_out = fullfile('/home/nsabet/TbCAPS/SavedData', ...
    sprintf('CAP_basis_RDoC_controls_pre_%s.mat', datestr(now,'yyyymmdd')));
save(basis_out, ...
     'CAP','idx','mask','brain_info', ...    
     'K_opt','Pp','Pn','n_rep','K_range','N','Pcc');
fprintf('Saved CAP basis to %s\n', basis_out);

% Save CAP maps as NIfTI
for k = 1:K_opt
    CAP_map = zeros(brain_info.dim);
    CAP_map(mask) = CAP(k,:);
    Vout = brain_info;
    Vout.fname = fullfile('/home/nsabet/TbCAPS/SavedData', ...
        sprintf('CAP_controls_pre_%s_%d.nii', datestr(now,'yyyymmdd'), k));
    spm_write_vol(Vout, CAP_map);
    fprintf('Saved CAP %d map as NIfTI\n', k);
end

% Save z-scored CAPs
CAP_z = CAP_Zscore(CAP);
z_out = fullfile('/home/nsabet/TbCAPS/SavedData', ...
    sprintf('CAP_z_RDoC_controls_pre_%s.mat', datestr(now,'yyyymmdd')));
save(z_out,'CAP_z','K_opt','mask','brain_info');
fprintf('Saved z-scored CAPs to %s\n', z_out);


%% Compute metrics

TR = 2;  % adjust if needed
[ExpressionMap, Counts, Entries, Avg_Duration, Duration, TransitionProbabilities, ...
    From_Baseline, To_Baseline, Baseline_resilience, Resilience, Betweenness, ...
    InDegree, OutDegree, SubjectEntries] = ...
    Compute_Metrics_simpler(idx, Indices.kept.active, Indices.scrubbedandactive, K_opt, TR);

metrics_out = fullfile('/home/nsabet/TbCAPS/SavedData', ...
    sprintf('CAP_metrics_RDoC_controls_pre_%s.mat', datestr(now,'yyyymmdd')));
save(metrics_out, ...
    'ExpressionMap','Counts','Entries','Avg_Duration','Duration', ...
    'TransitionProbabilities','From_Baseline','To_Baseline', ...
    'Baseline_resilience','Resilience','Betweenness','InDegree', ...
    'OutDegree','SubjectEntries','TR','K_opt');
fprintf('Saved metrics to %s\n', metrics_out);

%% QC Figure
figure('Position',[100 100 1400 900]);

% Temporal barcode
subplot(3,2,[1 2]);
nSubj = numel(Xon);
startIdx = 1;
for s = 1:nSubj
    nFramesSubj = size(Xon{s},2);
    subj_idx = idx(startIdx:startIdx+nFramesSubj-1);
    startIdx = startIdx + nFramesSubj;
    imagesc(1:length(subj_idx), [s-1 s], subj_idx'); hold on;
end
colormap(jet(K_opt));
caxis([1 K_opt]);
ylabel('Subjects'); xlabel('Frames');
title('Temporal barcode (CAP sequence)');

% Dwell time violin plots
subplot(3,2,3);
ah = gca;
dur_data = Avg_Duration(:,1:K_opt)'; 
Lab = arrayfun(@(k) sprintf('CAP%d',k), 1:K_opt, 'UniformOutput', false);
Color = {lines(K_opt)*0.5 + 0.5, lines(K_opt)};
MakeViolin(dur_data, ah, Lab, 'Duration (sec)', Color, 1, K_opt);
title('Dwell time per CAP');

% Transition probability heatmap
subplot(3,2,4);
Tmat = squeeze(mean(TransitionProbabilities(1:K_opt,1:K_opt,:),3));
Tmat = Tmat ./ max(sum(Tmat,2),eps);  
imagesc(Tmat,[0 1]);
colormap(gca,hot); colorbar; axis square;
xlabel('To CAP'); ylabel('From CAP');
title('Row-normalized transition probabilities');

% CAP counts
subplot(3,2,[5 6]);
cap_counts = nansum(Entries(:,1:K_opt),1);
bar(1:K_opt, cap_counts, 'FaceColor',[0.3 0.5 0.8]);
xticks(1:K_opt); xticklabels(Lab);
ylabel('# Entries');
title('Number of CAP entries across subjects');

fig_out = fullfile('/home/nsabet/TbCAPS/SavedData', ...
    sprintf('CAP_summary_RDoC_controls_pre_%s.png', datestr(now,'yyyymmdd')));
saveas(gcf, fig_out);
fprintf('Saved QC figure to %s\n', fig_out);
