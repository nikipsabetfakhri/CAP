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

%% Beging

fprintf('\n=== CAP clustering FULL RUN (39 subjects) ===\n');
fprintf('Xon size = %d frames × %d voxels\n', size(Xon,1), size(Xon,2));
fprintf('Mask3D voxels = %d | Flattened mask = %d\n', sum(mask3D(:)), sum(mask));

% ---------------------------
% Logging
% ---------------------------
logfile = fullfile('/home/nsabet/TbCAPs/SavedData', ...
    sprintf('CAP_clustering_log_%s.txt', datestr(now,'yyyymmdd_HHMM')));
diary(logfile);
diary on

% ---------------------------
% Clustering parameters
% ---------------------------
K_opt = 4;        % or your chosen K (from consensus elbow)
Pp    = 100;      % % positive voxels retained
Pn    = 100;      % % negative voxels retained
n_rep = 20;       % number of k-means repetitions (bump to 50+ if HPC)

tic;
[CAP,~,~,idx] = Run_Clustering(Xon', K_opt, mask3D, brain_info, Pp, Pn, n_rep, [], []);
toc;
fprintf('Clustering complete: %d CAPs, %d frames clustered\n', K_opt, length(idx));

% ---------------------------
% Z-score CAP maps (optional, for visualization/interpretation)
% ---------------------------
CAP_z = CAP_Zscore(CAP);

% Save CAP maps as NIfTIs
out_dir = '/home/nsabet/TbCAPs/SavedData';
for k = 1:K_opt
    Vout = brain_info;
    V = zeros(brain_info.dim); 
    V(mask3D) = CAP(k,:);   % map cluster k back to brain volume
    Vout.fname = fullfile(out_dir, ...
        sprintf('CAP_full_controls_pre_%s_k%d.nii', datestr(now,'yyyymmdd'), k));
    spm_write_vol(Vout, V);
end

% ---------------------------
% Compute metrics
% ---------------------------
TR = 2;   % seconds (set your true TR here!)
[ExpressionMap, Counts, Entries, Avg_Duration, Duration, TransitionProbabilities, ...
    From_Baseline, To_Baseline, Baseline_resilience, Resilience, Betweenness, ...
    InDegree, OutDegree, SubjectEntries] = ...
    Compute_Metrics_simpler(idx, Indices.kept.active, Indices.scrubbed, K_opt, TR);

% ---------------------------
% QC for transitions
% ---------------------------
disp('Transition probability matrix (first CAP):');
disp(TransitionProbabilities(:,:,1));

if all(TransitionProbabilities(:) == 0)
    warning('All transition probabilities are zero! Something may be wrong.');
else
    fprintf('Transitions look valid (not all zeros).\n');
end

% ---------------------------
% Save outputs (checkpoint + results)
% ---------------------------
cap_out = fullfile(out_dir, ...
    sprintf('CAP_basis_full_controls_pre_%s.mat', datestr(now,'yyyymmdd')));
save(cap_out, 'CAP','CAP_z','idx','mask','mask3D','brain_info', ...
    'K_opt','Pp','Pn','n_rep','-v7.3');
fprintf('Saved CAP basis to %s\n', cap_out);

metrics_out = fullfile(out_dir, ...
    sprintf('CAP_metrics_full_controls_pre_%s.mat', datestr(now,'yyyymmdd')));
save(metrics_out, ...
    'ExpressionMap','Counts','Entries','Avg_Duration','Duration', ...
    'TransitionProbabilities','From_Baseline','To_Baseline', ...
    'Baseline_resilience','Resilience','Betweenness','InDegree','OutDegree', ...
    'SubjectEntries','TR','K_opt','idx','Indices','-v7.3');
fprintf('Saved metrics to %s\n', metrics_out);

fprintf('=== FULL RUN complete ===\n\n');
diary off

% Tranfer to mac via terminal
% scp nsabet@128.248.41.45:/home/nsabet/TbCAPS/SavedData/CAP_basis_RDoC_controls_pre_20250929.mat \
%     nsabet@128.248.41.45:/home/nsabet/TbCAPS/SavedData/CAP_metrics_RDoC_controls_pre_20250929.mat \
%     nsabet@128.248.41.45:/home/nsabet/TbCAPS/SavedData/CAP_z_RDoC_controls_pre_20250929.mat \
%     'nsabet@128.248.41.45:/home/nsabet/TbCAPS/SavedData/CAP_controls_pre_20250929_*.nii' \
%     nsabet@128.248.41.45:/home/nsabet/TbCAPS/SavedData/CAP_summary_RDoC_controls_pre_20250929.png \
%     /Users/npsabet/TbCAPS/SavedData/


