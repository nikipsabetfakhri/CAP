%% CAP discovery script (seed-free, TbCAPs, no GUI)

% Add TbCAPs + SPM + data to the path

clear; close all; clc; 
% Set seed
rng('default');
rng(0);

% ---------- Location toggle ----------
% Loc = 1 → Mac
% Loc = 2 → Server
Loc = 1;

if Loc == 1
    addpath('/Users/npsabet/TbCAPS');   
    addpath('/Users/npsabet/spm12');
    data_dir = '/Users/npsabet/TbCAPs/MacTest_Control_Pre_3/';
    save_dir = '/Users/npsabet/TbCAPS/SavedData';
elseif Loc == 2
    addpath('/home/nsabet/TbCAPS');   
    addpath('/home/nsabet/spm12');
    data_dir = '/home/nsabet/RDoC_CAP/RDoC_Controls_Pre/';
    save_dir = '/home/nsabet/TbCAPS/SavedData';
else
    error('Invalid Loc value. Use 1=Mac, 2=Server.');
end

% Define population colors (like the GUI) 
PopColor{1} = [255,255,180; 219,224,252; 188,252,188; 230,230,230] / 255; % light colors (background)
PopColor{2} = [130,48,48; 51,75,163; 59,113,86; 0,0,0] / 255;             % dark colors (foreground)

%% 1. Loading the data files 

% Get subject folders
subj_dirs = dir(fullfile(data_dir,'swau*'));
subj_dirs = subj_dirs([subj_dirs.isdir]);
nSubj = numel(subj_dirs);

% Track dataset count (like GUI)
n_datasets = 1;                 % this is the first dataset
ReferencePopulation = 1;        % by default, dataset 1 is the reference

% Inspect first subject's folder content
subj_dir = fullfile(data_dir, subj_dirs(1).name);
% disp(subj_dir);

% % Show what spm_select returns
% files = spm_select('FPList', subj_dir, '.*\.nii$');

% Load subject header to use as reference space (first frame)
ref_img  = fullfile(subj_dir, filesep, [subj_dirs(1).name '_00001.nii']);
Vref     = spm_vol(ref_img);

% Load TbCAPS default GM mask
mask_in  = fullfile(fileparts(which('CAP_TB')), 'DefaultData', 'Default_mask.nii');
Vmask    = spm_vol(mask_in);
mask_vol = spm_read_vols(Vmask);

% Threshold at 0.9 (GUI behavior)
mask_vol(mask_vol < 0.9) = 0;
mask_vol(mask_vol >= 0.9) = 1;

% Resample mask into subject space
fprintf('Resampling mask to match subject space...\n');
mask_resampled = CAP_V2V(mask_vol, Vmask.dim, Vmask.mat, Vref.dim, Vref.mat);

% Store mask for this dataset (logical 1D vector)
mask3D = logical(mask_resampled);
mask{n_datasets} = mask3D(:);
clear mask_vol mask_resampled mask3D Vmask

% Store brain_info for this dataset
brain_info{n_datasets} = Vref;

% Sanity check
fprintf('Mask retains %d voxels out of total %d in volume grid.\n', ...
        sum(mask{n_datasets}), prod(Vref.dim));

% Load underlay + brain 
Underlay = load_nii(fullfile(fileparts(which('CAP_TB')), 'DefaultData', 'Underlay.nii'));
Underlay_mat = [Underlay.hdr.hist.srow_x;
                Underlay.hdr.hist.srow_y;
                Underlay.hdr.hist.srow_z;
                0 0 0 1];
Underlay_dim = Underlay.hdr.dime.dim(2:4);
Underlay_info.dim = Underlay_dim;
Underlay_info.mat = Underlay_mat;
clear Underlay Underlay_mat Underlay_dim

load(fullfile(fileparts(which('CAP_TB')), 'DefaultData', 'brain.mat')); % loads 'brain'
brain_ref = brain; clear brain

fprintf('Dataset %d loaded: %d subjects, mask voxels = %d\n', ...
        n_datasets, nSubj, sum(mask{n_datasets}));

%% 2) Build TC and FD 

% Containers follow GUI style:
%   TC{ds}{subj}   -> [nTP × nMaskedVox]
%   FD{ds}         -> [nTP × nSubj]
%   n_subjects{ds} -> number of subjects in dataset ds
TC          = cell(1, n_datasets);
FD          = cell(1, n_datasets);
n_subjects  = cell(1, n_datasets);

% Use same file prefix logic as GUI (handles.prefix = 'sw')
prefix = '^sw.*\.nii$';

for ds = 1:n_datasets
    fprintf('\n=== Processing dataset %d/%d ===\n', ds, n_datasets);

    n_subjects{ds} = nSubj;               
    TC{ds} = cell(1, n_subjects{ds});     

    for s = 1:n_subjects{ds}
        subj_dir = fullfile(data_dir, subj_dirs(s).name);
        fprintf('\n[Dataset %d | Subject %d/%d] %s\n', ...
            ds, s, n_subjects{ds}, subj_dirs(s).name);

        % --- Load functional volumes (all frames) ---
        files = spm_select('FPList', subj_dir, prefix);
        if isempty(files)
            error('No functional files found for %s with prefix %s', ...
                subj_dirs(s).name, prefix);
        end
        V = spm_vol(files);

        % Read all frames, convert to single
        Y = single(spm_read_vols(V));       % [X Y Z T]

        % Reshape to 2D (frames × voxels), apply mask immediately
        Y2D = reshape(Y, [], size(Y,4))';   % [nTP × nVox]
        tmp_data = Y2D(:, mask{ds});        % [nTP × nMaskedVox]

        % Free heavy intermediates early
        clear Y Y2D

        % --- Z-score per voxel across time ---
        tmp_data = single(zscore(double(tmp_data), 0, 1));
        tmp_data(isnan(tmp_data)) = 0;

        % Store in TC (GUI expects this cell structure)
        TC{ds}{s} = tmp_data;

        nTP = size(tmp_data,1);

        % --- Motion FD ---
        motfile = dir(fullfile(subj_dir, 'rp_*.txt'));
        if isempty(motfile)
            warning('No motion file for %s. Using zeros.', subj_dirs(s).name);
            thisFD = zeros(nTP,1,'single');
        else
            thisFD = CAP_ComputeFD(fullfile(subj_dir, motfile(1).name));
            thisFD = single(thisFD);   % keep consistent
            if length(thisFD) ~= nTP
                error('FD length (%d) ≠ nTP (%d) for subject %s', ...
                    length(thisFD), nTP, subj_dirs(s).name);
            end
        end

        % Build FD matrix
        if s == 1
            FD{ds} = zeros(nTP, n_subjects{ds}, 'single');
        end
        FD{ds}(:, s) = thisFD;

        % QC print
        fprintf('   TC = %d×%d | FD = %d frames | mean FD = %.3f mm\n', ...
            size(TC{ds}{s},1), size(TC{ds}{s},2), nTP, mean(thisFD));
    end

    % --- Store dataset-wide sizes (once) ---
    if ds == 1
        SubjSize.TP  = size(TC{ds}{1}, 1);
        SubjSize.VOX = size(TC{ds}{1}, 2);
    end

    % --- Resample underlay once ---
    if ds == 1
        brain = CAP_V2V(brain_ref, ...
            Underlay_info.dim, Underlay_info.mat, ...
            brain_info{1}.dim, brain_info{1}.mat);
    end

    % --- Consistency check ---
    [IsOK, Problems] = CAP_IsDataOK(TC{ds}, FD{ds}, mask{ds}, brain_info{ds});
    if ~IsOK
        error('Data consistency check failed: %s', Problems);
    else
        fprintf('Data check passed: %s\n', Problems);
    end

    % --- Voxel-count check ---
    nMasked = sum(mask{ds});
    for s = 1:n_subjects{ds}
        if size(TC{ds}{s},2) ~= nMasked
            error('Mismatch: Subject %d has %d voxels, mask has %d', ...
                s, size(TC{ds}{s},2), nMasked);
        end
    end
    fprintf('All subjects consistent: each TC has %d masked voxels.\n', nMasked);
end

%% 3. Selecting the frames to analyse (seed-free, GUI-faithful)

Tmot = 0.5;    % motion scrubbing threshold [mm]
fprintf('\n=== Frame selection (seed-free; Tmot = %.2f mm) ===\n', Tmot);

% Containers
Xonp        = cell(1, n_datasets);   % per-dataset retained frames (cell of subj-cells)
RetainedPct = cell(1, n_datasets);   % per-dataset retained percentages
FrameIdx    = cell(1, n_datasets);   % per-dataset frame indices

for ds = 1:n_datasets
    fprintf('\n--- Dataset %d/%d ---\n', ds, n_datasets);

    % Run TbCAPs helper (same as GUI’s SeedFreeButton_Callback)
    [Xonp{ds}, p_ds, Indices_ds] = CAP_find_activity_SeedFree(TC{ds}, FD{ds}, Tmot);

    % Convert each subject’s matrix to single right away to save memory
    for s = 1:numel(Xonp{ds})
        Xonp{ds}{s} = single(Xonp{ds}{s});
    end

    % Store outputs
    RetainedPct{ds} = p_ds(3,:);   
    FrameIdx{ds}    = Indices_ds;  

    % QC: per-subject frame retention
    nTP = size(FD{ds},1);
    for s = 1:n_subjects{ds}
        nKept  = sum(Indices_ds.kept.active(:,s));
        nScrub = sum(Indices_ds.scrubbed(:,s));
        fprintf('  Subj %2d: kept %3d / %3d (%.1f%%), scrubbed %3d\n', ...
            s, nKept, nTP, 100*nKept/nTP, nScrub);
    end
    fprintf('  Dataset %d total retained frames = %d\n', ...
        ds, sum(cellfun(@(x) size(x,2), Xonp{ds})));
end

% =============================
% Pool frames across reference group ONLY
% =============================
fprintf('\nPooling retained frames from ReferencePopulation only (dataset %d)...\n', ReferencePopulation);

% Instead of giant horzcat in double precision → use single
Xon_ref = single(horzcat(Xonp{ReferencePopulation}{:})');   % [frames × voxels]

fprintf('Pooled Xon_ref size = %d frames × %d voxels\n', size(Xon_ref,1), size(Xon_ref,2));

% QC: voxel count vs mask
assert(size(Xon_ref,2) == sum(mask{ReferencePopulation}), ...
    'Mismatch: Xon_ref has %d voxels, mask has %d', ...
    size(Xon_ref,2), sum(mask{ReferencePopulation}));

% Handle NaNs / Infs efficiently
badCols = find(any(~isfinite(Xon_ref),1));
if ~isempty(badCols)
    fprintf('Warning: %d voxel columns had NaN/Inf; setting to 0.\n', numel(badCols));
    Xon_ref(:,badCols) = 0;
end

fprintf('Final Xon_ref (for clustering) = %d frames × %d voxels\n', ...
    size(Xon_ref,1), size(Xon_ref,2));

% =============================
% Frame retention violin plots
% =============================
% tmp_toplot = ConcatMat(RetainedPct, n_datasets, 1, n_subjects, 'FD');
% group_labels = arrayfun(@(ds) sprintf('Group %d', ds), 1:n_datasets, 'UniformOutput', false);
% 
% figure(1);
% ax = gca;  
% [~,~,TPViolin] = MakeViolin(tmp_toplot, ax, group_labels, ...
%     'Frames ret. [%]', PopColor, n_datasets, 1);
% 
% set(TPViolin,'Visible','on');
% title(ax, 'Frame retention across subject groups');


 %% QC BEFORE CONSENSUS
% fprintf('\n================= QC CHECKS (pre-consensus) =================\n');
% 
% for ds = 1:n_datasets
%     fprintf('\n-- Dataset %d --\n', ds);
% 
%     % 1) Dimensions of TC{ds}{s} and FD{ds}
%     for s = 1:n_subjects{ds}
%         fprintf('  Subj %2d: TC = %3d frames × %6d voxels | mean FD = %.3f mm\n', ...
%             s, size(TC{ds}{s},1), size(TC{ds}{s},2), mean(FD{ds}(:,s)));
%     end
%     fprintf('  FD{ds} size = %d frames × %d subjects\n', size(FD{ds},1), size(FD{ds},2));
% 
%     % 2) Xonp{ds} structure (per subject, vox × keptFrames)
%     fprintf('  Xonp{%d} is a %s with %d subjects\n', ds, class(Xonp{ds}), numel(Xonp{ds}));
%     for s = 1:min(5,numel(Xonp{ds}))
%         fprintf('    Subj %2d retained: %6d vox × %4d frames\n', ...
%             s, size(Xonp{ds}{s},1), size(Xonp{ds}{s},2));
%     end
% 
%     % 3) Mask alignment
%     vox_per_subj = unique(cellfun(@(x) size(x,1), Xonp{ds}));
%     fprintf('  Unique voxel counts in Xonp{%d} = %s | Mask voxels = %d\n', ...
%         ds, mat2str(vox_per_subj), sum(mask{ds}));
%     assert(all(vox_per_subj == sum(mask{ds})), ...
%         'Voxel mismatch: Xonp{%d} voxels differ from mask.', ds);
% 
%     % 4) Retention summary
%     nTP = size(FD{ds},1);
%     total_kept = sum(cellfun(@(x) size(x,2), Xonp{ds}));
%     fprintf('  Total retained frames = %d (%.1f%% of all frames)\n', ...
%         total_kept, 100*total_kept/(nTP*n_subjects{ds}));
% 
%     % 5) Sanity: recompute pooled rows = sum kept
%     Xon_ds = horzcat(Xonp{ds}{:})';  % [frames × voxels]
%     assert(size(Xon_ds,1) == total_kept, 'Pooled rows ≠ sum of kept frames (ds=%d).', ds);
% 
%     % 6) NaN/Inf check on pooled (should be already clean)
%     badCols = find(any(~isfinite(Xon_ds),1));
%     if ~isempty(badCols)
%         fprintf('  Warning: %d NaN/Inf voxel columns in pooled ds=%d → set to 0.\n', numel(badCols), ds);
%         Xon_ds(:,badCols) = 0; 
%     end
% end
% 
% % 7) Reference group quick echo
% fprintf('\nReferencePopulation = dataset %d\n', ReferencePopulation);
% Xon_ref = horzcat(Xonp{ReferencePopulation}{:})';
% fprintf('Xon_ref = %d frames × %d voxels | Mask voxels = %d\n', ...
%     size(Xon_ref,1), size(Xon_ref,2), sum(mask{ReferencePopulation}));
% assert(size(Xon_ref,2) == sum(mask{ReferencePopulation}), 'Xon_ref cols ≠ mask voxels');
% 
% fprintf('================= END QC =================\n\n');

%% Memory check before consensus 

S = whos;
fprintf('Workspace memory use: %.2f GB\n', sum([S.bytes]) / 1e9);


%% 4. Consensus clustering (determine optimum K)

% Log console output to file
logfile = fullfile(save_dir, sprintf('CAP_log_%s.txt', datestr(now,'yyyymmdd_HHMM')));
diary(logfile);
diary on

% Parameters (TEST run)
Pcc     = 100;      % % of frames per fold
N       = 2;        % number of folds
K_range = 2:3;      % test range (expand later)

fprintf('Subjects in dataset = %d\n', n_subjects{ReferencePopulation});
fprintf('K_range = %s | Folds = %d | Subsample fraction = %.2f\n', ...
    mat2str(K_range), N, Pcc/100);

tStart = tic;
Consensus = [];

for k = K_range
    fprintf('\n>>> Running consensus clustering for K = %d ...\n', k);
    tK = tic;

    % Like GUI
    [Ck] = CAP_ConsensusClustering( ...
        Xonp{ReferencePopulation}, k, 'items', Pcc/100, N, 'correlation');

    Ck = single(Ck);   % <<< force single precision to halve memory

    % Allocate and store
    if isempty(Consensus)
        Consensus = zeros(size(Ck,1), size(Ck,2), numel(K_range), 'single');
    end
    Consensus(:,:,K_range==k) = Ck;

    fprintf('...done (%.2f sec)\n', toc(tK));
    clear Ck
end
fprintf('\nAll consensus clustering finished in %.2f seconds.\n', toc(tStart));

% Compute quality metrics
[~, Qual] = ComputeClusteringQuality(Consensus, K_range);

% % Plot quality curves
figure(2);
plot(K_range, Qual(2,1:numel(K_range)), '-o', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('K (# clusters)'); ylabel('Stability (1 - Lorena metric)');
title('Consensus clustering quality'); grid on;
 
% figure(3);
% bar(K_range, 1 - Qual(2,1:numel(K_range)), 'FaceColor',[0.3 0.6 0.9]);
% xlabel('K (# clusters)'); ylabel('Stability (1 - Lorena metric)');
% title('Consensus clustering quality (GUI style)');
% ylim([0 1]); grid on; box off;

% Save consensus results only (lightweight)
% outfile = fullfile(save_dir, sprintf('Consensus_%s.mat', datestr(now,'yyyymmdd_HHMM')));
% save(outfile, 'Consensus','K_range','N','Pcc','-v7.3');
% fprintf('\nConsensus results saved to: %s\n', outfile);
diary off;


%% 5) CAP clustering 

% Choose K from your consensus curve (using a small test value here)
K_opt = 2;          % <- later: pick from Qual
Pp    = 100;
Pn    = 100;
n_rep = 5;          % small for test; increase later

% IMPORTANT: use the same orientation the GUI uses: [n_vox × n_frames]
Xon_forClust = cell2mat(Xonp{ReferencePopulation});     % [vox × frames]
fprintf('Clustering on %d vox × %d frames\n', size(Xon_forClust,1), size(Xon_forClust,2));

% Run clustering exactly like the GUI ClusterButton does for seed-free:
[CAP, Disp, STDCAP, idx, CorrDist, sfrac] = Run_Clustering( ...
    Xon_forClust, ...                                      % [vox × frames]
    K_opt, ...
    mask{ReferencePopulation}, ...
    brain_info{ReferencePopulation}, ...
    Pp, Pn, n_rep, ...
    1, ...                                                 % idx_sep_seeds = 1 for seed-free
    'SeedFree');

fprintf('Clustering complete: %d CAPs, %d frames clustered\n', K_opt, numel(idx));

% QC after clustering (mirrors what the GUI enables you to inspect)
fprintf('\n=== QC after Run_Clustering ===\n');

% Frames per CAP
nFrames_pooled = size(Xon_forClust,2);
fprintf('idx length = %d | pooled frames = %d\n', numel(idx), nFrames_pooled);
assert(numel(idx) == nFrames_pooled, 'Mismatch: idx length ≠ pooled frames');

for k = 1:K_opt
    fprintf('CAP %d: %d frames (%.1f%%)\n', ...
        k, sum(idx==k), 100*sum(idx==k)/nFrames_pooled);
end

% Subject-level contributions (same logic as GUI uses with FrameIndices)
frame_counts = zeros(n_subjects{ReferencePopulation}, K_opt);
offset = 0;
for s = 1:n_subjects{ReferencePopulation}
    nKept = sum(FrameIdx{ReferencePopulation}.kept.active(:,s));
    subj_idx = idx(offset+1 : offset+nKept);
    offset = offset + nKept;
    for k = 1:K_opt
        frame_counts(s,k) = sum(subj_idx==k);
    end
end
disp('Frames per subject × CAP:'); disp(frame_counts);

disp(brain_info{ReferencePopulation}.dim);
disp(size(brain)); % should match brain_info dims
disp(sum(mask{ReferencePopulation}));

% Z-score CAPs and save to NIfTI 
fprintf('\nZ-scoring CAPs and saving to NIfTI...\n');
CAP_z = CAP_Zscore(CAP);

CAPToNIFTI(      CAP,   mask{ReferencePopulation}, brain_info{ReferencePopulation}, save_dir, ...
    sprintf('CAP_%s', datestr(now,'yyyymmdd_HHMM')));

CAPToNIFTI(CAP_z,       mask{ReferencePopulation}, brain_info{ReferencePopulation}, save_dir, ...
    sprintf('CAP_Zscored_%s', datestr(now,'yyyymmdd_HHMM')));

fprintf('Saved CAP and CAP_Zscored NIfTIs to %s\n', out_dir);

% Plot GUI-style slices for each CAP

fprintf('\nPlotting CAP slices (GUI)...\n');

% GUI defaults
Tv   = 0.5;   % threshold for transparency
maxC = 1.5;   % color saturation

% Slice coordinates like manual screenshot (X=0, Y=0, Z=30 mm in MNI)
slice_coords = [0 0 30];  

figure('Name','CAP slices (GUI)','Color','w');
for k = 1:K_opt
    for dim = 1:3
        subplot(K_opt,3,(k-1)*3+dim);
        switch dim
            case 1
                plot_slice(CAP_z(k,:), Tv, maxC, ...
                           mask{ReferencePopulation}, brain, ...
                           brain_info{ReferencePopulation}, ...
                           'X', slice_coords(1), gca);
                title(sprintf('CAP %d (X=%d)',k,slice_coords(1)));
            case 2
                plot_slice(CAP_z(k,:), Tv, maxC, ...
                           mask{ReferencePopulation}, brain, ...
                           brain_info{ReferencePopulation}, ...
                           'Y', slice_coords(2), gca);
                title(sprintf('CAP %d (Y=%d)',k,slice_coords(2)));
            case 3
                plot_slice(CAP_z(k,:), Tv, maxC, ...
                           mask{ReferencePopulation}, brain, ...
                           brain_info{ReferencePopulation}, ...
                           'Z', slice_coords(3), gca);
                title(sprintf('CAP %d (Z=%d)',k,slice_coords(3)));
        end
    end
end

% Add GUI-style colorbar
figure('Name','CAP colormap','Color','w');
Create_CAP_colorbar(-maxC,maxC,0.5,Tv,'Activation (Z)',gca,...
                    'Horizontal','div','RdBu',1000);


%% 6) CAP Metrics 
fprintf('\n=== CAP metrics ===\n');

TR = 2;  % <-- set to actual TR in seconds

[TPM, Counts, Number, Avg_Duration, Duration, TM, ...
    From_Baseline, To_Baseline, Baseline_resilience, Resilience, ...
    Betweenness, kin, kout, SubjectEntries] = ...
    Compute_Metrics_simpler(idx, ...
        FrameIdx{ReferencePopulation}.kept.active, ...
        FrameIdx{ReferencePopulation}.scrubbed, ...
        K_opt, TR);

% Quick inspection of TPM
disp('Unique TPM values (should be -1, 0, or 1..K):');
disp(unique(TPM(:)));

% Count per category
fprintf('\nCounts across all subjects:\n');
fprintf('  Scrubbed   = %d\n', sum(TPM(:) == -1));
fprintf('  Baseline   = %d\n', sum(TPM(:) == 0));
for k = 1:K_opt
    fprintf('  CAP %d      = %d\n', k, sum(TPM(:) == k));
end

% Check row by row
for s = 1:size(TPM,1)
    fprintf('\nSubject %d:\n', s);
    fprintf('  Scrubbed   = %d\n', sum(TPM(s,:) == -1));
    fprintf('  Baseline   = %d\n', sum(TPM(s,:) == 0));
    for k = 1:K_opt
        fprintf('  CAP %d      = %d\n', k, sum(TPM(s,:) == k));
    end
end

% Cross-check against SubjectEntries
disp('SubjectEntries (first few rows):');
disp(SubjectEntries(1:min(5,end),:));


%% Plot CAP metrics 

% Plot raster plot 
fprintf('\n=== Plotting CAP raster ===\n');

[nSubj, nFrames] = size(TPM);

% Define colors exactly as GUI: scrubbed=gray, baseline=white, CAPs=distinct
colors = [
    0.5 0.5 0.5;   % -1 scrubbed (gray)
    1   1   1;     % 0 baseline (white) -- won't appear in seed-free
    1   0   0;     % CAP1 red
    0   0.5 1;     % CAP2 blue
    0   0.8 0;     % CAP3 green
    0.6 0   0.8;   % CAP4 purple
    1   0.6 0;     % CAP5 orange
    0   0   0;     % CAP6 black
];

% Map TPM values (-1,0,1..K) into indices for colors
vals = unique(TPM(:));
mapVals = containers.Map([-1,0,1:K_opt], 1:(K_opt+2));

% Build RGB image for plotting
RasterImg = ones(nSubj, nFrames, 3); % init white
for s = 1:nSubj
    for f = 1:nFrames
        v = TPM(s,f);
        if isKey(mapVals,v)
            RasterImg(s,f,:) = colors(mapVals(v),:);
        end
    end
end

% Flip vertically to match GUI (subject 1 at top)
RasterImg = flipud(RasterImg);

% Plot
figure('Name','CAP dynamics (raster)','Color','w');
image(RasterImg);
title('State time courses');
xlabel('Time [frames]');
ylabel('Subjects');
set(gca,'YTick',1:nSubj,'YTickLabel',nSubj:-1:1);


% Plot transition matrices (GUI-style)

% Inputs you should already have from Compute_Metrics_simpler:
%   TM = [nStates x nStates x nSubjects] transition probabilities
%   K_opt = number of CAPs
%   ReferencePopulation = dataset index
%   subjID = subject index to show (e.g., 1, 2, ...)

subjID = 1;  % <-- change this to whichever subject you want to display

% 1) Group-averaged transition matrix
TM_group = squeeze(mean(TM,3));   % average across subjects

% 2) Single-subject transition matrix
TM_subj  = squeeze(TM(:,:,subjID));

% --- GUI crops away scrubbed/unassigned rows/cols:
TM_group_crop = TM_group(3:end-1, 3:end-1);
TM_subj_crop  = TM_subj(3:end-1, 3:end-1);

% --- Plot
figure('Name','Transition matrices (GUI-style)','Color','w');

subplot(2,1,1);
imagesc(TM_group_crop);
title('Group-averaged transition matrix');
axis square off;
colormap(flipud(gray));
caxis([0 0.03]);   % same as GUI
colorbar;

subplot(2,1,2);
imagesc(TM_subj_crop);
title(sprintf('Subject %d transition matrix', subjID));
axis square off;
colormap(flipud(gray));
caxis([0 0.03]);   % keep same scale
colorbar;


% Plot cumulative state distributions (GUI-faithful)

state_to_plot = 1;   % CAP number (1..K_opt)
TR = 2;              % replace with your TR (seconds)
[nSubj, nFrames] = size(TPM);

time_axis = (0:nFrames-1) * TR;   % time in seconds

figure('Name','Cumulative state distribution','Color','w'); hold on;

all_cum = zeros(nSubj, nFrames);

for s = 1:nSubj
    % 1 where subject is in this CAP, else 0
    binvec = (TPM(s,:) == state_to_plot);
    
    % cumulative sum over time
    all_cum(s,:) = cumsum(binvec);
    
    % thin subject line (light yellow, like GUI)
    plot(time_axis, all_cum(s,:), 'Color', [1 0.9 0.4]); 
end

% group average curve (thicker, darker brownish)
mean_cum = mean(all_cum,1);
plot(time_axis, mean_cum, 'Color', [0.5 0.2 0.2], 'LineWidth', 2);

xlabel('Time [s]');
ylabel('Cumul. sum [-]');
title(sprintf('Cumulative distribution for CAP %d', state_to_plot));
box off;


% CAP metrics

metrics = {'Counts','Number','Resilience','Betweenness','kin','kout'};
labels  = {'Raw counts','Number','Resilience','Betweenness','kin','kout'};
titles  = {'Raw counts [-]','Entries [-]','Resilience [-]', ...
           'Betweenness [-]','kin [-]','kout [-]'};

figure('Name','CAP metrics (GUI-style)','Color','w');

for m = 1:numel(metrics)
    subplot(2,3,m);

    % Collect data
    tmpM = cell(1,n_datasets);
    for ds = 1:n_datasets
        try
            tmpM{ds} = eval(metrics{m}); % e.g. Counts, Number, kin...
        catch
            tmpM{ds} = [];
        end
    end

    try
        tmp_toplot = ConcatMat(tmpM, n_datasets, K_opt, n_subjects, labels{m});
        MakeViolin(tmp_toplot, gca, ...
            arrayfun(@num2str,1:K_opt,'UniformOutput',false), ...
            titles{m}, PopColor, n_datasets, K_opt);
    catch
        title([labels{m} ' (not available)']);
    end
end