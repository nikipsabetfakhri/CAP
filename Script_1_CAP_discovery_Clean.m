%% CAP discovery script (seed-free, TbCAPs, no GUI)

% Add TbCAPs + SPM to the path
% Mac
% addpath('/Users/npsabet/TbCAPS');   
% addpath('/Users/npsabet/spm12');
% addpath('/Users/npsabet/TbCAPs/control_pre_test_subj/');

% Dell 
addpath('/home/nsabet/TbCAPs');   
addpath('/home/nsabet/spm12');
addpath('/home/nsabet/RDoC_CAP/RDoC_Controls_Pre/');

%% 1. Loading the data files

% Set input folder
data_dir = '/home/nsabet/RDoC_CAP/RDoC_Controls_Pre/';

% Get subject folders
subj_dirs = dir(fullfile(data_dir,'swau*'));
subj_dirs = subj_dirs([subj_dirs.isdir]);
nSubj = numel(subj_dirs);

% Inspect first subject's folder content
subj_dir = fullfile(data_dir, subj_dirs(1).name);
disp(subj_dir); 
ls(subj_dir);  % should list swauRDHC001_pre_00001.nii etc.

% Show what spm_select returns
files = spm_select('FPList', subj_dir, '.*\.nii$');
disp(files(1:min(size(files,1),5), :));  % first few paths

% Load subject header to use as reference space (first frame)
ref_img  = '/home/nsabet/RDoC_CAP/RDoC_Controls_Pre/swauRDHC001_pre/swauRDHC001_pre_00001.nii';
Vref     = spm_vol(ref_img);

% Load TbCAPS default GM mask and resample into subject space
mask_in  = fullfile(fileparts(which('CAP_TB')), 'DefaultData', 'Default_mask.nii');
Vmask    = spm_vol(mask_in);
mask_vol = spm_read_vols(Vmask);

fprintf('Resampling mask to match subject space...\n');
mask_resampled = CAP_V2V(mask_vol, Vmask.dim, Vmask.mat, Vref.dim, Vref.mat);

% Binarize
mask3D = mask_resampled > 0.5;
mask   = mask3D(:);   % flatten for indexing

% Sanity check
fprintf('Mask retains %d voxels out of total %d in volume grid.\n', ...
        sum(mask), prod(Vref.dim));

% Header info for clustering
brain_info = Vref;


%% 2. Build TC and FD (timecourses and framewise displacement)

TC = cell(1,nSubj);   % each cell: [nTP × nMaskedVox]
FD = cell(1,nSubj);   % each cell: [nTP × 1]

for s = 1:nSubj
    subj_dir = fullfile(data_dir, subj_dirs(s).name);
    fprintf('\n[Subject %d/%d] %s\n', s, nSubj, subj_dirs(s).name);

    % --- Load subject NIfTIs ---
    files = spm_select('FPList', subj_dir, '^swau.*\.nii$');
    V = spm_vol(files);
    Y = spm_read_vols(V);                % dims: [X Y Z T]
    Y2D = reshape(Y, [], size(Y,4))';    % [nTP × nVox]
    TC{s} = Y2D(:, mask);                % apply mask
    nTP = size(TC{s},1);

    % --- Load motion file → FD{s} ---
    motfile = dir(fullfile(subj_dir, 'rp_*.txt'));
    if isempty(motfile)
        warning('No motion file for %s. Using zeros.', subj_dirs(s).name);
        FD{s} = zeros(nTP,1);
    else
        motfile_name = fullfile(subj_dir, motfile(1).name);
        FD{s} = CAP_ComputeFD(motfile_name);   % [nTP × 1]
        if length(FD{s}) ~= nTP
            error('FD length (%d) ≠ nTP (%d) for subject %s', ...
                   length(FD{s}), nTP, subj_dirs(s).name);
        end
    end

    % QC print per subject
    fprintf('Subject %d: TC = %d×%d, FD = %d frames\n', ...
        s, size(TC{s},1), size(TC{s},2), length(FD{s}));
end

% Collapse FD correctly into numeric matrix [nTP × nSubj]
FDmat = cat(2, FD{:});
fprintf('\nFDmat size = %d × %d\n', size(FDmat,1), size(FDmat,2));

% Sanity checks
for s = 1:nSubj
    fprintf('Subject %d: mean FD = %.3f mm\n', s, mean(FDmat(:,s)));
end

[IsOK, Problems] = CAP_IsDataOK(TC, FDmat, mask, brain_info);
if ~IsOK
    error('Data consistency check failed: %s', Problems);
else
    fprintf('Data check passed: %s\n', Problems);
end

nMasked = sum(mask);
for s = 1:nSubj
    if size(TC{s},2) ~= nMasked
        error('Mismatch: Subject %d has %d voxels, mask has %d voxels', ...
              s, size(TC{s},2), nMasked);
    end
end
fprintf('All subjects consistent: each TC has %d masked voxels.\n', nMasked);

clear Y Y2D


%% 3. Selecting the frames to analyse (pool + normalize + QC)

Tmot = 0.5;    % scrubbing threshold (mm)

% Your CAP_find_activity_SeedFree returns per-subject cells (Xonp)
[Xonp,~,Indices,~] = CAP_find_activity_SeedFree(TC, FDmat, Tmot);

% Pool retained frames across subjects → [totalKeptFrames × nVox]
Xon = horzcat(Xonp{:})';   % concat columns, then transpose
fprintf('\nPooled Xon size = %d frames × %d voxels\n', size(Xon,1), size(Xon,2));

% QC: voxel count should match mask
assert(size(Xon,2) == nMasked, ...
    'Mismatch: Xon has %d voxels, mask has %d.', size(Xon,2), nMasked);

% Per-subject kept/scrubbed frames
nTP = size(FDmat,1);
fprintf('\n=== Frame Retention Summary (Tmot = %.2f mm) ===\n', Tmot);
for s = 1:nSubj
    nKept  = sum(Indices.kept.active(:,s));
    nScrub = sum(Indices.scrubbed(:,s));
    fprintf('Subj %2d: kept %3d / %3d (%.1f%%), scrubbed %3d\n', ...
        s, nKept, nTP, 100*nKept/nTP, nScrub);
end
fprintf('Total retained frames = %d (%.1f%% of all frames)\n', ...
    size(Xon,1), 100*size(Xon,1)/(nTP*nSubj));

% Normalize voxel columns across frames for clustering
Xon = zscore(Xon, 0, 1);

% QC: handle NaN/Inf columns (zero-variance voxels)
badCols = find(any(~isfinite(Xon),1));
if ~isempty(badCols)
    fprintf('Warning: %d voxel columns had NaN/Inf after z-score; setting to 0.\n', numel(badCols));
    Xon(:,badCols) = 0;
end

fprintf('Final Xon size (after z-score) = %d frames × %d voxels\n', ...
    size(Xon,1), size(Xon,2));

% (Optional) checkpoint save before consensus clustering
% save(fullfile('/home/nsabet/TbCAPS/SavedData', ...
%     sprintf('PreConsensus_%s.mat', datestr(now,'yyyymmdd_HHMM'))), ...
%     'Xon','Xonp','Indices','FDmat','mask','mask3D','brain_info','Tmot','-v7.3');

%% QC block for CAP preprocessing (after Section 3)

fprintf('\n================= QC CHECKS =================\n');

% 1. Dimensions of TC per subject
for s = 1:nSubj
    fprintf('Subj %2d: TC = %3d frames × %6d voxels\n', ...
        s, size(TC{s},1), size(TC{s},2));
end

% 2. FD matrix size
fprintf('\nFDmat size = %d frames × %d subjects\n', size(FDmat,1), size(FDmat,2));

% 3. Xon pooled size
fprintf('Pooled Xon size = %d frames × %d voxels\n', size(Xon,1), size(Xon,2));

% 4. Mask voxel count check
fprintf('Mask voxels = %d | Xon columns = %d\n', sum(mask), size(Xon,2));

% 5. Per-subject frame retention check
nTP = size(FDmat,1);
for s = 1:nSubj
    nKept  = sum(Indices.kept.active(:,s));
    nScrub = sum(Indices.scrubbed(:,s));
    fprintf('Subj %2d: kept %3d / %3d (%.1f%%), scrubbed %3d\n', ...
        s, nKept, nTP, 100*nKept/nTP, nScrub);
end
fprintf('Total retained frames in Xon = %d\n', size(Xon,1));

% 6. Check alignment with raw TC (before normalization)
[Xonp_raw,~,Indices_raw,~] = CAP_find_activity_SeedFree(TC, FDmat, Tmot);
Xon_raw = horzcat(Xonp_raw{:})';   % pooled [frames × voxels]

frame1 = TC{1}(Indices_raw.kept.active(:,1), :);  % subj 1 kept frames
if isequal(frame1(1,:), Xon_raw(1,:))
    fprintf('Alignment check: PASS (first retained frame matches pooled Xon_raw)\n');
else
    fprintf('Alignment check: FAIL (frame mismatch before normalization)\n');
end

% 7. Check normalization statistics
mean_val = mean(Xon(:));
std_val  = std(Xon(:));
fprintf('Xon global mean ~ %.3f, std ~ %.3f (expected ~0, ~1)\n', mean_val, std_val);

fprintf('================= END QC =================\n\n');

%% 4. Consensus clustering (determine optimum K)

% --- Log console output to file ---
dairyfile = fullfile('/home/nsabet/TbCAPS/SavedData', ...
    sprintf('CAP_log_%s.txt', datestr(now,'yyyymmdd_HHMM')));
diary(dairyfile);
diary on

% Parameters
Pcc = 80;      % % of frames per fold (use 100 for small test runs)
N   = 20;      % number of folds
K_range = 2:8; % cluster range to test

% --- Run consensus clustering (overnight job) ---
% NOTE: Comment out while testing. Run only when ready.
[Consensus] = CAP_ConsensusClustering(Xon, K_range, 'items', Pcc/100, N, 'correlation');
[~, Qual]   = ComputeClusteringQuality(Consensus,[]);

% --- Inspect Qual after run ---
% figure;
% plot(K_range, Qual,'-o','LineWidth',2);
% xlabel('K'); ylabel('Quality metric');
% title('Consensus clustering quality');

% --- CAP_Zscore note ---
% You only need CAP_Zscore *after clustering*, when you have the CAP maps.
% Right now you don’t have CAP maps (that comes from k-means inside consensus).
% So keep this line commented out until then:
CAPs_z = CAP_Zscore(CAPs);

% --- Save checkpoint (metadata only, lighter than Xon) ---
save(fullfile('/home/nsabet/TbCAPs/SavedData', ...
    sprintf('PreConsensus_%s.mat', datestr(now,'yyyymmdd_HHMM'))), ...
    'Indices','FDmat','mask','brain_info','Tmot','-v7.3');

% When you run consensus overnight, you can save Consensus as well:
% save(fullfile('/home/nsabet/TbCAPs/SavedData', ...
%     sprintf('Consensus_%s.mat', datestr(now,'yyyymmdd_HHMM'))), ...
%     'Consensus','K_range','N','Pcc','-v7.3');

