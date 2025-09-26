%% CAP analysis script (seed-free, TbCAPs, no GUI)

% Add TbCAPs + SPM to the path
addpath('/home/nsabet/TbCAPS');   
addpath('/home/nsabet/spm12');
addpath('/home/nsabet/RDoC_CAP');

%% -1. Resample default mask 

% Load subject header to use as reference space
ref_img  = '/home/nsabet/RDoC_CAP/RDoC_Controls_Pre/swauRDHC001_pre/swauRDHC001_pre_00001.nii';
Vref     = spm_vol(ref_img);

% Load TbCAPS default GM mask
mask_in  = fullfile(fileparts(which('CAP_TB')), 'DefaultData', 'Default_mask.nii');
Vmask    = spm_vol(mask_in);

% Check if dimensions match
if ~isequal(Vref.dim, Vmask.dim)
    fprintf('Mask dimensions (%s) do not match subject (%s). Resampling...\n', ...
        mat2str(Vmask.dim), mat2str(Vref.dim));
    
    % Reslice mask to subject resolution
    spm_reslice({ref_img, mask_in}, struct('mean',false,'which',1,'interp',0));
    
    % spm_reslice creates "rDefault_mask.nii"
    [mask_path,~,~] = fileparts(mask_in);
    mask_out = fullfile(mask_path, 'Default_mask_resliced.nii');
    movefile(fullfile(mask_path, ['r' spm_file(mask_in,'filename')]), mask_out);
    
    Vmask = spm_vol(mask_out); % reload the resampled mask
else
    mask_out = mask_in; % no resampling needed
end

% Finalize mask
mask_vol  = spm_read_vols(Vmask);
mask      = mask_vol > 0.5;  % binarize
brain_info = Vref;           % use subject header for clustering

%% 1. Loading the data files

% -------------------------------------------------------------------------
% TbCAPs expects:
%   TC   = cell array of [timepoints x voxels] per subject
%   mask = logical vector of voxels to keep (e.g., GM mask)
%   brain_info = spm_vol struct from a representative NIfTI
%   FD   = [timepoints x subjects] framewise displacement matrix
% -------------------------------------------------------------------------

% Set input folder
data_dir = '/home/nsabet/RDoC_CAP/RDoC_Controls_Pre';

% Get subject folders
subj_dirs = dir(fullfile(data_dir,'swau*'));
subj_dirs = subj_dirs([subj_dirs.isdir]);

% Load default GM mask shipped with TbCAPS
mask_file = fullfile(fileparts(which('CAP_TB')), 'DefaultData', 'Default_mask.nii');
Vmask     = spm_vol(mask_file);
mask_vol  = spm_read_vols(Vmask);
mask      = mask_vol > 0.5;  % binarize
brain_info = Vmask;          % header info for clustering

% Inspect first subject's folder content
subj_dir = fullfile(data_dir, subj_dirs(1).name);
disp(subj_dir); 
ls(subj_dir);  % should list swauRDHC001_pre_00001.nii etc.

% Show what spm_select returns
files = spm_select('FPList', subj_dir, '.*\.nii$');
disp(files(1:min(size(files,1),5), :));  % first few paths

% Load subject header to use as reference space
ref_img  = '/home/nsabet/RDoC_CAP/RDoC_Controls_Pre/swauRDHC001_pre/swauRDHC001_pre_00001.nii';
Vref     = spm_vol(ref_img);

% Load TbCAPS default GM mask
mask_in  = fullfile(fileparts(which('CAP_TB')), 'DefaultData', 'Default_mask.nii');
Vmask    = spm_vol(mask_in);

% Check if dimensions match
if ~isequal(Vref.dim, Vmask.dim)
    fprintf('Mask dimensions (%s) do not match subject (%s). Resampling...\n', ...
        mat2str(Vmask.dim), mat2str(Vref.dim));
    
    % Reslice mask to subject resolution
    spm_reslice({ref_img, mask_in}, struct('mean',false,'which',1,'interp',0));
    
    % spm_reslice creates "rDefault_mask.nii"
    [mask_path,~,~] = fileparts(mask_in);
    mask_out = fullfile(mask_path, 'Default_mask_resliced.nii');
    movefile(fullfile(mask_path, ['r' spm_file(mask_in,'filename')]), mask_out);
    
    Vmask = spm_vol(mask_out); % reload the resampled mask
else
    mask_out = mask_in; % no resampling needed
end

mask_vol  = spm_read_vols(Vmask);
mask      = mask_vol > 0.5;  % binarize
brain_info = Vref;           % use subject header for clustering

%% 2. Build TC and FD with dimension checks
%% 2. Build TC and FD correctly (cell arrays, one per subject)
TC = cell(1, numel(subj_dirs));
FD = cell(1, numel(subj_dirs));

for s = 1:numel(subj_dirs)
    subj_dir = fullfile(data_dir, subj_dirs(s).name);
    fprintf('\n[Subject %d/%d] %s\n', s, numel(subj_dirs), subj_dirs(s).name);

    % ---------------------------
    % Load subject NIfTIs
    % ---------------------------
    nii_list = dir(fullfile(subj_dir, '*.nii'));
    if isempty(nii_list)
        warning('No NIfTI files found for %s. Skipping subject.', subj_dirs(s).name);
        TC{s} = [];
        FD{s} = [];
        continue;
    end
    
    V = spm_vol(char(fullfile(subj_dir, {nii_list.name})));
    Y = spm_read_vols(V);
    if ndims(Y) ~= 4
        warning('Subject %s NIfTI not 4D. Skipping.', subj_dirs(s).name);
        TC{s} = [];
        FD{s} = [];
        continue;
    end
    Y2D = reshape(Y, [], size(Y,4))';   % [timepoints × voxels]
    TC{s} = Y2D(:, mask(:));

    % ---------------------------
    % Load motion file → FD{s}
    % ---------------------------
    motfile = dir(fullfile(subj_dir, 'rp_*.txt'));
    if isempty(motfile)
        warning('No motion file for %s. Skipping FD!', subj_dirs(s).name);
        FD{s} = [];
        continue;
    end
    FDsubj = load(fullfile(subj_dir, motfile(1).name));
    FDsubj = [zeros(1,6); diff(FDsubj)];   % frame-to-frame diffs
    FD{s}  = sqrt(sum(FDsubj(:,1:3).^2,2));  % Euclidean displacement per frame
end

% Drop empties
valid_idx = ~cellfun(@isempty, TC);
TC = TC(valid_idx);
FD = FD(valid_idx);

fprintf('\nFinished loading. Retained %d subjects.\n', sum(valid_idx));



%% 2. Specifying the main parameters

Tmot = 0.5;    % scrubbing threshold
Pp   = 100;    % % positive voxels to retain
Pn   = 100;    % % negative voxels to retain
n_rep = 10;    % K-means replicates
Pcc   = 80;    % % items per consensus fold
N     = 10;    % number of consensus folds

save_dir = '/home/nsabet/TbCAPS/SavedData/Control_test';
if ~exist(save_dir,'dir'), mkdir(save_dir); end

%% 3. Selecting the frames to analyse    

[Xon,~,Indices] = CAP_find_activity_SeedFree(TC,FD,Tmot);


%% 4. Consensus clustering (parallel)

% Choose cluster range
K_range = 2:6;

% Start parallel pool if not already running
numWorkers = 20;  % adjust based on system
if isempty(gcp('nocreate'))
    parpool('local', numWorkers);
end

Consensus = cell(1, numel(K_range));
Qual      = cell(1, numel(K_range));

parfor i = 1:numel(K_range)
    K = K_range(i);
    fprintf('Running consensus clustering for K=%d\n', K);

    % Run clustering
    C = CAP_ConsensusClustering(Xon, K, 'items', Pcc/100, N, 'correlation');
    [~,Q] = ComputeClusteringQuality(C, []);

    % Store results
    Consensus{i} = C;
    Qual{i}      = Q;

    % Save checkpoint per K
    save(sprintf('Consensus_K%d.mat', K), 'C', 'Q', 'K', '-v7.3');
end

% Convert to arrays
Consensus = cat(3, Consensus{:});
Qual      = cat(2, Qual{:});

% Save all results together
save('Consensus_All.mat', 'Consensus', 'Qual', 'K_range', '-v7.3');
fprintf('All consensus clustering finished and saved.\n');


% Inspect Qual to pick best K

K_opt = 3;   % <-- manually set once you decide

%% 5. Clustering into CAPs

[CAP,~,~,idx] = Run_Clustering(cell2mat(Xon),...
    K_opt,mask,brain_info,Pp,Pn,n_rep,[],[]);

%% 6. Computing metrics

TR = 2;   % seconds

[ExpressionMap,Counts,Entries,Avg_Duration,Duration,TransitionProbabilities,...
 From_Baseline,To_Baseline,Baseline_resilience,Resilience,Betweenness,...
 InDegree,OutDegree,SubjectEntries] = Compute_Metrics_simpler(idx,...
    Indices.kept.active,Indices.scrubbedandactive,K_opt,TR);

disp('CAP analysis complete.')
