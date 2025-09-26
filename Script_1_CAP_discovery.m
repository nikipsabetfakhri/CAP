%% CAP discovery script (seed-free, TbCAPs, no GUI)

% Add TbCAPs + SPM to the path

% Mac
% addpath('/Users/npsabet/TbCAPs');   
% addpath('/Users/npsabet/spm12');
% addpath('/Users/npsabet/TbCAPs/control_pre_test_subj/');

% Dell 
addpath('/home/nsabet/TbCAPs');   
addpath('/home/nsabet/spm12');
addpath('/home/nsabet/RDoC_CAP/RDoC_Controls_Pre/');


%% 1. Loading the data files

% Set input folder
% data_dir = '/Users/npsabet/TbCAPs/control_pre_test_subj/';
data_dir = '/home/nsabet/RDoC_CAP/RDoC_Controls_Pre/';

% Get subject folders
subj_dirs = dir(fullfile(data_dir,'swau*'));
subj_dirs = subj_dirs([subj_dirs.isdir]);

% Inspect first subject's folder content
subj_dir = fullfile(data_dir, subj_dirs(1).name);
disp(subj_dir); 
ls(subj_dir);  % should list swauRDHC001_pre_00001.nii etc.

% Show what spm_select returns
files = spm_select('FPList', subj_dir, '.*\.nii$');
disp(files(1:min(size(files,1),5), :));  % first few paths

% Load subject header to use as reference space
% ref_img  = '/Users/npsabet/TbCAPs/control_pre_test_subj/swauRDHC003_pre/swauRDHC003_pre_00001.nii';
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
        [mask_path,mask_name,mask_ext] = fileparts(mask_in);
        resliced = fullfile(mask_path, ['r' mask_name mask_ext]);
        mask_out = fullfile(mask_path, [mask_name '_resliced' mask_ext]);
        movefile(resliced, mask_out, 'f');
        Vmask = spm_vol(mask_out);
else
    mask_out = mask_in; % no resampling needed
end

% Finalize mask
mask_vol  = spm_read_vols(Vmask);

% Mask: n_voxels x 1 logical vector
mask      = mask_vol > 0.5;  % binarize

% Header: the header (obtained by spm_vol) of one NIFTI file with proper
% data dimension and .mat information
brain_info = Vref;           % use subject header for clustering


%% 2. Build TC and FD (timecourses and framewise displacement)

% Time Course Data: cell array, each cell of size n_TP x n_masked_voxels
% Framewise displacement: a n_TP x n_subj matrix with framewise
% displacement information

nSubj = numel(subj_dirs);
TC    = cell(1,nSubj);              % cell array: each [timepoints × voxels]
FD    = zeros(240,nSubj);           % matrix: [timepoints × subjects] 
                                    % (replace 240 with your nTP if different)

for s = 1:nSubj
    subj_dir = fullfile(data_dir, subj_dirs(s).name);
    fprintf('\n[Subject %d/%d] %s\n', s, nSubj, subj_dirs(s).name);

    % ---------------------------
    % Load subject NIfTIs
    % ---------------------------
    files = spm_select('FPList', subj_dir, '^swau.*\.nii$');
    V = spm_vol(files);
    Y = spm_read_vols(V);               % dims: [X Y Z T]
    Y2D = reshape(Y, [], size(Y,4))';   % → [timepoints × voxels]
    TC{s} = Y2D(:, mask(:));            % apply mask

    % ---------------------------
    % Load motion file → FD(:,s)
    % ---------------------------
    motfile = dir(fullfile(subj_dir, 'rp_*.txt'));
    if isempty(motfile)
        warning('No motion file for %s. Using zeros.', subj_dirs(s).name);
        FD(:,s) = zeros(size(TC{s},1),1);
    else
        R = load(fullfile(subj_dir, motfile(1).name));   % [nTP × 6]
        dR = [zeros(1,6); diff(R)];
        trans = dR(:,1:3);               % translations in mm
        rot   = dR(:,4:6) * 50;          % rotations → mm (50 mm radius)
        FD(:,s) = sqrt(sum(trans.^2,2) + sum(rot.^2,2));
    end
end

fprintf('\nFinished loading TC and FD for %d subjects.\n', nSubj);

% Sanity checks
for s = 1:nSubj
    fprintf('Subject %d: TC = %d×%d, FD = %d frames, mean FD = %.3f mm\n', ...
        s, size(TC{s},1), size(TC{s},2), length(FD(:,s)), mean(FD(:,s)));
end

% Optional quick plots
% figure;
% for s = 1:nSubj
%     subplot(nSubj,1,s);
%     plot(FD(:,s),'k-'); 
%     ylabel(sprintf('FD subj %d',s));
%     if s==1, title('Framewise displacement (mm)'); end
% end
% xlabel('Frame');

% Clear workspace
clear Y Y2D


%% 3. Selecting the frames to analyse    

% Threshold of FD above which to scrub out the frame and also the t-1 and
% t+1 frames (if you want another scrubbing setting, directly edit the
% code)

Tmot = 0.5;    % scrubbing threshold

% Xon will contain the retained frames, and Indices will tag the time
% points associated to these frames, for each subject (it contains a
% subfield for retained frames and a subfield for scrubbed frames)
% Note that since this is a seed-free analysis, we do not need any
% seed-related information to be given as argument

[Xon,~,Indices] = CAP_find_activity_SeedFree(TC,FD,Tmot);

% Check # of frames retained per subject
nTotal = size(TC{1},1);  % should be 240
for s = 1:size(Indices.kept.active,2)
    nKept  = sum(Indices.kept.active(:,s));
    nScrub = sum(Indices.scrubbed(:,s));
    fprintf('Subj %d: kept %d / %d frames (%.1f%%), scrubbed %d\n', ...
        s, nKept, nTotal, 100*nKept/nTotal, nScrub);
end


%% 4. Consensus clustering (if wished to determine the optimum K)

dairyfile = fullfile('/home/nsabet/TbCAPs/SavedData', ...
    sprintf('CAP_log_%s.txt', datestr(now,'yyyymmdd_HHMM')));
diary(dairyfile);   % save console output
diary on

% Percentage of frames to use in each fold of consensus clustering
Pcc = 80; % If doing a small sample test, set this to 100
% so that it doesn't crash otherwise some of the denominators will be 0

% Number of folds we run consensus clustering for
N = 20;

% This specifies the range of values over which to perform consensus
% clustering: if you want to run parallel consensus clustering processes,
% you should feed in different ranges to each call of the function

K_range = 2:8;  
[Consensus] = CAP_ConsensusClustering(Xon, K_range, 'items', Pcc/100, N, 'correlation');

% Calculates the quality metrics
[~, Qual]   = ComputeClusteringQuality(Consensus,[]);

% Qual should be inspected to determine the best cluster number(s)

% Plot quick quality curve
% figure;
% plot(K_range, Qual,'-o','LineWidth',2);
% xlabel('K'); ylabel('Quality metric');
% title('Test run: consensus clustering quality');

save(fullfile('/home/nsabet/TbCAPs/SavedData', ...
    sprintf('Consensus_%s.mat', datestr(now,'yyyymmdd_HHMM'))), ...
    'Consensus','K_range','Pcc','N','-v7.3');

