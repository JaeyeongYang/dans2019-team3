function run_l2_analysis_accept_reject(datapath, outdir)

    if nargin < 2
        outdir = 'proc_2nd';
    end

    subjdirs = dir(fullfile(datapath, 'sub-*'));
    dir_flags = [subjdirs.isdir];
    subjdirs = subjdirs(dir_flags);

    outpath = fullfile(datapath, outdir);
    cond_name = 'stimulus_onset';

    % covars = {'age', 'male'};
    covars = {};

    all_conds = {'gain', 'loss'};
    subjids = {subjdirs.name};
    n_subjs = length(subjids);

    %% Initialise SPM defaults
    spm('Defaults', 'fMRI');
    spm_jobman('initcfg'); % SPM8 only

    for i = 1:length(all_conds)
        cond = all_conds{i}
        dir_cond = fullfile(outpath, [cond '_N' num2str(n_subjs)])
        file_spmmat = fullfile(dir_cond, 'SPM.mat')

        if ~exist(dir_cond)
            mkdir(dir_cond);
        end

        if exist(file_spmmat)
            fprintf('\n SPM.mat exists in this directory. Overwriting SPM.mat file! \n\n');
            delete(file_spmmat);
        end

        % Contrast number
        switch cond
            case 'gain'
                contrast_num = '0001';
            case 'loss'
                contrast_num = '0002';
        end

        %% specification
        batch = [];
        scanFiles = strcat(datapath, '/', subjids', '/proc_1st/con_', contrast_num, '.nii,1') % if using subjids

        batch{1}.spm.stats.factorial_design.dir = {dir_cond};
        batch{1}.spm.stats.factorial_design.des.t1.scans = scanFiles;
        batch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        batch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        batch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        batch{1}.spm.stats.factorial_design.masking.im = 1;
        batch{1}.spm.stats.factorial_design.masking.em = {''};
        batch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
        batch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        batch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

        spm_jobman('run', batch)
        disp('2nd Level model is specified');

        %% Estimation
        batch = [];
        batch{1}.spm.stats.fmri_est.spmmat = {file_spmmat};
        batch{1}.spm.stats.fmri_est.method.Classical = 1;
        spm_jobman('run', batch)

        %% T-contrast (one-step t-test)
        batch = [];
        batch{1}.spm.stats.con.spmmat = {file_spmmat};
        batch{1}.spm.stats.con.consess{1}.tcon.name = [cond '_N' num2str(n_subjs)];
        batch{1}.spm.stats.con.consess{1}.tcon.convec = 1;
        batch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        batch{1}.spm.stats.con.delete = 1;
        spm_jobman('run', batch)

        disp([[cond '_N' num2str(n_subjs)] ' contrast is created'])
        disp('All done')

    end
