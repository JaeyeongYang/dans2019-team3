function run_l1_analysis_accept_reject(datapath, behavpath, outdir, smoothing)

    if nargin < 3
        outdir = 'proc_1st_accept_reject_2';
    end

    if nargin < 4
        smoothing = [8 8 8];
    end

    subjdirs = dir(fullfile(datapath, 'sub-*'));
    dir_flags = [subjdirs.isdir];
    subjdirs = subjdirs(dir_flags);

    % Run 1st-level analysis in parallel.
    for i = 1:length(subjdirs)
        subjdir = subjdirs(i);
        subjid = subjdir.name;
        subjpath = fullfile(datapath, subjdir.name);
        funcpath = fullfile(subjpath, 'func');
        outpath = fullfile(subjpath, outdir);
        file_spmmat = fullfile(outpath, 'SPM.mat');

        % Create a directory where new data will be saved
        if ~exist(outpath)
            mkdir(outpath);
        end

        % delete SPM.mat file if it exists already
        if exist(file_spmmat)
            fprintf('\n SPM.mat exists in this directory. Overwriting SPM.mat file! \n\n')
            delete(file_spmmat)
        end

        % Unzip all the nii.gz files into nii files.
        gunzip(fullfile(funcpath, '*.gz'));

        %% Initialise SPM defaults
        spm('defaults', 'FMRI');
        spm_jobman('initcfg'); % SPM8 only

        batch = [];
        files_movement = {};
        files_movement_new = {};
        b = 1;

        for i = 1:3
            % Load confounding factors
            file_confounds = fullfile(funcpath, [subjid '_task-mixedgamblestask_run-0' num2str(i) '_bold_confounds.tsv']);
            [confdata, ~, ] = tsvread(file_confounds);

            % Remove the first row, 26-31 columns --> movement regressors
            R = confdata(2:end, (end - 5):end);
            file_movement = fullfile(funcpath, ['movement_regressors_for_epi_0' num2str(i) '.mat']);
            save(file_movement, 'R');
            files_movement = [files_movement, file_movement];

            % Remove the first row, 7-18 columns --> new movement regressors
            R = confdata(2:end, 7:18);
            file_movement_new = fullfile(funcpath, ['movement_regressors_for_epi_0' num2str(i) '_new.mat']);
            save(file_movement_new, 'R');
            files_movement_new = [files_movement_new, file_movement_new];

            % Scan files
            tmp_niis = dir(fullfile(funcpath, ['sub-*run-0' num2str(i) '*preproc.nii']));
            tmp_hdr = spm_vol(fullfile(funcpath, tmp_niis.name));

            for j = 1:size(tmp_hdr, 1)
                file_scans{j, 1} = fullfile(funcpath, [tmp_niis.name ',' num2str(j)]);
            end

            %batch{b}.spm.spatial.smooth.data = file_scans;
            %batch{b}.spm.spatial.smooth.fwhm = smoothing;
            %batch{b}.spm.spatial.smooth.dtype = 0;
            %batch{b}.spm.spatial.smooth.im = 0;
            %batch{b}.spm.spatial.smooth.prefix = 's';
            %b = b + 1;
        end

        %spm_jobman('run', batch);
        %disp('Smoothing is complete.')

        %% Model specification
        batch = [];

        batch{1}.spm.stats.fmri_spec.dir = {outpath};
        batch{1}.spm.stats.fmri_spec.timing.units = 'secs';
        batch{1}.spm.stats.fmri_spec.timing.RT = 2;
        batch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
        batch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

        for r = 1:3
            % Scan files
            tmp_niis = dir(fullfile(funcpath, ['ssub-*run-0' num2str(r) '*preproc.nii']));
            tmp_hdr = spm_vol(fullfile(funcpath, tmp_niis.name));

            for j = 1:size(tmp_hdr, 1)
                file_scans{j, 1} = fullfile(funcpath, [tmp_niis.name ',' num2str(j)]);
            end

            % Load run data
            runpath = fullfile(behavpath, [subjid '_task-mixedgamblestask_run-0' num2str(r) '_events.tsv']);
            [rundata, ~, ] = tsvread(runpath);

            batch{1}.spm.stats.fmri_spec.sess(r).scans = file_scans;
            batch{1}.spm.stats.fmri_spec.sess(r).cond.name = 'stimulus_onset';
            batch{1}.spm.stats.fmri_spec.sess(r).cond.onset = rundata(2:end, 1);
            batch{1}.spm.stats.fmri_spec.sess(r).cond.duration = 0;
            batch{1}.spm.stats.fmri_spec.sess(r).cond.tmod = 0;

            
            % PM1: parametric gain
            batch{1}.spm.stats.fmri_spec.sess(r).cond.pmod(1).name = 'gain';
            batch{1}.spm.stats.fmri_spec.sess(r).cond.pmod(1).param = rundata(2:end, 5);
            batch{1}.spm.stats.fmri_spec.sess(r).cond.pmod(1).poly = 1;

            % PM2: parametric loss
            batch{1}.spm.stats.fmri_spec.sess(r).cond.pmod(2).name = 'loss';
            batch{1}.spm.stats.fmri_spec.sess(r).cond.pmod(2).param = rundata(2:end, 3);
            batch{1}.spm.stats.fmri_spec.sess(r).cond.pmod(2).poly = 1;

            % PM3: indifference
            batch{1}.spm.stats.fmri_spec.sess(r).cond.pmod(3).name = 'indifference';
            batch{1}.spm.stats.fmri_spec.sess(r).cond.pmod(3).param = rundata(2:end, 4);
            batch{1}.spm.stats.fmri_spec.sess(r).cond.pmod(3).poly = 1;

            % PM4: accept
            batch{1}.spm.stats.fmri_spec.sess(r).cond.pmod(4).name = 'accept';
            batch{1}.spm.stats.fmri_spec.sess(r).cond.pmod(4).param = 1 * (rundata(2:end, 10) == 1);
            batch{1}.spm.stats.fmri_spec.sess(r).cond.pmod(4).poly = 1;

            % PM5: reject
            batch{1}.spm.stats.fmri_spec.sess(r).cond.pmod(5).name = 'reject';
            batch{1}.spm.stats.fmri_spec.sess(r).cond.pmod(5).param = 1 * (rundata(2:end, 10) == 0);
            batch{1}.spm.stats.fmri_spec.sess(r).cond.pmod(5).poly = 1;
            
            batch{1}.spm.stats.fmri_spec.sess(r).cond.orth = 0;
            batch{1}.spm.stats.fmri_spec.sess(r).multi = {''};
            batch{1}.spm.stats.fmri_spec.sess(r).regress = struct('name', {}, 'val', {});
            batch{1}.spm.stats.fmri_spec.sess(r).multi_reg = {files_movement_new{r}};
            batch{1}.spm.stats.fmri_spec.sess(r).hpf = 128;
        end

        batch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        batch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        batch{1}.spm.stats.fmri_spec.volt = 1;
        batch{1}.spm.stats.fmri_spec.global = 'None';
        batch{1}.spm.stats.fmri_spec.mthresh = 0.2;
        batch{1}.spm.stats.fmri_spec.mask = {''};
        batch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

        %% Categorical model estimation
        batch{2}.spm.stats.fmri_est.spmmat = {file_spmmat};
        batch{2}.spm.stats.fmri_est.method.Classical = 1;

        %% Create contrasts
        % parametric modulation of gain & loss
        batch{3}.spm.stats.con.spmmat = {file_spmmat};

        con_acc_rej = [0 0 0 0 1  -1 0 0 0 0  0 0 0 0 0  0 0 0];
        con_gain = [0 1 0 0 0  0 0 0 0 0  0 0 0 0 0  0 0 0];
        con_loss = [0 0 1 0 0  0 0 0 0 0  0 0 0 0 0  0 0 0];
        con_gain_loss = con_gain - con_loss;

        batch{3}.spm.stats.con.consess{1}.tcon.name = 'accept_reject_PM';
        batch{3}.spm.stats.con.consess{1}.tcon.convec = [con_acc_rej con_acc_rej con_acc_rej];
        batch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

        batch{3}.spm.stats.con.consess{2}.tcon.name = 'gain_PM';
        batch{3}.spm.stats.con.consess{2}.tcon.convec = [con_gain con_gain con_gain];
        batch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

        batch{3}.spm.stats.con.consess{3}.tcon.name = 'loss_PM';
        batch{3}.spm.stats.con.consess{3}.tcon.convec = [con_gain con_gain con_gain];
        batch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';

        batch{3}.spm.stats.con.consess{4}.tcon.name = 'gain_loss_PM';
        batch{3}.spm.stats.con.consess{4}.tcon.convec = [con_gain_loss con_gain_loss con_gain_loss];
        batch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';

        batch{3}.spm.stats.con.delete = 0;

        spm_jobman('run', batch);
        disp(['Job is done for ' subjid '.']);
    end
