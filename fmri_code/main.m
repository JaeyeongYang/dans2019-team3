% Before running this code, whole data of Tom et al. (2007) should be stored
% on the /data folder (relative to the root of the repository).
clear all

%% Set paths
datapath = fullfile(fileparts(pwd), 'data', 'fmriprep');
behavpath = fullfile(fileparts(pwd), 'behav_data');

%% First level analysis
run_l1_analysis(datapath, behavpath)

%% Second level analysis
run_l2_analysis(datapath)