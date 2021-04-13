%-----------------------------------------------------------------------
% Job created by Gianpaolo & Davide & Luca Bondi
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
% 12/04/2021
%-----------------------------------------------------------------------
%%

%clear all
%spm_jobman('initcfg');
spm('defaults', 'FMRI');
global defaults;

BasePath = ('C:\Users\CNPL\Desktop\switch_logs\pipeline');
% Warn user if there is no such folder.
if ~exist(BasePath, 'dir')
  message = sprintf('This folder does not exist:\n%s', BasePath);
  uiwait(errordlg(message));
  return;
end
% Get a list of all files, including folders.
cd (BasePath)
DirList  = dir(BasePath);
% Extract only the folders, not regular files.
DirList  = DirList([DirList.isdir]);  % Folders only
% Get rid of first two folders: dot and dot dot.
DirList = DirList(3:end);
% Warn user if there are no subfolders.
if isempty(DirList)
  message = sprintf('This folder does not contain any subfolders:\n%s', BasePath);
  uiwait(errordlg(message));
  return;
end
% Count the number of subfolders.
numberOfFolders = numel(DirList);
% Loop over all subfolders, processing each one.

%-------------------------------------------------------------------------------
%for k = 1: numberOfFolders %1
k = 2 
thisDir = fullfile(BasePath, DirList(k).name);
  subdirname = dir(thisDir);
  subdirname = subdirname(3:end);
  ENG_fold = fullfile(BasePath, DirList(k).name, subdirname(1).name); %Switch_fold 
  ITA_fold = fullfile(BasePath, DirList(k).name, subdirname(2).name); %stroop_fold
  Switch_fold = fullfile(BasePath, DirList(k).name, subdirname(3).name);
  T1w_fold = fullfile(BasePath, DirList(k).name, subdirname(4).name);
%   Volume_fold = fullfile(BasePath, DirList(k).name, 'Volume');
%   Volume_switch = fullfile(Volume_fold, 'Switch');
%   Volume_Ita = fullfile(Volume_fold, 'Ita');
%   Volume_Eng = fullfile(Volume_fold, 'Eng');
  
  % %%%% VOLUME FOLDERS%%%%%
% mkdir(fullfile(BasePath, DirList(k).name, 'Volume'));
 Volume_fold = fullfile(BasePath, DirList(k).name, 'Volume');
% mkdir(fullfile(Volume_fold, 'Switch'));
 Volume_switch = fullfile(Volume_fold, 'Switch');
% mkdir(fullfile(Volume_fold, 'Ita'));
 Volume_Ita = fullfile(Volume_fold, 'Ita');
% mkdir(fullfile(Volume_fold, 'Eng'));
 Volume_Eng = fullfile(Volume_fold, 'Eng');


  fprintf('Processing subject %d of %d: %s\n', ...
  k, numberOfFolders, thisDir);
  pathparts = strsplit(thisDir,filesep);
  %subname = (pathparts{1,4}); for H: hard disk
  subname = (pathparts{1,7});
  
  T1 = dir(fullfile(T1w_fold, '*_3D_T1_0.7_ax_*.nii'));
  nameselect_0 = strcat(T1.name, ',1');
  
  %% image expansion
  switch_run_1 = dir(fullfile(Switch_fold, '*_fMRI_Switch_A*.nii')); 
  nameselect_1 = strcat(switch_run_1.name, ',1');
  matlabbatch{1}.spm.util.exp_frames.files = {(fullfile(Switch_fold, nameselect_1))};
  matlabbatch{1}.spm.util.exp_frames.frames = Inf;
  switch_run_2 = dir(fullfile(Switch_fold, '*_fMRI_Switch_B*.nii'));
  nameselect_2 = strcat(switch_run_2.name, ',1');
  matlabbatch{2}.spm.util.exp_frames.files = {(fullfile(Switch_fold, nameselect_2))};
  matlabbatch{2}.spm.util.exp_frames.frames = Inf;

%%
%%slice time correction
matlabbatch{3}.spm.temporal.st.scans{1}(1) = cfg_dep('Expand image frames: Expanded filename list.', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{3}.spm.temporal.st.scans{2}(1) = cfg_dep('Expand image frames: Expanded filename list.', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{3}.spm.temporal.st.nslices = 35;
matlabbatch{3}.spm.temporal.st.tr = 2;
matlabbatch{3}.spm.temporal.st.ta = 1.94285714285714; 
matlabbatch{3}.spm.temporal.st.so = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34];
matlabbatch{3}.spm.temporal.st.refslice = 1;
matlabbatch{3}.spm.temporal.st.prefix = 'a';

%%% realign and unwarp
matlabbatch{4}.spm.spatial.realignunwarp.data(1).scans(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 1)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{4}.spm.spatial.realignunwarp.data(1).pmscan = '';
matlabbatch{4}.spm.spatial.realignunwarp.data(2).scans(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 2)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{2}, '.','files'));
matlabbatch{4}.spm.spatial.realignunwarp.data(2).pmscan = '';
matlabbatch{4}.spm.spatial.realignunwarp.eoptions.quality = 0.9;
matlabbatch{4}.spm.spatial.realignunwarp.eoptions.sep = 4;
matlabbatch{4}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
matlabbatch{4}.spm.spatial.realignunwarp.eoptions.rtm = 0;
matlabbatch{4}.spm.spatial.realignunwarp.eoptions.einterp = 2;
matlabbatch{4}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
matlabbatch{4}.spm.spatial.realignunwarp.eoptions.weight = '';
matlabbatch{4}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
matlabbatch{4}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
matlabbatch{4}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
matlabbatch{4}.spm.spatial.realignunwarp.uweoptions.jm = 0;
matlabbatch{4}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
matlabbatch{4}.spm.spatial.realignunwarp.uweoptions.sot = [];
matlabbatch{4}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
matlabbatch{4}.spm.spatial.realignunwarp.uweoptions.rem = 1;
matlabbatch{4}.spm.spatial.realignunwarp.uweoptions.noi = 5;
matlabbatch{4}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
matlabbatch{4}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
matlabbatch{4}.spm.spatial.realignunwarp.uwroptions.rinterp = 4;
matlabbatch{4}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
matlabbatch{4}.spm.spatial.realignunwarp.uwroptions.mask = 1;
matlabbatch{4}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';
%%% check FIACH R package

%%%segmentation and surface estimation
matlabbatch{5}.spm.tools.cat.estwrite.data = {(fullfile(T1w_fold, nameselect_0))};
matlabbatch{5}.spm.tools.cat.estwrite.data_wmh = {''};
matlabbatch{5}.spm.tools.cat.estwrite.nproc = 2;
matlabbatch{5}.spm.tools.cat.estwrite.useprior = '';
matlabbatch{5}.spm.tools.cat.estwrite.opts.tpm = {'C:\Users\CNPL\Documents\spm12\tpm\TPM.nii'};
matlabbatch{5}.spm.tools.cat.estwrite.opts.affreg = 'mni';
matlabbatch{5}.spm.tools.cat.estwrite.opts.biasstr = 0.5;
matlabbatch{5}.spm.tools.cat.estwrite.opts.accstr = 0.5;
matlabbatch{5}.spm.tools.cat.estwrite.extopts.segmentation.APP = 1070;
matlabbatch{5}.spm.tools.cat.estwrite.extopts.segmentation.NCstr = -Inf;
matlabbatch{5}.spm.tools.cat.estwrite.extopts.segmentation.spm_kamap = 0;
matlabbatch{5}.spm.tools.cat.estwrite.extopts.segmentation.LASstr = 0.5;
matlabbatch{5}.spm.tools.cat.estwrite.extopts.segmentation.gcutstr = 2;
matlabbatch{5}.spm.tools.cat.estwrite.extopts.segmentation.cleanupstr = 0.5;
matlabbatch{5}.spm.tools.cat.estwrite.extopts.segmentation.BVCstr = 0.5;
matlabbatch{5}.spm.tools.cat.estwrite.extopts.segmentation.WMHC = 1;
matlabbatch{5}.spm.tools.cat.estwrite.extopts.segmentation.SLC = 0;
matlabbatch{5}.spm.tools.cat.estwrite.extopts.segmentation.mrf = 1;
matlabbatch{5}.spm.tools.cat.estwrite.extopts.segmentation.restypes.optimal = [1 0.1];
matlabbatch{5}.spm.tools.cat.estwrite.extopts.registration.shooting.shootingtpm = {'C:\Users\CNPL\Documents\spm12\toolbox\cat12\templates_volumes\Template_0_IXI555_MNI152_GS.nii'};
matlabbatch{5}.spm.tools.cat.estwrite.extopts.registration.shooting.regstr = 0.5;
matlabbatch{5}.spm.tools.cat.estwrite.extopts.vox = 1.5;
matlabbatch{5}.spm.tools.cat.estwrite.extopts.surface.pbtres = 0.5;
matlabbatch{5}.spm.tools.cat.estwrite.extopts.surface.pbtmethod = 'pbt2x';
matlabbatch{5}.spm.tools.cat.estwrite.extopts.surface.pbtlas = 0;
matlabbatch{5}.spm.tools.cat.estwrite.extopts.surface.collcorr = 0;
matlabbatch{5}.spm.tools.cat.estwrite.extopts.surface.reduce_mesh = 1;
matlabbatch{5}.spm.tools.cat.estwrite.extopts.surface.vdist = 1.33333333333333;
matlabbatch{5}.spm.tools.cat.estwrite.extopts.surface.scale_cortex = 0.7;
matlabbatch{5}.spm.tools.cat.estwrite.extopts.surface.add_parahipp = 0.1;
matlabbatch{5}.spm.tools.cat.estwrite.extopts.surface.close_parahipp = 0;
matlabbatch{5}.spm.tools.cat.estwrite.extopts.admin.experimental = 0;
matlabbatch{5}.spm.tools.cat.estwrite.extopts.admin.new_release = 0;
matlabbatch{5}.spm.tools.cat.estwrite.extopts.admin.lazy = 0;
matlabbatch{5}.spm.tools.cat.estwrite.extopts.admin.ignoreErrors = 1;
matlabbatch{5}.spm.tools.cat.estwrite.extopts.admin.verb = 2;
matlabbatch{5}.spm.tools.cat.estwrite.extopts.admin.print = 2;
matlabbatch{5}.spm.tools.cat.estwrite.output.surface = 1;
matlabbatch{5}.spm.tools.cat.estwrite.output.surf_measures = 1;
matlabbatch{5}.spm.tools.cat.estwrite.output.ROImenu.atlases.neuromorphometrics = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.ROImenu.atlases.lpba40 = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.ROImenu.atlases.cobra = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.ROImenu.atlases.hammers = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.ROImenu.atlases.ibsr = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.ROImenu.atlases.aal3 = 1;
matlabbatch{5}.spm.tools.cat.estwrite.output.ROImenu.atlases.mori = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.ROImenu.atlases.anatomy = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.ROImenu.atlases.julichbrain = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.ROImenu.atlases.Schaefer2018_100Parcels_17Networks_order = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.ROImenu.atlases.Schaefer2018_200Parcels_17Networks_order = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.ROImenu.atlases.Schaefer2018_400Parcels_17Networks_order = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.ROImenu.atlases.Schaefer2018_600Parcels_17Networks_order = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.ROImenu.atlases.ownatlas = {''};
matlabbatch{5}.spm.tools.cat.estwrite.output.GM.native = 1;
matlabbatch{5}.spm.tools.cat.estwrite.output.GM.warped = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.GM.mod = 1;
matlabbatch{5}.spm.tools.cat.estwrite.output.GM.dartel = 2;
matlabbatch{5}.spm.tools.cat.estwrite.output.WM.native = 1;
matlabbatch{5}.spm.tools.cat.estwrite.output.WM.warped = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.WM.mod = 1;
matlabbatch{5}.spm.tools.cat.estwrite.output.WM.dartel = 2;
matlabbatch{5}.spm.tools.cat.estwrite.output.CSF.native = 1;
matlabbatch{5}.spm.tools.cat.estwrite.output.CSF.warped = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.CSF.mod = 1;
matlabbatch{5}.spm.tools.cat.estwrite.output.CSF.dartel = 2;
matlabbatch{5}.spm.tools.cat.estwrite.output.ct.native = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.ct.warped = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.ct.dartel = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.pp.native = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.pp.warped = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.pp.dartel = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.WMH.native = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.WMH.warped = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.WMH.mod = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.WMH.dartel = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.SL.native = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.SL.warped = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.SL.mod = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.SL.dartel = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.TPMC.native = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.TPMC.warped = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.TPMC.mod = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.TPMC.dartel = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.atlas.native = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.atlas.warped = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.atlas.dartel = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.label.native = 1;
matlabbatch{5}.spm.tools.cat.estwrite.output.label.warped = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.label.dartel = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.labelnative = 1;
matlabbatch{5}.spm.tools.cat.estwrite.output.bias.native = 1;
matlabbatch{5}.spm.tools.cat.estwrite.output.bias.warped = 1;
matlabbatch{5}.spm.tools.cat.estwrite.output.bias.dartel = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.las.native = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.las.warped = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.las.dartel = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.jacobianwarped = 0;
matlabbatch{5}.spm.tools.cat.estwrite.output.warps = [1 1];
matlabbatch{5}.spm.tools.cat.estwrite.output.rmat = 0;

%%skull stripping
matlabbatch{6}.spm.util.imcalc.input(1) = cfg_dep('CAT12: Segmentation (current release): Native Bias Corr. Image', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','biascorr', '()',{':'}));
matlabbatch{6}.spm.util.imcalc.input(2) = cfg_dep('CAT12: Segmentation (current release): p1 Image', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','p', '()',{':'}));
matlabbatch{6}.spm.util.imcalc.input(3) = cfg_dep('CAT12: Segmentation (current release): p2 Image', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','p', '()',{':'}));
matlabbatch{6}.spm.util.imcalc.input(4) = cfg_dep('CAT12: Segmentation (current release): p3 Image', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{3}, '.','p', '()',{':'}));
matlabbatch{6}.spm.util.imcalc.output = 'Brain';
matlabbatch{6}.spm.util.imcalc.outdir = {T1w_fold};
matlabbatch{6}.spm.util.imcalc.expression = '(i2 + i3 + i4) .* i1';
matlabbatch{6}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{6}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{6}.spm.util.imcalc.options.mask = 0;
matlabbatch{6}.spm.util.imcalc.options.interp = 1;
matlabbatch{6}.spm.util.imcalc.options.dtype = 4;

%%% coregistration
matlabbatch{7}.spm.spatial.coreg.estimate.ref(1) = cfg_dep('Image Calculator: ImCalc Computed Image: Brain', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{7}.spm.spatial.coreg.estimate.source(1) = cfg_dep('Realign & Unwarp: Unwarped Mean Image', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','meanuwr'));
matlabbatch{7}.spm.spatial.coreg.estimate.other(1) = cfg_dep('Realign & Unwarp: Unwarped Images (Sess 1)', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','uwrfiles'));
matlabbatch{7}.spm.spatial.coreg.estimate.other(2) = cfg_dep('Realign & Unwarp: Unwarped Images (Sess 2)', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{2}, '.','uwrfiles'));
matlabbatch{7}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{7}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{7}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{7}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];


%%%% SWITCH - SURFACE %%%%

%%% FIRST LEVEL Switching A
matlabbatch{8}.spm.stats.fmri_spec.dir = {Switch_fold};
matlabbatch{8}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{8}.spm.stats.fmri_spec.timing.RT = 2;
matlabbatch{8}.spm.stats.fmri_spec.timing.fmri_t = 35;
matlabbatch{8}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
matlabbatch{8}.spm.stats.fmri_spec.sess(1).scans(1) = cfg_dep('Realign & Unwarp: Unwarped Images (Sess 1)', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','uwrfiles'));
%%
matlabbatch{8}.spm.stats.fmri_spec.sess(1).cond(1).name = 'ITA_switch';
ITA_switch_A = dir(fullfile(Switch_fold, '*_switch_A_ita_switch.txt'));
ITA_switch_A = load(fullfile(Switch_fold,ITA_switch_A.name));
matlabbatch{8}.spm.stats.fmri_spec.sess(1).cond(1).onset = ITA_switch_A;
matlabbatch{8}.spm.stats.fmri_spec.sess(1).cond(1).duration = 0;
matlabbatch{8}.spm.stats.fmri_spec.sess(1).cond(1).tmod = 0;
matlabbatch{8}.spm.stats.fmri_spec.sess(1).cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{8}.spm.stats.fmri_spec.sess(1).cond(1).orth = 1;
%%
matlabbatch{8}.spm.stats.fmri_spec.sess(1).cond(2).name = 'ENG_switch';
ENG_switch_A = dir(fullfile(Switch_fold, '*_switch_A_eng_switch.txt'));
ENG_switch_A = load(fullfile(Switch_fold,ENG_switch_A.name));
matlabbatch{8}.spm.stats.fmri_spec.sess(1).cond(2).onset = ENG_switch_A;
matlabbatch{8}.spm.stats.fmri_spec.sess(1).cond(2).duration = 0;
matlabbatch{8}.spm.stats.fmri_spec.sess(1).cond(2).tmod = 0;
matlabbatch{8}.spm.stats.fmri_spec.sess(1).cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{8}.spm.stats.fmri_spec.sess(1).cond(2).orth = 1;
%%
matlabbatch{8}.spm.stats.fmri_spec.sess(1).cond(3).name = 'ITA_stay';
ITA_stay_A = dir(fullfile(Switch_fold, '*_switch_A_ita_stay.txt'));
ITA_stay_A = load(fullfile(Switch_fold,ITA_stay_A.name));
matlabbatch{8}.spm.stats.fmri_spec.sess(1).cond(3).onset = ITA_stay_A;
matlabbatch{8}.spm.stats.fmri_spec.sess(1).cond(3).duration = 0;
matlabbatch{8}.spm.stats.fmri_spec.sess(1).cond(3).tmod = 0;
matlabbatch{8}.spm.stats.fmri_spec.sess(1).cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{8}.spm.stats.fmri_spec.sess(1).cond(3).orth = 1;

%%
matlabbatch{8}.spm.stats.fmri_spec.sess(1).cond(4).name = 'ENG_stay';
ENG_stay_A = dir(fullfile(Switch_fold, '*_switch_A_eng_stay.txt'));
ENG_stay_A = load(fullfile(Switch_fold,ENG_stay_A.name));
matlabbatch{8}.spm.stats.fmri_spec.sess(1).cond(4).onset = ENG_stay_A;
matlabbatch{8}.spm.stats.fmri_spec.sess(1).cond(4).duration = 0;
matlabbatch{8}.spm.stats.fmri_spec.sess(1).cond(4).tmod = 0;
matlabbatch{8}.spm.stats.fmri_spec.sess(1).cond(4).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{8}.spm.stats.fmri_spec.sess(1).cond(4).orth = 1;


%%
matlabbatch{8}.spm.stats.fmri_spec.sess(1).multi = {''};
matlabbatch{8}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});
matlabbatch{8}.spm.stats.fmri_spec.sess(1).multi_reg = cfg_dep('Realign & Unwarp: Realignment Param File (Sess 1)', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','rpfile'));;
matlabbatch{8}.spm.stats.fmri_spec.sess(1).hpf = 128;
%%


%%% FIRST LEVEL Switching B
matlabbatch{8}.spm.stats.fmri_spec.sess(2).scans(1) = cfg_dep('Realign & Unwarp: Unwarped Images (Sess 2)', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{2}, '.','uwrfiles'));

%%
matlabbatch{8}.spm.stats.fmri_spec.sess(2).cond(1).name = 'ITA_switch';
ITA_switch_B = dir(fullfile(Switch_fold, '*_switch_B_ita_switch.txt'));
ITA_switch_B = load(fullfile(Switch_fold,ITA_switch_B.name));
matlabbatch{8}.spm.stats.fmri_spec.sess(2).cond(1).onset = ITA_switch_B;
matlabbatch{8}.spm.stats.fmri_spec.sess(2).cond(1).duration = 0;
matlabbatch{8}.spm.stats.fmri_spec.sess(2).cond(1).tmod = 0;
matlabbatch{8}.spm.stats.fmri_spec.sess(2).cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{8}.spm.stats.fmri_spec.sess(2).cond(1).orth = 1;
%%
matlabbatch{8}.spm.stats.fmri_spec.sess(2).cond(2).name = 'ENG_switch';
ENG_switch_B = dir(fullfile(Switch_fold, '*_switch_B_eng_switch.txt'));
ENG_switch_B = load(fullfile(Switch_fold,ENG_switch_B.name));
matlabbatch{8}.spm.stats.fmri_spec.sess(2).cond(2).onset = ENG_switch_B;
matlabbatch{8}.spm.stats.fmri_spec.sess(2).cond(2).duration = 0;
matlabbatch{8}.spm.stats.fmri_spec.sess(2).cond(2).tmod = 0;
matlabbatch{8}.spm.stats.fmri_spec.sess(2).cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{8}.spm.stats.fmri_spec.sess(2).cond(2).orth = 1;
%%
matlabbatch{8}.spm.stats.fmri_spec.sess(2).cond(3).name = 'ITA_stay';
ITA_stay_B = dir(fullfile(Switch_fold, '*_switch_B_ita_stay.txt'));
ITA_stay_B  = load(fullfile(Switch_fold,ITA_stay_B.name));
matlabbatch{8}.spm.stats.fmri_spec.sess(2).cond(3).onset = ITA_stay_B ;
matlabbatch{8}.spm.stats.fmri_spec.sess(2).cond(3).duration = 0;
matlabbatch{8}.spm.stats.fmri_spec.sess(2).cond(3).tmod = 0;
matlabbatch{8}.spm.stats.fmri_spec.sess(2).cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{8}.spm.stats.fmri_spec.sess(2).cond(3).orth = 1;
%%
matlabbatch{8}.spm.stats.fmri_spec.sess(2).cond(4).name = 'ENG_stay';
ENG_stay_B = dir(fullfile(Switch_fold, '*_switch_B_eng_stay.txt'));
ENG_stay_B = load(fullfile(Switch_fold,ENG_stay_B.name));
matlabbatch{8}.spm.stats.fmri_spec.sess(2).cond(4).onset = ENG_stay_B ;
matlabbatch{8}.spm.stats.fmri_spec.sess(2).cond(4).duration = 0;
matlabbatch{8}.spm.stats.fmri_spec.sess(2).cond(4).tmod = 0;
matlabbatch{8}.spm.stats.fmri_spec.sess(2).cond(4).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{8}.spm.stats.fmri_spec.sess(2).cond(4).orth = 1;
%%
matlabbatch{8}.spm.stats.fmri_spec.sess(2).multi = {''};
matlabbatch{8}.spm.stats.fmri_spec.sess(2).regress = struct('name', {}, 'val', {});
matlabbatch{8}.spm.stats.fmri_spec.sess(2).multi_reg = cfg_dep('Realign & Unwarp: Realignment Param File (Sess 2)', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{2}, '.','rpfile'));
matlabbatch{8}.spm.stats.fmri_spec.sess(2).hpf = 128;

%% FIRST LEVEL ESTIMATION and CONTRASTS
matlabbatch{8}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{8}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{8}.spm.stats.fmri_spec.volt = 1;
matlabbatch{8}.spm.stats.fmri_spec.global = 'None';
matlabbatch{8}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{8}.spm.stats.fmri_spec.mask = {''};
matlabbatch{8}.spm.stats.fmri_spec.cvi = 'AR(1)';
%%
matlabbatch{9}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{8}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{9}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{9}.spm.stats.fmri_est.method.Classical = 1;
%%
matlabbatch{10}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{9}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{10}.spm.stats.con.consess{1}.tcon.name = 'Main_ITA_switch';
matlabbatch{10}.spm.stats.con.consess{1}.tcon.weights = [1 0 0 0];
matlabbatch{10}.spm.stats.con.consess{1}.tcon.sessrep = 'repl';
matlabbatch{10}.spm.stats.con.consess{2}.tcon.name = 'Main_ENG_switch';
matlabbatch{10}.spm.stats.con.consess{2}.tcon.weights = [0 1 0 0];
matlabbatch{10}.spm.stats.con.consess{2}.tcon.sessrep = 'repl';
matlabbatch{10}.spm.stats.con.consess{3}.tcon.name = 'Main_ITA_stay';
matlabbatch{10}.spm.stats.con.consess{3}.tcon.weights = [0 0 1 0];
matlabbatch{10}.spm.stats.con.consess{3}.tcon.sessrep = 'repl';
matlabbatch{10}.spm.stats.con.consess{4}.tcon.name = 'Main_ENG_stay';
matlabbatch{10}.spm.stats.con.consess{4}.tcon.weights = [0 0 0 1];
matlabbatch{10}.spm.stats.con.consess{4}.tcon.sessrep = 'repl';
matlabbatch{10}.spm.stats.con.consess{5}.tcon.name = 'Main_switch';
matlabbatch{10}.spm.stats.con.consess{5}.tcon.weights = [1 1 0 0];
matlabbatch{10}.spm.stats.con.consess{5}.tcon.sessrep = 'repl';
matlabbatch{10}.spm.stats.con.consess{6}.tcon.name = 'Main_stay';
matlabbatch{10}.spm.stats.con.consess{6}.tcon.weights = [0 0 1 1];
matlabbatch{10}.spm.stats.con.consess{6}.tcon.sessrep = 'repl';
matlabbatch{10}.spm.stats.con.consess{7}.tcon.name = 'Main_ITA';
matlabbatch{10}.spm.stats.con.consess{7}.tcon.weights = [1 0 1 0];
matlabbatch{10}.spm.stats.con.consess{7}.tcon.sessrep = 'repl';
matlabbatch{10}.spm.stats.con.consess{8}.tcon.name = 'Main_ENG';
matlabbatch{10}.spm.stats.con.consess{8}.tcon.weights = [0 1 0 1];
matlabbatch{10}.spm.stats.con.consess{8}.tcon.sessrep = 'repl';
matlabbatch{10}.spm.stats.con.consess{9}.tcon.name = 'switch_gt_stay';
matlabbatch{10}.spm.stats.con.consess{9}.tcon.weights = [1 1 -1 -1];
matlabbatch{10}.spm.stats.con.consess{9}.tcon.sessrep = 'repl';
matlabbatch{10}.spm.stats.con.consess{10}.tcon.name = 'stay_gt_switch';
matlabbatch{10}.spm.stats.con.consess{10}.tcon.weights = [-1 -1 1 1];
matlabbatch{10}.spm.stats.con.consess{10}.tcon.sessrep = 'repl';
matlabbatch{10}.spm.stats.con.consess{11}.tcon.name = 'ITA_gt_ENG';
matlabbatch{10}.spm.stats.con.consess{11}.tcon.weights = [1 -1 1 -1];
matlabbatch{10}.spm.stats.con.consess{11}.tcon.sessrep = 'repl';
matlabbatch{10}.spm.stats.con.consess{12}.tcon.name = 'ENG_gt_ITA';
matlabbatch{10}.spm.stats.con.consess{12}.tcon.weights = [-1 1 -1 1];
matlabbatch{10}.spm.stats.con.consess{12}.tcon.sessrep = 'repl';
matlabbatch{10}.spm.stats.con.consess{13}.tcon.name = 'ITA_switch_gt_ENG_switch';
matlabbatch{10}.spm.stats.con.consess{13}.tcon.weights = [1 -1 0 0];
matlabbatch{10}.spm.stats.con.consess{13}.tcon.sessrep = 'repl';
matlabbatch{10}.spm.stats.con.consess{14}.tcon.name = 'ENG_switch_gt_ITA_switch';
matlabbatch{10}.spm.stats.con.consess{14}.tcon.weights = [-1 1 0 0];
matlabbatch{10}.spm.stats.con.consess{14}.tcon.sessrep = 'repl';
matlabbatch{10}.spm.stats.con.consess{15}.tcon.name = 'Inter_1';
matlabbatch{10}.spm.stats.con.consess{15}.tcon.weights = [-1 1 1 -1];
matlabbatch{10}.spm.stats.con.consess{15}.tcon.sessrep = 'repl';
matlabbatch{10}.spm.stats.con.consess{16}.tcon.name = 'Inter_2';
matlabbatch{10}.spm.stats.con.consess{16}.tcon.weights = [1 -1 -1 1];
matlabbatch{10}.spm.stats.con.consess{16}.tcon.sessrep = 'repl';
matlabbatch{10}.spm.stats.con.delete = 0;

%% MAPPING TO SURFACE
matlabbatch{11}.spm.tools.cat.stools.vol2surf.data_vol(1) = cfg_dep('Contrast Manager: All Con Images', substruct('.','val', '{}',{10}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','con'));
matlabbatch{11}.spm.tools.cat.stools.vol2surf.data_mesh_lh(1) = cfg_dep('CAT12: Segmentation: Left Central Surface', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','lhcentral', '()',{':'}));
matlabbatch{11}.spm.tools.cat.stools.vol2surf.sample = {'maxabs'};
matlabbatch{11}.spm.tools.cat.stools.vol2surf.interp = {'linear'};
matlabbatch{11}.spm.tools.cat.stools.vol2surf.datafieldname = 'Switch';
matlabbatch{11}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.class = 'GM';
matlabbatch{11}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.startpoint = -0.5;
matlabbatch{11}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.steps = 7;
matlabbatch{11}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.endpoint = 0.5;
matlabbatch{12}.spm.tools.cat.stools.surfresamp.sample{1}.data_surf(1) = cfg_dep('Map Volume (Native Space) to Individual Surface: Left mapped values', substruct('.','val', '{}',{11}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','lh'));
matlabbatch{12}.spm.tools.cat.stools.surfresamp.merge_hemi = 1;
matlabbatch{12}.spm.tools.cat.stools.surfresamp.mesh32k = 1;
matlabbatch{12}.spm.tools.cat.stools.surfresamp.fwhm_surf = 3;
matlabbatch{12}.spm.tools.cat.stools.surfresamp.lazy = 0;
matlabbatch{12}.spm.tools.cat.stools.surfresamp.nproc = 4;
matlabbatch{13}.spm.tools.cat.stools.surfresamp.sample{1}.data_surf(1) = cfg_dep('Map Volume (Native Space) to Individual Surface: Left mapped values', substruct('.','val', '{}',{11}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','lh'));
matlabbatch{13}.spm.tools.cat.stools.surfresamp.merge_hemi = 0;
matlabbatch{13}.spm.tools.cat.stools.surfresamp.mesh32k = 1;
matlabbatch{13}.spm.tools.cat.stools.surfresamp.fwhm_surf = 3;
matlabbatch{13}.spm.tools.cat.stools.surfresamp.lazy = 0;
matlabbatch{13}.spm.tools.cat.stools.surfresamp.nproc = 4;

%%ITA_mono%%
  ITA_mono = dir(fullfile(ITA_fold, '*fMRI_Switch_mono_IT_*.nii')); 
  ITAselect = strcat(ITA_mono.name, ',1');
%% image expansion
matlabbatch{14}.spm.util.exp_frames.files = {(fullfile(ITA_fold, ITAselect))};
matlabbatch{14}.spm.util.exp_frames.frames = Inf;
%% slice time correction
matlabbatch{15}.spm.temporal.st.scans{1}(1) = cfg_dep('Expand image frames: Expanded filename list.', substruct('.','val', '{}',{14}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{15}.spm.temporal.st.nslices = 35;
matlabbatch{15}.spm.temporal.st.tr = 2;
matlabbatch{15}.spm.temporal.st.ta = 1.94285714285714;
matlabbatch{15}.spm.temporal.st.so = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34];
matlabbatch{15}.spm.temporal.st.refslice = 1;
matlabbatch{15}.spm.temporal.st.prefix = 'a';
%% realign and unwarp
matlabbatch{16}.spm.spatial.realignunwarp.data(1).scans(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 1)', substruct('.','val', '{}',{15}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{16}.spm.spatial.realignunwarp.data(1).pmscan = '';
matlabbatch{16}.spm.spatial.realignunwarp.eoptions.quality = 0.9;
matlabbatch{16}.spm.spatial.realignunwarp.eoptions.sep = 4;
matlabbatch{16}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
matlabbatch{16}.spm.spatial.realignunwarp.eoptions.rtm = 0;
matlabbatch{16}.spm.spatial.realignunwarp.eoptions.einterp = 2;
matlabbatch{16}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
matlabbatch{16}.spm.spatial.realignunwarp.eoptions.weight = '';
matlabbatch{16}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
matlabbatch{16}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
matlabbatch{16}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
matlabbatch{16}.spm.spatial.realignunwarp.uweoptions.jm = 0;
matlabbatch{16}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
matlabbatch{16}.spm.spatial.realignunwarp.uweoptions.sot = [];
matlabbatch{16}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
matlabbatch{16}.spm.spatial.realignunwarp.uweoptions.rem = 1;
matlabbatch{16}.spm.spatial.realignunwarp.uweoptions.noi = 5;
matlabbatch{16}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
matlabbatch{16}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
matlabbatch{16}.spm.spatial.realignunwarp.uwroptions.rinterp = 4;
matlabbatch{16}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
matlabbatch{16}.spm.spatial.realignunwarp.uwroptions.mask = 1;
matlabbatch{16}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';
%% coregistration
matlabbatch{17}.spm.spatial.coreg.estimate.ref(1) = cfg_dep('Image Calculator: ImCalc Computed Image: Brain', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{17}.spm.spatial.coreg.estimate.source(1) = cfg_dep('Realign & Unwarp: Unwarped Mean Image', substruct('.','val', '{}',{16}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','meanuwr'));
matlabbatch{17}.spm.spatial.coreg.estimate.other(1) = cfg_dep('Realign & Unwarp: Unwarped Images (Sess 1)', substruct('.','val', '{}',{16}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','uwrfiles'));
matlabbatch{17}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{17}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{17}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{17}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
%% FIRST LEVEL ITA MONO
matlabbatch{18}.spm.stats.fmri_spec.dir = {ITA_fold};
matlabbatch{18}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{18}.spm.stats.fmri_spec.timing.RT = 2;
matlabbatch{18}.spm.stats.fmri_spec.timing.fmri_t = 35;
matlabbatch{18}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
matlabbatch{18}.spm.stats.fmri_spec.sess(1).scans(1) = cfg_dep('Realign & Unwarp: Unwarped Images (Sess 1)', substruct('.','val', '{}',{16}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','uwrfiles'));

%%
matlabbatch{18}.spm.stats.fmri_spec.sess(1).cond(1).name = 'ITA_mono';
ITA_mono_onsets = dir(fullfile(ITA_fold, '*_ITA_stay.txt'));
ITA_mono_onsets = load(fullfile(ITA_fold,ITA_mono_onsets.name));
matlabbatch{18}.spm.stats.fmri_spec.sess(1).cond(1).onset = ITA_mono_onsets;
matlabbatch{18}.spm.stats.fmri_spec.sess(1).cond(1).duration = 0;
matlabbatch{18}.spm.stats.fmri_spec.sess(1).cond(1).tmod = 0;
matlabbatch{18}.spm.stats.fmri_spec.sess(1).cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{18}.spm.stats.fmri_spec.sess(1).cond(1).orth = 1;

%%
matlabbatch{18}.spm.stats.fmri_spec.sess(1).multi = {''};
matlabbatch{18}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});
matlabbatch{18}.spm.stats.fmri_spec.sess(1).multi_reg = cfg_dep('Realign & Unwarp: Realignment Param File (Sess 1)', substruct('.','val', '{}',{16}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','rpfile'));
matlabbatch{18}.spm.stats.fmri_spec.sess(1).hpf = 128;
matlabbatch{18}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{18}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{18}.spm.stats.fmri_spec.volt = 1;
matlabbatch{18}.spm.stats.fmri_spec.global = 'None';
matlabbatch{18}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{18}.spm.stats.fmri_spec.mask = {''};
matlabbatch{18}.spm.stats.fmri_spec.cvi = 'AR(1)';
%% %% FIRST LEVEL EST AND CONTRASTS
matlabbatch{19}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{18}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{19}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{19}.spm.stats.fmri_est.method.Classical = 1;
%%
matlabbatch{20}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{19}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{20}.spm.stats.con.consess{1}.tcon.name = 'Mono_ITA';
matlabbatch{20}.spm.stats.con.consess{1}.tcon.weights = [1];
matlabbatch{20}.spm.stats.con.consess{1}.tcon.sessrep = 'repl';
matlabbatch{20}.spm.stats.con.delete = 0;
%% MAPPING TO SURFACE
matlabbatch{21}.spm.tools.cat.stools.vol2surf.data_vol(1) = cfg_dep('Contrast Manager: All Con Images', substruct('.','val', '{}',{20}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','con'));
matlabbatch{21}.spm.tools.cat.stools.vol2surf.data_mesh_lh(1) = cfg_dep('CAT12: Segmentation: Left Central Surface', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','lhcentral', '()',{':'}));
matlabbatch{21}.spm.tools.cat.stools.vol2surf.sample = {'maxabs'};
matlabbatch{21}.spm.tools.cat.stools.vol2surf.interp = {'linear'};
matlabbatch{21}.spm.tools.cat.stools.vol2surf.datafieldname = 'Mono_ITA';
matlabbatch{21}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.class = 'GM';
matlabbatch{21}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.startpoint = -0.5;
matlabbatch{21}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.steps = 7;
matlabbatch{21}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.endpoint = 0.5;
matlabbatch{22}.spm.tools.cat.stools.surfresamp.sample{1}.data_surf(1) = cfg_dep('Map Volume (Native Space) to Individual Surface: Left mapped values', substruct('.','val', '{}',{21}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','lh'));
matlabbatch{22}.spm.tools.cat.stools.surfresamp.merge_hemi = 1;
matlabbatch{22}.spm.tools.cat.stools.surfresamp.mesh32k = 1;
matlabbatch{22}.spm.tools.cat.stools.surfresamp.fwhm_surf = 3;
matlabbatch{22}.spm.tools.cat.stools.surfresamp.lazy = 0;
matlabbatch{22}.spm.tools.cat.stools.surfresamp.nproc = 4;
matlabbatch{23}.spm.tools.cat.stools.surfresamp.sample{1}.data_surf(1) = cfg_dep('Map Volume (Native Space) to Individual Surface: Left mapped values', substruct('.','val', '{}',{21}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','lh'));
matlabbatch{23}.spm.tools.cat.stools.surfresamp.merge_hemi = 0;
matlabbatch{23}.spm.tools.cat.stools.surfresamp.mesh32k = 1;
matlabbatch{23}.spm.tools.cat.stools.surfresamp.fwhm_surf = 3;
matlabbatch{23}.spm.tools.cat.stools.surfresamp.lazy = 0;
matlabbatch{23}.spm.tools.cat.stools.surfresamp.nproc = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%ENG_mono%%
  ENG_mono = dir(fullfile(ENG_fold, '*fMRI_Switch_mono_UK_*.nii')); 
  ENGselect = strcat(ENG_mono.name, ',1');
%% image expansion
matlabbatch{24}.spm.util.exp_frames.files = {(fullfile(ENG_fold, ENGselect))};
matlabbatch{24}.spm.util.exp_frames.frames = Inf;
%% slice time correction
matlabbatch{25}.spm.temporal.st.scans{1}(1) = cfg_dep('Expand image frames: Expanded filename list.', substruct('.','val', '{}',{24}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{25}.spm.temporal.st.nslices = 35;
matlabbatch{25}.spm.temporal.st.tr = 2;
matlabbatch{25}.spm.temporal.st.ta = 1.94285714285714;
matlabbatch{25}.spm.temporal.st.so = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34];
matlabbatch{25}.spm.temporal.st.refslice = 1;
matlabbatch{25}.spm.temporal.st.prefix = 'a';
%% realign and unwarp
matlabbatch{26}.spm.spatial.realignunwarp.data(1).scans(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 1)', substruct('.','val', '{}',{25}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{26}.spm.spatial.realignunwarp.data(1).pmscan = '';
matlabbatch{26}.spm.spatial.realignunwarp.eoptions.quality = 0.9;
matlabbatch{26}.spm.spatial.realignunwarp.eoptions.sep = 4;
matlabbatch{26}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
matlabbatch{26}.spm.spatial.realignunwarp.eoptions.rtm = 0;
matlabbatch{26}.spm.spatial.realignunwarp.eoptions.einterp = 2;
matlabbatch{26}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
matlabbatch{26}.spm.spatial.realignunwarp.eoptions.weight = '';
matlabbatch{26}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
matlabbatch{26}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
matlabbatch{26}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
matlabbatch{26}.spm.spatial.realignunwarp.uweoptions.jm = 0;
matlabbatch{26}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
matlabbatch{26}.spm.spatial.realignunwarp.uweoptions.sot = [];
matlabbatch{26}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
matlabbatch{26}.spm.spatial.realignunwarp.uweoptions.rem = 1;
matlabbatch{26}.spm.spatial.realignunwarp.uweoptions.noi = 5;
matlabbatch{26}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
matlabbatch{26}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
matlabbatch{26}.spm.spatial.realignunwarp.uwroptions.rinterp = 4;
matlabbatch{26}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
matlabbatch{26}.spm.spatial.realignunwarp.uwroptions.mask = 1;
matlabbatch{26}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';
%% coregistration
matlabbatch{27}.spm.spatial.coreg.estimate.ref(1) = cfg_dep('Image Calculator: ImCalc Computed Image: Brain', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{27}.spm.spatial.coreg.estimate.source(1) = cfg_dep('Realign & Unwarp: Unwarped Mean Image', substruct('.','val', '{}',{26}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','meanuwr'));
matlabbatch{27}.spm.spatial.coreg.estimate.other(1) = cfg_dep('Realign & Unwarp: Unwarped Images (Sess 1)', substruct('.','val', '{}',{26}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','uwrfiles'));
matlabbatch{27}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{27}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{27}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{27}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
%% FIRST LEVEL ENG
matlabbatch{28}.spm.stats.fmri_spec.dir = {ENG_fold};
matlabbatch{28}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{28}.spm.stats.fmri_spec.timing.RT = 2;
matlabbatch{28}.spm.stats.fmri_spec.timing.fmri_t = 35;
matlabbatch{28}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
matlabbatch{28}.spm.stats.fmri_spec.sess(1).scans(1) = cfg_dep('Realign & Unwarp: Unwarped Images (Sess 1)', substruct('.','val', '{}',{26}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','uwrfiles'));

%%
matlabbatch{28}.spm.stats.fmri_spec.sess(1).cond(1).name = 'ENG_mono';
ENG_mono_onsets = dir(fullfile(ENG_fold, '*_ENG_stay.txt'));
ENG_mono_onsets = load(fullfile(ENG_fold,ENG_mono_onsets.name));
matlabbatch{28}.spm.stats.fmri_spec.sess(1).cond(1).onset = ENG_mono_onsets;
matlabbatch{28}.spm.stats.fmri_spec.sess(1).cond(1).duration = 0;
matlabbatch{28}.spm.stats.fmri_spec.sess(1).cond(1).tmod = 0;
matlabbatch{28}.spm.stats.fmri_spec.sess(1).cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{28}.spm.stats.fmri_spec.sess(1).cond(1).orth = 1;

%%
matlabbatch{28}.spm.stats.fmri_spec.sess(1).multi = {''};
matlabbatch{28}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});
matlabbatch{28}.spm.stats.fmri_spec.sess(1).multi_reg = cfg_dep('Realign & Unwarp: Realignment Param File (Sess 1)', substruct('.','val', '{}',{26}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','rpfile'));
matlabbatch{28}.spm.stats.fmri_spec.sess(1).hpf = 128;
matlabbatch{28}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{28}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{28}.spm.stats.fmri_spec.volt = 1;
matlabbatch{28}.spm.stats.fmri_spec.global = 'None';
matlabbatch{28}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{28}.spm.stats.fmri_spec.mask = {''};
matlabbatch{28}.spm.stats.fmri_spec.cvi = 'AR(1)';
%% %% FIRST LEVEL EST AND CONTRASTS
matlabbatch{29}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{28}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{29}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{29}.spm.stats.fmri_est.method.Classical = 1;
%%
matlabbatch{30}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{29}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{30}.spm.stats.con.consess{1}.tcon.name = 'Mono_ENG';
matlabbatch{30}.spm.stats.con.consess{1}.tcon.weights = [1];
matlabbatch{30}.spm.stats.con.consess{1}.tcon.sessrep = 'repl';
matlabbatch{30}.spm.stats.con.delete = 0;
%% MAPPING TO SURFACE
matlabbatch{31}.spm.tools.cat.stools.vol2surf.data_vol(1) = cfg_dep('Contrast Manager: All Con Images', substruct('.','val', '{}',{30}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','con'));
matlabbatch{31}.spm.tools.cat.stools.vol2surf.data_mesh_lh(1) = cfg_dep('CAT12: Segmentation: Left Central Surface', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','lhcentral', '()',{':'}));
matlabbatch{31}.spm.tools.cat.stools.vol2surf.sample = {'maxabs'};
matlabbatch{31}.spm.tools.cat.stools.vol2surf.interp = {'linear'};
matlabbatch{31}.spm.tools.cat.stools.vol2surf.datafieldname = 'Mono_ENG';
matlabbatch{31}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.class = 'GM';
matlabbatch{31}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.startpoint = -0.5;
matlabbatch{31}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.steps = 7;
matlabbatch{31}.spm.tools.cat.stools.vol2surf.mapping.rel_mapping.endpoint = 0.5;
matlabbatch{32}.spm.tools.cat.stools.surfresamp.sample{1}.data_surf(1) = cfg_dep('Map Volume (Native Space) to Individual Surface: Left mapped values', substruct('.','val', '{}',{31}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','lh'));
matlabbatch{32}.spm.tools.cat.stools.surfresamp.merge_hemi = 1;
matlabbatch{32}.spm.tools.cat.stools.surfresamp.mesh32k = 1;
matlabbatch{32}.spm.tools.cat.stools.surfresamp.fwhm_surf = 3;
matlabbatch{32}.spm.tools.cat.stools.surfresamp.lazy = 0;
matlabbatch{32}.spm.tools.cat.stools.surfresamp.nproc = 4;
matlabbatch{33}.spm.tools.cat.stools.surfresamp.sample{1}.data_surf(1) = cfg_dep('Map Volume (Native Space) to Individual Surface: Left mapped values', substruct('.','val', '{}',{31}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','lh'));
matlabbatch{33}.spm.tools.cat.stools.surfresamp.merge_hemi = 0;
matlabbatch{33}.spm.tools.cat.stools.surfresamp.mesh32k = 1;
matlabbatch{33}.spm.tools.cat.stools.surfresamp.fwhm_surf = 3;
matlabbatch{33}.spm.tools.cat.stools.surfresamp.lazy = 0;
matlabbatch{33}.spm.tools.cat.stools.surfresamp.nproc = 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% VOLUME
%SWITCH
matlabbatch{34}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('CAT12: Segmentation (current release): Deformation Field', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','fordef', '()',{':'}));
matlabbatch{34}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Realign & Unwarp: Unwarped Images (Sess 1)', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','uwrfiles'));
matlabbatch{34}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                           78 76 85];
matlabbatch{34}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
matlabbatch{34}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{34}.spm.spatial.normalise.write.woptions.prefix = 'w';
matlabbatch{35}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('CAT12: Segmentation (current release): Deformation Field', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','fordef', '()',{':'}));
matlabbatch{35}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Realign & Unwarp: Unwarped Images (Sess 2)', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{2}, '.','uwrfiles'));
matlabbatch{35}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                           78 76 85];
matlabbatch{35}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
matlabbatch{35}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{35}.spm.spatial.normalise.write.woptions.prefix = 'w';
matlabbatch{36}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{34}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{36}.spm.spatial.smooth.fwhm = [8 8 8];
matlabbatch{36}.spm.spatial.smooth.dtype = 0;
matlabbatch{36}.spm.spatial.smooth.im = 0;
matlabbatch{36}.spm.spatial.smooth.prefix = 's';
matlabbatch{37}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{35}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{37}.spm.spatial.smooth.fwhm = [8 8 8];
matlabbatch{37}.spm.spatial.smooth.dtype = 0;
matlabbatch{37}.spm.spatial.smooth.im = 0;
matlabbatch{37}.spm.spatial.smooth.prefix = 's';

%% FIRST LEVEL SWITCH VOLUME

matlabbatch{38}.spm.stats.fmri_spec.dir = {Volume_switch};
matlabbatch{38}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{38}.spm.stats.fmri_spec.timing.RT = 2;
matlabbatch{38}.spm.stats.fmri_spec.timing.fmri_t = 35;
matlabbatch{38}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
matlabbatch{38}.spm.stats.fmri_spec.sess(1).scans(1) = cfg_dep('Smooth: Smoothed Images', substruct('.','val', '{}',{36}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{38}.spm.stats.fmri_spec.sess(1).cond(1).name = 'ITA_switch';
%%
matlabbatch{38}.spm.stats.fmri_spec.sess(1).cond(1).onset = ITA_switch_A;
%%
matlabbatch{38}.spm.stats.fmri_spec.sess(1).cond(1).duration = 0;
matlabbatch{38}.spm.stats.fmri_spec.sess(1).cond(1).tmod = 0;
matlabbatch{38}.spm.stats.fmri_spec.sess(1).cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{38}.spm.stats.fmri_spec.sess(1).cond(1).orth = 1;
matlabbatch{38}.spm.stats.fmri_spec.sess(1).cond(2).name = 'ENG_switch';
%%
matlabbatch{38}.spm.stats.fmri_spec.sess(1).cond(2).onset = ENG_switch_A;
%%
matlabbatch{38}.spm.stats.fmri_spec.sess(1).cond(2).duration = 0;
matlabbatch{38}.spm.stats.fmri_spec.sess(1).cond(2).tmod = 0;
matlabbatch{38}.spm.stats.fmri_spec.sess(1).cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{38}.spm.stats.fmri_spec.sess(1).cond(2).orth = 1;
matlabbatch{38}.spm.stats.fmri_spec.sess(1).cond(3).name = 'ITA_stay';
%%
matlabbatch{38}.spm.stats.fmri_spec.sess(1).cond(3).onset = ITA_stay_A;
%%
matlabbatch{38}.spm.stats.fmri_spec.sess(1).cond(3).duration = 0;
matlabbatch{38}.spm.stats.fmri_spec.sess(1).cond(3).tmod = 0;
matlabbatch{38}.spm.stats.fmri_spec.sess(1).cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{38}.spm.stats.fmri_spec.sess(1).cond(3).orth = 1;
matlabbatch{38}.spm.stats.fmri_spec.sess(1).cond(4).name = 'ENG_stay';
%%
matlabbatch{38}.spm.stats.fmri_spec.sess(1).cond(4).onset = ENG_stay_A;
%%
matlabbatch{38}.spm.stats.fmri_spec.sess(1).cond(4).duration = 0;
matlabbatch{38}.spm.stats.fmri_spec.sess(1).cond(4).tmod = 0;
matlabbatch{38}.spm.stats.fmri_spec.sess(1).cond(4).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{38}.spm.stats.fmri_spec.sess(1).cond(4).orth = 1;
matlabbatch{38}.spm.stats.fmri_spec.sess(1).multi = {''};
matlabbatch{38}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});
matlabbatch{38}.spm.stats.fmri_spec.sess(1).multi_reg(1) = cfg_dep('Realign & Unwarp: Realignment Param File (Sess 1)', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','rpfile'));
matlabbatch{38}.spm.stats.fmri_spec.sess(1).hpf = 128;
matlabbatch{38}.spm.stats.fmri_spec.sess(2).scans(1) = cfg_dep('Smooth: Smoothed Images', substruct('.','val', '{}',{37}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{38}.spm.stats.fmri_spec.sess(2).cond(1).name = 'ITA_switch';
%%
matlabbatch{38}.spm.stats.fmri_spec.sess(2).cond(1).onset = ITA_switch_B;
%%
matlabbatch{38}.spm.stats.fmri_spec.sess(2).cond(1).duration = 0;
matlabbatch{38}.spm.stats.fmri_spec.sess(2).cond(1).tmod = 0;
matlabbatch{38}.spm.stats.fmri_spec.sess(2).cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{38}.spm.stats.fmri_spec.sess(2).cond(1).orth = 1;
matlabbatch{38}.spm.stats.fmri_spec.sess(2).cond(2).name = 'ENG_switch';
%%
matlabbatch{38}.spm.stats.fmri_spec.sess(2).cond(2).onset = ENG_switch_B;
%%
matlabbatch{38}.spm.stats.fmri_spec.sess(2).cond(2).duration = 0;
matlabbatch{38}.spm.stats.fmri_spec.sess(2).cond(2).tmod = 0;
matlabbatch{38}.spm.stats.fmri_spec.sess(2).cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{38}.spm.stats.fmri_spec.sess(2).cond(2).orth = 1;
matlabbatch{38}.spm.stats.fmri_spec.sess(2).cond(3).name = 'ITA_stay';
%%
matlabbatch{38}.spm.stats.fmri_spec.sess(2).cond(3).onset = ITA_stay_B;
%%
matlabbatch{38}.spm.stats.fmri_spec.sess(2).cond(3).duration = 0;
matlabbatch{38}.spm.stats.fmri_spec.sess(2).cond(3).tmod = 0;
matlabbatch{38}.spm.stats.fmri_spec.sess(2).cond(3).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{38}.spm.stats.fmri_spec.sess(2).cond(3).orth = 1;
matlabbatch{38}.spm.stats.fmri_spec.sess(2).cond(4).name = 'ENG_stay';
%%
matlabbatch{38}.spm.stats.fmri_spec.sess(2).cond(4).onset = ENG_stay_B;
%%
matlabbatch{38}.spm.stats.fmri_spec.sess(2).cond(4).duration = 0;
matlabbatch{38}.spm.stats.fmri_spec.sess(2).cond(4).tmod = 0;
matlabbatch{38}.spm.stats.fmri_spec.sess(2).cond(4).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{38}.spm.stats.fmri_spec.sess(2).cond(4).orth = 1;
matlabbatch{38}.spm.stats.fmri_spec.sess(2).multi = {''};
matlabbatch{38}.spm.stats.fmri_spec.sess(2).regress = struct('name', {}, 'val', {});
matlabbatch{38}.spm.stats.fmri_spec.sess(2).multi_reg(1) = cfg_dep('Realign & Unwarp: Realignment Param File (Sess 2)', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{2}, '.','rpfile'));
matlabbatch{38}.spm.stats.fmri_spec.sess(2).hpf = 128;
matlabbatch{38}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{38}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{38}.spm.stats.fmri_spec.volt = 1;
matlabbatch{38}.spm.stats.fmri_spec.global = 'None';
matlabbatch{38}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{38}.spm.stats.fmri_spec.mask = {''};
matlabbatch{38}.spm.stats.fmri_spec.cvi = 'AR(1)';
matlabbatch{39}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{38}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{39}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{39}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{40}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{39}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{40}.spm.stats.con.consess{1}.tcon.name = 'Main_ITA_switch';
matlabbatch{40}.spm.stats.con.consess{1}.tcon.weights = [1 0 0 0];
matlabbatch{40}.spm.stats.con.consess{1}.tcon.sessrep = 'repl';
matlabbatch{40}.spm.stats.con.consess{2}.tcon.name = 'Main_ENG_switch';
matlabbatch{40}.spm.stats.con.consess{2}.tcon.weights = [0 1 0 0];
matlabbatch{40}.spm.stats.con.consess{2}.tcon.sessrep = 'repl';
matlabbatch{40}.spm.stats.con.consess{3}.tcon.name = 'Main_ITA_stay';
matlabbatch{40}.spm.stats.con.consess{3}.tcon.weights = [0 0 1 0];
matlabbatch{40}.spm.stats.con.consess{3}.tcon.sessrep = 'repl';
matlabbatch{40}.spm.stats.con.consess{4}.tcon.name = 'Main_ENG_stay';
matlabbatch{40}.spm.stats.con.consess{4}.tcon.weights = [0 0 0 1];
matlabbatch{40}.spm.stats.con.consess{4}.tcon.sessrep = 'repl';
matlabbatch{40}.spm.stats.con.consess{5}.tcon.name = 'Main_switch';
matlabbatch{40}.spm.stats.con.consess{5}.tcon.weights = [1 1 0 0];
matlabbatch{40}.spm.stats.con.consess{5}.tcon.sessrep = 'repl';
matlabbatch{40}.spm.stats.con.consess{6}.tcon.name = 'Main_stay';
matlabbatch{40}.spm.stats.con.consess{6}.tcon.weights = [0 0 1 1];
matlabbatch{40}.spm.stats.con.consess{6}.tcon.sessrep = 'repl';
matlabbatch{40}.spm.stats.con.consess{7}.tcon.name = 'Main_ITA';
matlabbatch{40}.spm.stats.con.consess{7}.tcon.weights = [1 0 1 0];
matlabbatch{40}.spm.stats.con.consess{7}.tcon.sessrep = 'repl';
matlabbatch{40}.spm.stats.con.consess{8}.tcon.name = 'Main_ENG';
matlabbatch{40}.spm.stats.con.consess{8}.tcon.weights = [0 1 0 1];
matlabbatch{40}.spm.stats.con.consess{8}.tcon.sessrep = 'repl';
matlabbatch{40}.spm.stats.con.consess{9}.tcon.name = 'switch_gt_stay';
matlabbatch{40}.spm.stats.con.consess{9}.tcon.weights = [1 1 -1 -1];
matlabbatch{40}.spm.stats.con.consess{9}.tcon.sessrep = 'repl';
matlabbatch{40}.spm.stats.con.consess{10}.tcon.name = 'stay_gt_switch';
matlabbatch{40}.spm.stats.con.consess{10}.tcon.weights = [-1 -1 1 1];
matlabbatch{40}.spm.stats.con.consess{10}.tcon.sessrep = 'repl';
matlabbatch{40}.spm.stats.con.consess{11}.tcon.name = 'ITA_gt_ENG';
matlabbatch{40}.spm.stats.con.consess{11}.tcon.weights = [1 -1 1 -1];
matlabbatch{40}.spm.stats.con.consess{11}.tcon.sessrep = 'repl';
matlabbatch{40}.spm.stats.con.consess{12}.tcon.name = 'ENG_gt_ITA';
matlabbatch{40}.spm.stats.con.consess{12}.tcon.weights = [-1 1 -1 1];
matlabbatch{40}.spm.stats.con.consess{12}.tcon.sessrep = 'repl';
matlabbatch{40}.spm.stats.con.consess{13}.tcon.name = 'ITA_switch_gt_ENG_switch';
matlabbatch{40}.spm.stats.con.consess{13}.tcon.weights = [1 -1 0 0];
matlabbatch{40}.spm.stats.con.consess{13}.tcon.sessrep = 'repl';
matlabbatch{40}.spm.stats.con.consess{14}.tcon.name = 'ENG_switch_gt_ITA_switch';
matlabbatch{40}.spm.stats.con.consess{14}.tcon.weights = [-1 1 0 0];
matlabbatch{40}.spm.stats.con.consess{14}.tcon.sessrep = 'repl';
matlabbatch{40}.spm.stats.con.consess{15}.tcon.name = 'Inter_1';
matlabbatch{40}.spm.stats.con.consess{15}.tcon.weights = [-1 1 1 -1];
matlabbatch{40}.spm.stats.con.consess{15}.tcon.sessrep = 'repl';
matlabbatch{40}.spm.stats.con.consess{16}.tcon.name = 'Inter_2';
matlabbatch{40}.spm.stats.con.consess{16}.tcon.weights = [1 -1 -1 1];
matlabbatch{40}.spm.stats.con.consess{16}.tcon.sessrep = 'repl';
matlabbatch{40}.spm.stats.con.delete = 0;
%%%%%Normalize Brain
matlabbatch{41}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('CAT12: Segmentation (current release): Deformation Field', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','fordef', '()',{':'}));
matlabbatch{41}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Image Calculator: ImCalc Computed Image: Brain', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{41}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                           78 76 85];
matlabbatch{41}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
matlabbatch{41}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{41}.spm.spatial.normalise.write.woptions.prefix = 'w';
% %%%%%%% ITA MONO VOLUME
% 
matlabbatch{42}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('CAT12: Segmentation (current release): Deformation Field', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','fordef', '()',{':'}));
matlabbatch{42}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Realign & Unwarp: Unwarped Images (Sess 1)', substruct('.','val', '{}',{16}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','uwrfiles'));
matlabbatch{42}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                           78 76 85];
matlabbatch{42}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
matlabbatch{42}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{42}.spm.spatial.normalise.write.woptions.prefix = 'w';
matlabbatch{43}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{42}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{43}.spm.spatial.smooth.fwhm = [8 8 8];
matlabbatch{43}.spm.spatial.smooth.dtype = 0;
matlabbatch{43}.spm.spatial.smooth.im = 0;
matlabbatch{43}.spm.spatial.smooth.prefix = 's';
matlabbatch{44}.spm.stats.fmri_spec.dir = {Volume_Ita};
matlabbatch{44}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{44}.spm.stats.fmri_spec.timing.RT = 2;
matlabbatch{44}.spm.stats.fmri_spec.timing.fmri_t = 35;
matlabbatch{44}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
matlabbatch{44}.spm.stats.fmri_spec.sess.scans(1) = cfg_dep('Smooth: Smoothed Images', substruct('.','val', '{}',{43}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{44}.spm.stats.fmri_spec.sess.cond.name = 'ITA_mono';
%%
matlabbatch{44}.spm.stats.fmri_spec.sess.cond.onset = ITA_mono_onsets;
%%
matlabbatch{44}.spm.stats.fmri_spec.sess.cond.duration = 0;
matlabbatch{44}.spm.stats.fmri_spec.sess.cond.tmod = 0;
matlabbatch{44}.spm.stats.fmri_spec.sess.cond.pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{44}.spm.stats.fmri_spec.sess.cond.orth = 1;
matlabbatch{44}.spm.stats.fmri_spec.sess.multi = {''};
matlabbatch{44}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
matlabbatch{44}.spm.stats.fmri_spec.sess.multi_reg(1) = cfg_dep('Realign & Unwarp: Realignment Param File (Sess 1)', substruct('.','val', '{}',{16}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','rpfile'));
matlabbatch{44}.spm.stats.fmri_spec.sess.hpf = 128;
matlabbatch{44}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{44}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{44}.spm.stats.fmri_spec.volt = 1;
matlabbatch{44}.spm.stats.fmri_spec.global = 'None';
matlabbatch{44}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{44}.spm.stats.fmri_spec.mask = {''};
matlabbatch{44}.spm.stats.fmri_spec.cvi = 'AR(1)';
matlabbatch{45}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{44}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{45}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{45}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{46}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{45}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{46}.spm.stats.con.consess{1}.tcon.name = 'Mono_ITA';
matlabbatch{46}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{46}.spm.stats.con.consess{1}.tcon.sessrep = 'repl';
matlabbatch{46}.spm.stats.con.delete = 0;

%%% ENG MONO VOLUME

matlabbatch{47}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('CAT12: Segmentation (current release): Deformation Field', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','fordef', '()',{':'}));
matlabbatch{47}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Realign & Unwarp: Unwarped Images (Sess 1)', substruct('.','val', '{}',{26}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','uwrfiles'));
matlabbatch{47}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                           78 76 85];
matlabbatch{47}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
matlabbatch{47}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{47}.spm.spatial.normalise.write.woptions.prefix = 'w';
matlabbatch{48}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{47}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{48}.spm.spatial.smooth.fwhm = [8 8 8];
matlabbatch{48}.spm.spatial.smooth.dtype = 0;
matlabbatch{48}.spm.spatial.smooth.im = 0;
matlabbatch{48}.spm.spatial.smooth.prefix = 's';
matlabbatch{49}.spm.stats.fmri_spec.dir = {Volume_Eng};
matlabbatch{49}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{49}.spm.stats.fmri_spec.timing.RT = 2;
matlabbatch{49}.spm.stats.fmri_spec.timing.fmri_t = 35;
matlabbatch{49}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
matlabbatch{49}.spm.stats.fmri_spec.sess.scans(1) = cfg_dep('Smooth: Smoothed Images', substruct('.','val', '{}',{48}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{49}.spm.stats.fmri_spec.sess.cond.name = 'ENG_mono';
%%
matlabbatch{49}.spm.stats.fmri_spec.sess.cond.onset = ENG_mono_onsets;
%%
matlabbatch{49}.spm.stats.fmri_spec.sess.cond.duration = 0;
matlabbatch{49}.spm.stats.fmri_spec.sess.cond.tmod = 0;
matlabbatch{49}.spm.stats.fmri_spec.sess.cond.pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{49}.spm.stats.fmri_spec.sess.cond.orth = 1;
matlabbatch{49}.spm.stats.fmri_spec.sess.multi = {''};
matlabbatch{49}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
matlabbatch{49}.spm.stats.fmri_spec.sess.multi_reg(1) = cfg_dep('Realign & Unwarp: Realignment Param File (Sess 1)', substruct('.','val', '{}',{26}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','rpfile'));
matlabbatch{49}.spm.stats.fmri_spec.sess.hpf = 128;
matlabbatch{49}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{49}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{49}.spm.stats.fmri_spec.volt = 1;
matlabbatch{49}.spm.stats.fmri_spec.global = 'None';
matlabbatch{49}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{49}.spm.stats.fmri_spec.mask = {''};
matlabbatch{49}.spm.stats.fmri_spec.cvi = 'AR(1)';
matlabbatch{50}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{49}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{50}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{50}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{51}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{50}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{51}.spm.stats.con.consess{1}.tcon.name = 'Mono_ENG';
matlabbatch{51}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{51}.spm.stats.con.consess{1}.tcon.sessrep = 'repl';
matlabbatch{51}.spm.stats.con.delete = 0;

%%

save (['matlabbatch' subname '.mat'], 'matlabbatch');
    
%%%it's time to do it!
%spm_jobman('run',matlabbatch);

fprintf('Finished folder %d of %d: %s\n', ...
  k, numberOfFolders, thisDir);
%end
