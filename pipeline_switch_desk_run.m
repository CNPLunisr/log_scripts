% List of open inputs
nrun = 1; % enter the number of runs here
jobfile = {'C:\Users\CNPL\Desktop\switch_logs\pipeline\pipeline_switch_desk.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
cat12('expert');
spm_jobman('run', jobs, inputs{:});