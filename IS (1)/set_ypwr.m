% Calibrates 'D:\PhD\06 Experiments\yPRESS\PwrCalib.mat'
%            'D:\PhD\06 Experiments\yPhaseMap GE\PwrCalib.mat'
% Input: pw360 [us]
function set_ypwr(pw360)
load 'D:\PhD\06 Experiments\yPRESS\PwrCalib.mat'
pw90  = pw360/4;
f_max = 1/(pw90*4*1e-3);
tpwr  = 60;
save 'D:\PhD\06 Experiments\yPRESS\PwrCalib.mat' pw90 f_max tpwr

load 'D:\PhD\06 Experiments\yPhaseMap GE\PwrCalib.mat'
pw90  = pw360/4;
f_max = 1/(pw90*4*1e-3);
tpwr  = 60;
save 'D:\PhD\06 Experiments\yPhaseMap GE\PwrCalib.mat' pw90 f_max tpwr

return;
