

exp_path = ['C:/Users/Antonis/Desktop/myMatlabFiles/rsResults/exp_' num2str(expirement_count)];
cnq_path = 'C:/Users/Antonis/Desktop/myMatlabFiles/mathematica/cnqFiles/cnq';

N=1e6; 
Qmat = [4;16];
allM = [2];
% Qmat = [16];
% allM = [40];

step = 0.5;
PbitB_min = 00.0;
PbitB_max = 20.0;
% PbitB=[11.05;11.1] ;
PbitB=PbitB_min:step:PbitB_max;

% 
nn_pair = [255;127;63;];
kk_pair = [239 ;117;57;];
% nn_pair = [255];%511
% kk_pair = [239];%338
% nn_pair = [127];
% kk_pair = [117];
% 511,1023,2047
useInterleving = false; 

run_parallel = true;
max_cores =8;
cmode = 1;%
% channel = 2;% exp(-h)
% channel = 3;% gamma-gamma