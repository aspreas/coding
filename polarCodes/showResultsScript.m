addpath(genpath(pwd))

inputParams;
folder_name = 'results';
use_canonical = true;
resultmatrix2 = showResultsV2(EbN0dB,codeLens,codeRates,Qmat,Mmat,folder_name,use_canonical);


folder_name = 'results_4';
use_canonical = false;
resultmatrix1 = showResultsV2(EbN0dB,codeLens,codeRates,Qmat,Mmat,folder_name,use_canonical);

showResultsV3(EbN0dB,codeLens,codeRates,Qmat,Mmat,resultmatrix1,resultmatrix2)
