addpath(genpath(pwd))
inputParams
tt=tic;
simInput = initSimFunctionPolarPPM(EbN0dB,codeLens,codeRates,Qmat,Mmat,conflvl);
mainSimFunctionPolarPPM(simInput,transmissions,use_canonical_ppm)
% mainSimFunctionPolarPPM2(simInput,transmissions,use_canonical_ppm)
t2 = toc(tt);
fprintf('Sim took %f hours\n',t2/60/60)
% resultmatrix = showResults(EbN0dB,codeLens,codeRates,Qmat,Mmat,conflvl);
folder_name = 'results';
resultmatrix = showResultsV2(EbN0dB,codeLens,codeRates,Qmat,Mmat,folder_name,use_canonical_ppm);