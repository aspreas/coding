EbN0dB = 0:0.3:15;
% EbN0dB = ;

codeLens = 2^10;%2^11;
% codeRates = [1/4;1/2;3/4;7/8];
% Qmat = [4;16];
% Mmat = [2;40;200];
Qmat = [4];
Mmat = [2;40];
codeRates = [1/3;0.5;0.75;0.875;0.9375];
conflvl = 0.95;
transmissions = 1e4;
use_canonical_ppm = false;


%6.9.8.1.4.5.7.0.8.6 Ρια27