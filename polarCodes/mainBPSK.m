addpath(genpath(pwd))

codeLen = 2^3;
codeRate = 0.5;

conflvl = 0.99;
alphaconf = 1 - conflvl;
EbN0dB = 0:0.5:7;
bers = zeros(length(EbN0dB),3);
[GmatT,databitspos,frozenbitspos,frozenbitchecks,dataLen] = PolarConstr(codeLen,codeRate);
for kk =1:length(EbN0dB)
    disp(EbN0dB(kk))
    EbN0 =codeRate * 10^(EbN0dB(kk)/10);
    errors = 0;
    h1 = 1.0; h2=0.0;
    transmissions = 0;
    while (h1-h2)*transmissions*dataLen >= 0.1*errors

        databits = randi([0, 1], dataLen, 1);
        %         disp(databits')
        codebits = PolarEncoder(databits,codeLen, GmatT, databitspos,dataLen);
        LLR = (1 - 2 * codebits) + randn(codeLen, 1) / sqrt(2 * EbN0);
        decodebits = PolarDecoder(LLR, codeLen, GmatT, databitspos, frozenbitchecks);
        decodebitsdemoded =zeros(length(decodebits),1);
        decodebitsdemoded(decodebits<=0)=1;
        
        currentErrors = sum(mod(decodebitsdemoded+databits,2));
        currentErrors2 = 1-sum(decodebitsdemoded==databits)/dataLen;

        errors = errors + currentErrors;
             transmissions = transmissions+1;
        if  transmissions*dataLen > errors
            h1 = betaincinv(1 - alphaconf/2, errors + 1, transmissions*dataLen - errors);
        end
        if errors > 0
            h2 = betaincinv(alphaconf/2, errors, transmissions*dataLen - errors + 1);
        end
    end
    bers(kk,1) = errors;
    bers(kk,2) = errors/dataLen/transmissions;
    bers(kk,3) = transmissions;

end

%% 
figure;
semilogy(EbN0dB,bers(:,2))
hold on
grid on
semilogy(EbN0dB,0.5*erfc(sqrt(10.^(EbN0dB/10))))
title('BPSK polar codes')
ylabel('BEP')
xlabel('EbNo')
legend('polar','theory')
