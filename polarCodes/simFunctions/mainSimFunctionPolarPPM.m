function mainSimFunctionPolarPPM(simInput,transmissions,use_canonical_ppm)
totSims = length(simInput);
if exist('counter.csv','file')==2
   general_counter_init =  csvread('counter.csv');
else
    general_counter_init = [];
end
allSims = 1:totSims;
allSims(general_counter_init) = [];
totSims = length(allSims);
parfor ss = 1:totSims
    general_counter = allSims(ss);
    t1=tic;
    codeLen = simInput(general_counter).codeLen;
    codeRate = simInput(general_counter).codeRate;
    M = simInput(general_counter).M;
    Q = simInput(general_counter).Q;
    EbN0dB = simInput(general_counter).EbN0dB;
    alphaconf = simInput(general_counter).alphaconf;
    coefM = f1coeffs(M);
    %     bers = zeros(length(EbN0dB),3);
    ber = zeros(1,3);
    if use_canonical_ppm
        codeLen = codeLen/log2(Q);
        [GmatT,databitspos,~,frozenbitchecks,dataLen] = PolarConstr(codeLen,codeRate);
        dataLen = dataLen*log2(Q);
    else
        [GmatT,databitspos,~,frozenbitchecks,dataLen] = PolarConstr(codeLen,codeRate);
    end
    fid =fopen('simLog.txt','a+');
    fid_counter = fopen('counter.csv','a+');
    EbN0 = log2(Q) * codeRate * 10^(EbN0dB/10);
    errors = 0;
%     h1 = 1.0; h2=0.0;
%     transmissions = 1e4;
%     while (h1-h2)*transmissions*dataLen >= 0.1*errors
    
    for tt=1:transmissions
        %% Trasnmitter
        databits = randi([0, 1], dataLen, 1);
        if use_canonical_ppm
             encodedBits = canonicalPPMPolarEncoder(databits,codeLen, GmatT, databitspos,Q);
        else
            encodedBits = PolarEncoder(databits,codeLen, GmatT, databitspos);
        end

        [transmitted,zeroPos,onePos] = ppmModulation(Q,encodedBits);
        %% Channel
        noise = randn(size(transmitted,1),2*M)./sqrt(2*EbN0);
        x = [transmitted zeros(size(transmitted,1),2*M-1)];
        
        hx = x;
        hx2 = zeros(size(hx));
        for ii =1:2*M
            hx2(:,ii) = (abs(hx(:,ii)+noise(:,ii) )).^2;
        end
        y_channel = sum(hx2,2);
        %         y_channel =transmitted;
        %% Receiver
        if use_canonical_ppm
            PPMsymbols = codeLen;
        else
            PPMsymbols = codeLen/log2(Q);
        end
        LLR = ppmLLRs(y_channel,PPMsymbols, Q, EbN0, coefM, zeroPos, onePos);
        
        if use_canonical_ppm
            decodebitsdemoded = canonicalPPMPolarDecoder(LLR, codeLen, GmatT, databitspos, frozenbitchecks,dataLen,Q);
        else
            decodebitsdemoded = PolarDecoder(LLR, codeLen, GmatT, databitspos, frozenbitchecks)';
        end
        
        %% Error counting
        currentErrors = sum(mod(decodebitsdemoded+databits,2));
%         currentErrors2 = 1-sum(decodebitsdemoded==databits)/dataLen;
        
        errors = errors + currentErrors;
%         transmissions = transmissions+1;
%         if  transmissions*dataLen > errors
%             h1 = betaincinv(1 - alphaconf/2, errors + 1, transmissions*dataLen - errors);
%         end
%         if errors > 0
%             h2 = betaincinv(alphaconf/2, errors, transmissions*dataLen - errors + 1);
%         end
    end
    ber(1,1) = errors;
    ber(1,2) = 1-errors/dataLen/transmissions;
%     ber(1,3) = transmissions;
    
    res=matfile(sprintf('results/res_%d.mat', general_counter),'writable',true);
    res.ber = ber;
    res.Q = Q;
    res.M = M;
    res.codeRate = codeRate;
    res.codeLen = codeLen;
    res.EbN0dB = EbN0dB;
%     g=matfile('general_counter.mat','writable',true);
%     g.general_counter = [g.general_counter;general_counter];
    fprintf(fid_counter,"%d\n",general_counter);
    fclose(fid_counter);
    t2 = toc(t1);
    fprintf("Q:%d, M:%d, codeRate:%f, codeLen:2^%d, EbN0:%f, %d/%d ...errors:%f, BEP:%f, transmissions:%d - took %f mins\n"...
        ,Q,M,codeRate,log2(codeLen),EbN0dB,general_counter,totSims,ber(1,1), ber(1,2),ber(1,3),t2/60);
    fprintf(fid,"Q:%d, M:%d, codeRate:%f, codeLen:2^%d, EbN0:%f, %d/%d ...errors:%f, BEP:%f, transmissions:%d - took %f mins\n"...
        ,Q,M,codeRate,log2(codeLen),EbN0dB,general_counter,totSims, ber(1,1),ber(1,2),ber(1,3),t2/60);
    fclose(fid);
end
end