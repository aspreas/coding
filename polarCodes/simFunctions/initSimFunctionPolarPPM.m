function simInput = initSimFunctionPolarPPM(EbN0dB,codeLens,codeRates,Qmat,Mmat,conflvl)
alphaconf = 1 - conflvl;
general_counter = 0;
resultmatrix =zeros(length(codeLens),length(codeRates),length(Mmat),length(Qmat),length(EbN0dB));
for ll =1:length(codeLens)
    codeLen = codeLens(ll);
    for cc = 1:length(codeRates)
        codeRate = codeRates(cc);
        for mm =1:length(Mmat)
            M = Mmat(mm);
            for qq =1:length(Qmat)
                Q = Qmat(qq);
                for kk =1:length(EbN0dB)
                    general_counter = general_counter+1;
                    simInput(general_counter).codeLen = codeLen;
                    simInput(general_counter).codeRate = codeRate;
                    simInput(general_counter).M = M;
                    simInput(general_counter).Q = Q;
                    simInput(general_counter).EbN0dB = EbN0dB(kk);
                    simInput(general_counter).conflvl = conflvl;
                    simInput(general_counter).alphaconf = alphaconf;
                end
            end
        end
    end
end
end