function resultmatrix = showResults(EbN0dB,codeLens,codeRates,Qmat,Mmat,conflvl)
alphaconf = 1 - conflvl;
general_counter = 0;
resultmatrix =zeros(length(codeLens),length(codeRates),length(Mmat),length(Qmat),length(EbN0dB));
legendstr = {};a
figure;
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
                    fname = ['Copy_of_results/res_' num2str(general_counter) '.mat'];
                    if exist(fname)==2
                    res = load(fname);
                    resultmatrix(ll,cc,mm,qq,kk) = 1-res.ber(2);
                    else
                    resultmatrix(ll,cc,mm,qq,kk) =nan;
                    end
                end
                semilogy(EbN0dB,squeeze(resultmatrix(ll,cc,mm,qq,:)),'Linewidth', 2)
                hold on
                legendstr{end+1}=sprintf("codeLen: 2^{%d}, codeRate :%.1f, M :%d, Q: %d",log2(codeLen),codeRate,M,Q);
            end
        end
    end
end
                legend( legendstr,'location', 'best','ItemHitFcn',@cb_legend);
                grid on;

end
function cb_legend(~,evt)
if strcmp(evt.Peer.Visible,'on')
    evt.Peer.Visible = 'off';
else 
    evt.Peer.Visible = 'on';
end
end