function showResultsV3(EbN0dB,codeLens,codeRates,Qmat,Mmat,resultmatrix1,resultmatrix2)
figure;
legendstr={};
for ll =1:length(codeLens)
    codeLen = codeLens(ll);
    for cc = 1:length(codeRates)
        codeRate = codeRates(cc);
        for mm =1:length(Mmat)
            M = Mmat(mm);
            for qq =1:length(Qmat)
                Q = Qmat(qq);
                semilogy(EbN0dB,squeeze(resultmatrix1(ll,cc,mm,qq,:)),'Linewidth', 2)
                hold on
                semilogy(EbN0dB,squeeze(resultmatrix2(ll,cc,mm,qq,:)),'Linewidth', 2)

                legendstr{end+1}=sprintf("codeLen: 2^{%d}, codeRate :%.4f, M :%d, Q: %d",log2(codeLen),codeRate,M,Q);
                legendstr{end+1}=sprintf("canonical codeLen: 2^{%d}, codeRate :%.4f, M :%d, Q: %d",log2(codeLen),codeRate,M,Q);

                legend( legendstr,'location', 'best','ItemHitFcn',@cb_legend);
                grid on;
            end
        end
    end
end
end

function cb_legend(~,evt)
if strcmp(evt.Peer.Visible,'on')
    evt.Peer.Visible = 'off';
else 
    evt.Peer.Visible = 'on';
end
end