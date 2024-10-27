function resultmatrix = showResultsV2(EbN0dB,codeLens,codeRates,Qmat,Mmat,folder_name,use_canonical)
resultmatrix =zeros(length(codeLens),length(codeRates),length(Mmat),length(Qmat),length(EbN0dB));
legendstr = {};
d = dir([folder_name '/*.mat']);
fnames = {d.name};


for jj =1:length(fnames)
    fname = [folder_name '/res_' num2str(jj) '.mat'];
    if exist(fname,'file')
        res = load(fname);
        if use_canonical
            ll = find(codeLens==res.codeLen*log2(res.Q));
        else
            ll = find(codeLens==res.codeLen);
        end
        cc = find(codeRates==res.codeRate);
        mm = find(Mmat==res.M);
        qq = find(Qmat==res.Q);
        kk = find(EbN0dB==res.EbN0dB);
        if sum(isempty([cc;mm;qq;kk;ll]))==0
            resultmatrix(ll,cc,mm,qq,kk) = 1-res.ber(2);
        else
            resultmatrix(ll,cc,mm,qq,kk) =nan;
        end
    end
    
end

resulttable = table();
resulttable.EbN0dB = EbN0dB';
for ll =1:length(codeLens)
    codeLen = codeLens(ll);
    for cc = 1:length(codeRates)
        codeRate = codeRates(cc);
        for mm =1:length(Mmat)
            M = Mmat(mm);
            for qq =1:length(Qmat)
                Q=Qmat(qq);
                currname = "Q"+num2str(Q)+"M"+num2str(M)+"C"+num2str(codeRate,2)+"L"+num2str(codeLen);
                resulttable.(currname) = squeeze(resultmatrix(ll,cc,mm,qq,:));
            end
        end
    end
end

writetable(resulttable,'resultscsv/polarppm.csv')


figure;
for ll =1:length(codeLens)
    codeLen = codeLens(ll);
    for cc = 1:length(codeRates)
        codeRate = codeRates(cc);
        for mm =1:length(Mmat)
            M = Mmat(mm);
            for qq =1:length(Qmat)
                Q = Qmat(qq);
                semilogy(EbN0dB,squeeze(resultmatrix(ll,cc,mm,qq,:)),'Linewidth', 2)
                hold on
                legendstr{end+1}=sprintf("codeLen: 2^{%d}, codeRate :%.4f, M :%d, Q: %d",log2(codeLen),codeRate,M,Q);
                
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