
for iPair = 1:length(kk_pair)
    
    mm_small = mm_pair(iPair);
    kk_small = kk_pair(iPair);
    nn_small = nn_pair(iPair);
    
    codingrate = kk_small/nn_small;
    fprintf('Start RS Pair: %d,%d\n',nn_small,kk_small)

    for mm = 1:length(allM)
        M = allM(mm);
        for qq = 1:length(Qmat)
            
            Q = Qmat(qq);
            Qrs = Qrsmat(qq);
            
            berRs = zeros(SS,1);
            berRs_inter = zeros(SS,1);
            berRs_packet = zeros(SS,1);
            berRs_packet_inter  = zeros(SS,1);
            
            berVanila = zeros(SS,1);
            berRs_gamma = zeros(SS,1);
            berRs_packet_gamma= zeros(SS,1);
            bitsPerSymbol = zeros(N,SS);
            
            rsEncoder = comm.RSEncoder('BitInput',true,'CodewordLength',nn_small,'MessageLength',kk_small);
            rsDecoder = comm.RSDecoder('BitInput',true,'CodewordLength',nn_small,'MessageLength',kk_small);
            
            if run_parallel
                for i_pwr = 1:SS
                    powerp =Pbit_lin(i_pwr)*(log2(Qrs));%codingrate
                    fprintf('%f dB ', PbitB(i_pwr))
                    [berRs(i_pwr),~,~] = RSPart_V4(powerp,N,Qrs,M,mm_small,kk_small,rsEncoder,rsDecoder,cmode,useInterleving);
%                     berVanila(i_pwr) = VanilaPart(Pbit_lin(i_pwr)*(log2(Q)),1,Q,M,nn_small,kk_small,cmode);
                    fprintf('Done for %f dB\n', PbitB(i_pwr))
                end
            else
                for i_pwr = 1:SS
%                     codingrate
                    [berRs(i_pwr),berRs_packet(i_pwr),bitsPerSymbol(:,i_pwr)] = RSPart_V4(Pbit_lin(i_pwr)*(log2(Qrs)),N,Qrs,M,mm_small,kk_small,rsEncoder,rsDecoder,cmode,useInterleving);
                    berVanila(i_pwr) = VanilaPart(Pbit_lin(i_pwr)*(log2(Q)),N,Q,M,nn_small,kk_small,cmode);
                    fprintf('Done for %f dB\n', PbitB(i_pwr))
                end
            end
            simStruct.berVanila = berVanila;
            simStruct.berRs = berRs;
            
            simStruct.berRs_inter = berRs_inter;
            simStruct.berRs_packet_inter = berRs_packet_inter;

            simStruct.berRs_packet = berRs_packet;
            simStruct.berRs_gamma = berRs_gamma;
            simStruct.berRs_packet_gamma = berRs_packet_gamma;
            
            simStruct.qq = qq;
            simStruct.mm = mm;
            simStruct.nn_small = nn_small;
            simStruct.kk_small = kk_small;
            simStruct.mm_small = mm_small;
            simStruct.Pbit_lin = Pbit_lin;
            simStruct.cmode = cmode;
            simStruct.bitsPerSymbol = bitsPerSymbol;
            tt=table();
            tt.EbNo = PbitB';
            tt.rsBer = berRs;
            writetable(tt,[exp_path '/results/sim_' num2str(count) '_ber_Q' num2str(Q) 'M' num2str(M) 'n' num2str(nn_small) 'k' num2str(kk_small) '.csv']);
            berCell{qq,mm,iPair} =  simStruct;
            save([exp_path '/results/sim_' num2str(count) '_ber_Q' num2str(Q) 'M' num2str(M) 'n' num2str(nn_small) 'k' num2str(kk_small) '.mat'],'simStruct','Q','M','nn_small','kk_small')
            
            count = count +1;
        end
        
    end
    
    fprintf('Done RS Pair: %d,%d\n',nn_small,kk_small)
end
fprintf('--- Sim Done ---\n')
save([exp_path '/berCell.mat'],'berCell','Qmat','allM','nn_pair','kk_pair','mm_pair','step','PbitB','Pbit_lin','Qrsmat','N','expirement_count','exp_path')
