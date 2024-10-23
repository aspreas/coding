addpath(genpath(pwd))
fig22 = figure('units','normalized','outerposition',[0 0 1 1]);
legendS ={}; count=1;countM=1;

styleSim         ='-b*';
styleTheory      ='-rd';
styleSimRS       ='-k*';
styleSimRSG      ='-ro';
styleTheoryRS    ='-mo';
styleTheoryRSUPP ='-g*';

all_marks = {'o','+','*','.','x','s','d','^','v','>','<','p','h'};

for iPair = 1:length(kk_pair)
    
    mm_small = mm_pair(iPair);
    kk_small = kk_pair(iPair);
    nn_small = nn_pair(iPair);
    codingrate = kk_small/nn_small;
    
    for mm = 1:length(allM)
        M = allM(mm);
        for qq = 1:length(Qmat)
            Q = Qmat(qq);
            Qrs = Qrsmat(qq);
            
            berRs = berCell{qq,mm,iPair}.berRs;
            berRs_gamma = berCell{qq,mm,iPair}.berRs_gamma;
            berVanila = berCell{qq,mm,iPair}.berVanila;
            berRs_inter = berCell{qq,mm,iPair}.berRs_inter;
            
            berRs_packet = berCell{qq,mm,iPair}.berRs_packet;
            berRs_packet_inter = berCell{qq,mm,iPair}.berRs_packet_inter;
            
            
            [Pe2,Ser]     = ber_ppm_theory_v2(10*log10(Pbit_lin.*log2(Q)),Q,M,cnq_path);
            [Pe22,Ser2]   = ber_ppm_theory_v2(10*log10(Pbit_lin*log2(Qrs)*codingrate),Qrs,M,cnq_path);            

            Per_RS_Theory = ber_RS_theroyv5(Ser2,mm_small,kk_small,Qrs);
            Per_RS_Theory2 = ber_RS_theroyv7(Ser2,mm_small,kk_small,Qrs);

            fig = figure('units','normalized','outerposition',[0 0 1 1]);
            set(fig, 'Visible', 'off');
            semilogy(PbitB,berRs,styleSimRS)
            hold on
            semilogy(PbitB,berRs_inter,'-x')
            semilogy(PbitB,berVanila,styleSim)
            semilogy(PbitB,Per_RS_Theory,styleTheoryRS)
            semilogy(PbitB,Pe2,styleTheory)
            semilogy(PbitB,min(Pe2,Per_RS_Theory),'p-')
            semilogy(PbitB,Per_RS_Theory2,'d-')
            semilogy(PbitB,berRs_packet,'>-')
            semilogy(PbitB,berRs_packet_inter,'<-')
            
            legendStyle = cell(7,1);
            legendStyle{1} =['Sim RS no inter Q = ' num2str(Qrs) ' - M = ' num2str(M)];
            legendStyle{2} =['Sim RS inter Q = ' num2str(Qrs) ' - M = ' num2str(M)];

            legendStyle{3} =['Sim Q = ' num2str(Q) ' - M = ' num2str(M)];
            legendStyle{4} =['Theory RS (' num2str(nn_small) ',' num2str(kk_small)  ') Q = ' num2str(Qrs) ' - M = ' num2str(M)];
            legendStyle{5} =['Theory Q = ' num2str(Q) ' - M = ' num2str(M)];
            legendStyle{6} =['Two folded RS (' num2str(nn_small) ',' num2str(kk_small)  ') Q = ' num2str(Qrs) ' - M = ' num2str(M)];
            legendStyle{7} =['Theory RS Bound (' num2str(nn_small) ',' num2str(kk_small)  ') Q = ' num2str(Qrs) ' - M = ' num2str(M)];
            legendStyle{8} =['Sim RS packer no inter Q = ' num2str(Qrs) ' - M = ' num2str(M)];
            legendStyle{9} =['Sim RS packet inter Q = ' num2str(Qrs) ' - M = ' num2str(M)];
            
            xlabel('EbNo')
            ylabel('Bit Error Rate')
            title(['Number of bits = ', num2str(N*nn_small*nn_small,'%.e') ,', ' num2str(Q), '-PPM, M = ' num2str(M) ])
            legend (legendStyle);
            grid on
            axis([min(PbitB) max(PbitB) 1e-6 1])
            xticks(PbitB(1:1/step:end))
            set(fig, 'Visible', 'on');
            
            
            figname=[exp_path '/figs/Q' num2str(Q) '-M' num2str(M) 'rs(' num2str(nn_small) ',' num2str(kk_small) ').jpg'];
            fignamefig=[exp_path '/figs/Q' num2str(Q) '-M' num2str(M) 'rs(' num2str(nn_small) ',' num2str(kk_small) ').fig'];
            
            saveas(fig,figname)
            saveas(fig,fignamefig)
            
            
        end
    end
end


figname=[exp_path '/figs/Q' num2str(Q) '-M' num2str(M) 'rs(' num2str(nn_small) ',' num2str(kk_small) ')_distance.jpg'];
fignamefig=[exp_path '/figs/Q' num2str(Q) '-M' num2str(M) 'rs(' num2str(nn_small) ',' num2str(kk_small) ')_distance.fig'];

saveas(fig,figname)
saveas(fig,fignamefig)