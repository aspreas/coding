clear;close all;clc;
addpath(genpath(pwd))
warning off;
expirement_count=16;
config_sim;
general_counter = 0;
resPath='C:/Users/Antonis/Desktop/myMatlabFiles/rsResults/exp_16/results';
d = dir(resPath);
d(1:2)=[];
structNames={d.name};

for qq = 1:length(Qmat)
    Q =Qmat(qq);
    legendStyle={};
    for mm = 1:length(allM)
        M=allM(mm);
                fig = figure('units','normalized','outerposition',[0 0 1 1]);
                for i_rs_pair =1:length(nn_pair)
                    general_counter = general_counter+1;
                    kk_small = kk_pair(i_rs_pair);
                    nn_small = nn_pair(i_rs_pair);
                 
                    
                    currentName = ['ber_Q' num2str(Q) ...
                        'M' num2str(M) ...
                        'n' num2str(nn_small) 'k' num2str(kk_small) '.mat'];
                    indx = find(contains(structNames,currentName));
                    loadPathStruct = [ resPath '/' structNames{indx}];
                    %     keyboard
                    %                     if general_counter~=2
                    load(loadPathStruct)
                    %%
                    
                    %set(fig, 'Visible', 'off');
                    PbitB = 10*log10(simStruct.Pbit_lin);
                    
                    if i_rs_pair==1
                        semilogy(PbitB,simStruct.berVanila,'-k','Linewidth',1.1)
                        legendStyle{end+1} =['Q = ' num2str(Q) ' - M = ' num2str( M)];
                    end
                        hold on
                        semilogy(PbitB,simStruct.berRs,'-d','Linewidth',1.1)
                        
                        legendStyle{end+1} =['RS Q = ' num2str(Q) ' - M = ' num2str( M) 'RS(' num2str(nn_small) ',' num2str(kk_small) ')'];
                      
                  
                end
                axis([0 28 1e-6 1])
                legend (legendStyle);
                grid on
                xlabel('Eb/No')
                ylabel ('Bit Error Rate')
                legendStyle={};
                title(['Q=' num2str(Q)  ', M=' num2str(M) ])
%                 saveas(fig,filename,'jpg')
%                 saveas(fig,filename,'fig')
            end
            % %%
        end

