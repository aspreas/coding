
function [BERRC,BERRenc,bitsPerSymbol] = RSPart_V4(Pbit,N,Q,M,mm_small,kk_small,rsEncoder,rsDecoder,cmode,useInterleving)
%% TRANSMITTER
nn_small = 2^mm_small-1;
BERRCMat = zeros(N,1);
BERRencMat = zeros(N,1);
bitsPerSymbol = zeros(N,1);
parfor iBatch=1:N
%     symbols =[];
%     disp(iBatch)
    tranmitted_bits = randi([0 1],mm_small*kk_small,1);
    encodedBits = step(rsEncoder, tranmitted_bits);
    l_enc_bits = mm_small*nn_small;
 
    redundant_bits = ceil( l_enc_bits/mm_small/(log2(Q)) )*mm_small*(log2(Q)) - l_enc_bits;
    encodedBits = [encodedBits; zeros(redundant_bits,1)];
    
    rs_block_tx=(reshape(encodedBits,mm_small,[]));

    if useInterleving
        groupbits = reshape(encodedBits,[],log2(Q));
    else
        groupbits = reshape(encodedBits,log2(Q),[])';
    end
    symbols = zeros(size(groupbits,1)*Q,1);
    for jj =1:size(groupbits,1)
        PPM_position_index=bi2de(groupbits(jj,:),'left-msb');% converting to decimal value
        PPM_position=zeros(Q,1);
        PPM_position(PPM_position_index+1)=1;
%         symbols=[symbols; PPM_position]; % PPM symbol
        symbols(Q*(jj-1)+1: jj*Q) = PPM_position; % PPM symbol
    end
    
    transmitted = symbols;
    %% CHANNEL
    noise = randn(size(transmitted,1),2*M)./sqrt(2*Pbit);
    x = [transmitted zeros(size(transmitted,1),2*M-1)];
    chan = channel(cmode,size(x,1),size(x,2));
    
    hx = x;
    hx2 = zeros(size(hx));
    for ii =1:2*M
        hx2(:,ii) = (abs(hx(:,ii)+noise(:,ii) )).^2;
    end
    y_channel = sum(hx2,2)*chan;
    %% RECEIVER -- soft

    [~,ind] = max(reshape(y_channel,Q,[]));
    received_bits = de2bi(ind-1,log2(Q),'left-msb');

    if useInterleving
        received_bits = received_bits(:);
    else
        received_bits =  reshape(received_bits',[],1);
    end
    
%     BERRencMat(iBatch) = sum(encodedBits~=received_bits)./length(encodedBits);
    rs_block_rx = reshape(received_bits,mm_small,[]);
    BERRencMat(iBatch) = mean(mean(rs_block_tx==rs_block_rx,2));
    bitsPerSymbol(iBatch) = mean(sum(rs_block_tx~=rs_block_rx,1));
    received_bits(end-redundant_bits+1:end)=[];
    
    decodedHard = step(rsDecoder, received_bits);
    BERRCMat(iBatch) = sum(tranmitted_bits~=decodedHard)./length(decodedHard);
end
BERRenc = mean(BERRencMat);
BERRC = mean(BERRCMat);

end
