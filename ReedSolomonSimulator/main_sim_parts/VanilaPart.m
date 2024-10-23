function BERVanila = VanilaPart(Pbit,N,Q,M,mm_small,kk_small,cmode)


%% TRANSMITTER
symbols =[];
tranmitted_bits = randi([0 1],N*log2(Q)*mm_small*kk_small,1);
groupbits = reshape(tranmitted_bits,[],log2(Q));

NN = size(groupbits,1);

for jj=1:NN
    PPM_position_index=bi2de(groupbits(jj,:),'left-msb');% converting to decimal value
    PPM_position=zeros(Q,1);
    PPM_position(PPM_position_index+1)=1;
    symbols=[symbols; PPM_position]; % PPM symbol 
end
transmitted = symbols;
%% CHANNEL
noise = randn(size(transmitted,1),2*M)./sqrt(2*Pbit);
x = [transmitted zeros(size(transmitted,1),2*M-1)];
chan = channel(cmode,size(x,1),size(x,2));

hx = x;
hx2 =zeros(size(hx));
for ii =1:2*M
    hx2(:,ii) = (abs(hx(:,ii)+noise(:,ii) )).^2;
end
y_channel = sum(hx2.*chan,2);
%% RECEIVER -- soft
[~,ind] = max(reshape(y_channel,Q,[]));
received_bits = de2bi(ind-1,'left-msb');
received_bits = reshape(received_bits,1,[]);

BERVanila= sum(tranmitted_bits~=received_bits')./length(received_bits);
