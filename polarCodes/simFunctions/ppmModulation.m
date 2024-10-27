function [symbols,zeroPos,onePos] = ppmModulation(Q,encodedBits)
groupbits = reshape(encodedBits,log2(Q),[])';
symbols = zeros(size(groupbits,1)*Q,1);
for jj =1:size(groupbits,1)
    PPM_position_index=bi2de(groupbits(jj,:),'left-msb');% converting to decimal value
%    PPM_position_index=bi2de(groupbits(jj,:));% converting to decimal value

    PPM_position=zeros(Q,1);
    PPM_position(PPM_position_index+1)=1;
    symbols(Q*(jj-1)+1: jj*Q) = PPM_position; % PPM symbol
end
bits = de2bi(0:Q-1,log2(Q),'left-msb');
%bits = de2bi(0:Q-1,log2(Q));

onePos=zeros(Q/2,log2(Q));
zeroPos=zeros(Q/2,log2(Q));
for kk =1:log2(Q)
    onePos(:,kk)=find(bits(:,kk)) ;
    zeroPos(:,kk)=find(bits(:,kk)==0);
end
