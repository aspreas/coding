function encodedBits = canonicalPPMPolarEncoder(databits,codeLen, GmatT, databitspos,Q)
encodedBits=zeros(codeLen*log2(Q),1);
databits_par = buffer(databits,log2(Q),0)';
for ii =1:log2(Q)
    encodedBits_buff = PolarEncoder(databits_par(:,ii),codeLen, GmatT, databitspos);
    encodedBits(ii:log2(Q):end) = encodedBits_buff;
end

end