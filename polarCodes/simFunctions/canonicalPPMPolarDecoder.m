function decodebitsdemoded = canonicalPPMPolarDecoder(LLR, codeLen, GmatT, databitspos, frozenbitchecks,dataLen,Q)

decodebitsdemoded = zeros(dataLen,1);
LLR_buff = buffer(LLR,log2(Q),0)';
for ii =1:log2(Q)
    decodebitsdemoded_buff = PolarDecoder(LLR_buff(:,ii), codeLen, GmatT, databitspos, frozenbitchecks)';
    decodebitsdemoded(ii:log2(Q):end) = decodebitsdemoded_buff;
end
end