function encodedBits = PolarEncoder(databits, codeLen, GmatT, databitspos)
    
    paddedword = zeros(1, codeLen);
    for kk = 1:length(databits)
        paddedword(databitspos(kk)) = databits(kk);
    end
    
    encodedBits = mod(GmatT * paddedword', 2);
end