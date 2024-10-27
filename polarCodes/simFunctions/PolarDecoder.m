function decodedBits = PolarDecoder(recvec, codelen, GmatT, databitspos, frozenbitchecks)
    [decodedbits, llr] = PolarDec(recvec, 0, codelen, frozenbitchecks);
    decodedBits = mod(GmatT * decodedbits', 2)';
    decodedBits = decodedBits(databitspos);
end