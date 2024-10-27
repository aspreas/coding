function [decodedbits, llr] = PolarDec(recvec, id, codelen, frozenbitchecks)

    currentM = length(recvec);
    if currentM == 1
        currentbit = id - codelen + 2;
        decodedbits = frozenbitchecks(currentbit) * round((1 - sign(recvec)) / 2);
        llr = recvec;
        return;
    end
    
    avec = recvec(1:currentM/2);
    bvec = recvec(end-(currentM/2-1):end);

    lvec = sign(avec) .* sign(bvec) .* min(abs(avec), abs(bvec));
    [lbits, lllr] = PolarDec(lvec, 2*id + 1, codelen, frozenbitchecks);
    
    rvec = bvec + avec .* (1 - 2 * lbits);
    [rbits, rllr] = PolarDec(rvec, 2*id + 2, codelen, frozenbitchecks);
    
    decodedbits = [mod(lbits + rbits, 2), rbits];
    llr = [lllr, rllr];
end