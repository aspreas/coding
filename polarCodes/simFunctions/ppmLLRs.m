function recloglikehood = ppmLLRs(PPMslotblockInput,PPMsymbols, Q, lamda, coefM, zeropos, onespos)
PPMslotblock = reshape(PPMslotblockInput, Q,PPMsymbols)';
PPMslotblock0f1 = zeros(PPMsymbols, Q);
for i = 1:PPMsymbols
    for slotidx = 1:Q
        PPMslotblock0f1(i, slotidx) =  exp(-lamda)*myf1(lamda * PPMslotblock(i,slotidx),coefM);
%         PPMslotblock0f1(i, slotidx) = myf1(PPMslotblock(i,slotidx),coefM);
    end
end

recloglikehood = zeros(PPMsymbols , log2(Q));
for i = 1:PPMsymbols
    for j = 1:log2(Q)
        recloglikehood(i, j) = log(sum(PPMslotblock0f1(i,zeropos(:,j))) / sum(PPMslotblock0f1(i, onespos(:,j))));
    end
end

recloglikehood = reshape(recloglikehood', 1, []);
% recloglikehood = reshape(recloglikehood, [], codelen);
end