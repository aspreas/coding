function [GmatT,databitspos,frozenbitspos,frozenbitchecks,dataLen] = PolarConstr(codeLen,codeRate)
dataLen = round(codeLen * codeRate);

% z  = calculate_channel_polarization( 1, codeLen);
% z = sort(z,'descend');


zz = 1/exp(1);
zz=[2*zz-zz.^2,zz^2];
G=[1, 0 ; 1, 1];
Gmat=G;
for i = 1:log2(codeLen)-1
    zz = [2 * zz - zz.^2, zz.^2];
    zz = reshape([zz(1:end/2); zz(end/2 + 1:end)], [], 1);
    Gmat = kron(Gmat,G);
%     Gmat = [Gmat, zeros(size(Gmat)); Gmat, Gmat];
end
GmatT = Gmat';
[~, sIndx] = sort(zz,'descend');
frozenbitspos = sort(sIndx(1:end - dataLen));
databitspos = sort(sIndx(end-dataLen+1:end));
frozenbitchecks = ones(codeLen,1);
frozenbitchecks(frozenbitspos) = zeros(codeLen-dataLen,1);
end

