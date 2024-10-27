% Example call to the main function
Kblock = 1280;
CodeRate = 1/3;
conflvl = 0.99;
conflvlf = 0.9;
Q = 4; % or 16
noisemodes = 40; %2,40,200 to 3
main_simulation(Kblock, CodeRate, conflvl, conflvlf, Q, noisemodes);

function res = Tabmult(A, x)
[matrows, matcols] = size(A);
[vecrows, veccols] = size(x);
if matcols ~= vecrows
    error('Wrong Dimensions in Tabmult!');
end
res = [];
for i = 1:matrows
    sum = zeros(1, veccols);
    for j = 1:matcols
        if A(i, j) >= 0
            sum = sum + circshift(x(j, :), [0, A(i, j)]);
        end
    end
    res = [res; mod(sum, 2)];
end
end

function res = Vecmult(y, x)
[matcols, ~] = size(y);
[vecrows, veccols] = size(x);
if matcols ~= vecrows
    error('Wrong Dimensions in Vecmult!');
end
res = zeros(1, veccols);
for j = 1:matcols
    if y(j) >= 0
        res = res + circshift(x(j, :), [0, y(j)]);
    end
end
res = mod(res, 2);
end

function codestream = LDPCEncoder(bitstream,Zc,CoreCheckmat,CorePartitymat,ExtendedCheckmat,checkcol,coderows,codecols,datacols,BG)

pbitstream = reshape(bitstream, [], Zc);

% Find p1
pvec = Tabmult(CoreCheckmat, pbitstream);
p = circshift(mod(sum(pvec, 1), 2), [0, checkcol]);
codestream = [pbitstream; p];

% Find p2, p3, p4
for i = 1:3
    newvec = [CoreCheckmat(i, :), CorePartitymat(i,1:i)]';
    p = Vecmult(newvec, codestream);
    codestream = [codestream; p];
end

% Find p5 ...
if coderows > 4
    pvec = Tabmult(ExtendedCheckmat, codestream);
    codestream = [codestream; pvec];
end
end

function decodebits = LDPCDecoder(bitstream,Hmat,Zc,coderows,codecols,datacols,BG)

tries = 0;
decodebits = bitstream;
checkpath = {};

syndrome = Tabmult(Hmat, decodebits);
oldsyndrome = mod(syndrome + syndrome, 2);
synpos = [];
totalTallyunf = {};
while sum(syndrome(:)) ~= 0 && tries < 100
    tries = tries + 1;
    
    % Remove Syndrome Positions
    delta = syndrome - oldsyndrome;
    deltaminus = find(delta == -1);
    elems = length(deltaminus);
    j = 1;
    k = 1;
    while j <= elems
        if deltaminus(j) == synpos(k)
            synpos(k) = [];
            totalTallyunf(k) = [];
            j = j + 1;
        else
            k = k + 1;
        end
    end
    
    % Add Syndrome Positions
    deltaplus = find(delta == 1);
    elems = length(deltaplus);
    i = 1;
    k = 1;
    while i <= elems
        if k > length(synpos)
            synpos = [synpos(1:k-1), deltaplus(i), synpos(k:end)];
            totalTallyunf = [totalTallyunf(1:k-1), ...
                arrayfun(@(j) ...
                ifelse(Hmat(synpos(k)) ~= -1, ...
                (j-1)*Zc + mod(Hmat(synpos(k)) + synpos(i) - 1, Zc), ...
                []), ...
                1:codecols), ...
                totalTallyunf(k:end)];
            k = k + 1;
            i = i + 1;
            continue;
        end
        while k <= length(synpos) && deltaplus(i) > synpos(k)
            k = k + 1;
        end
        while k <= length(synpos) && deltaplus(i) == synpos(k) && deltaplus(i) > synpos(k)
            k = k + 1;
        end
        synpos = [synpos(1:k-1), deltaplus(i), synpos(k:end)];
        totalTallyunf = [totalTallyunf(1:k-1), ...
            arrayfun(@(j) ...
            ifelse(Hmat(synpos(k)) ~= -1, ...
            (j-1)*Zc + mod(Hmat(synpos(k)) + synpos(i) - 1, Zc), ...
            []), ...
            1:codecols), ...
            totalTallyunf(k:end)];
        k = k + 1;
        i = i + 1;
    end
    
    totalTally = cell2mat(totalTallyunf);
    chosencol = 0;
    while chosencol == 0
        maxpos = datasample(mode(totalTally), 1);
        synerrors = [];
        for k = 1:length(maxpos)
            testdecodebits = decodebits;
            bitblock = 1 + floor((maxpos(k)-1) / Zc);
            bitpos = 1 + mod(maxpos(k)-1, Zc);
            testdecodebits(bitblock, bitpos) = mod(1 + testdecodebits(bitblock, bitpos), 2);
            
            syndrome = Tabmult(Hmat, testdecodebits);
            synerrorstmp = sum(syndrome(:));
            
            if synerrorstmp == 0
                decodebits = testdecodebits;
                return;
            end
            
            if syndrome(1) >= 5
                decodebits = testdecodebits;
                return;
            end
            
            if any(ismember(checkpath, [maxpos(k), synerrorstmp]))
                totalTally(totalTally == maxpos(k)) = [];
                continue;
            end
            
            synerrors = [synerrors, synerrorstmp];
            if synerrorstmp == min(synerrors)
                savedecodebits = testdecodebits;
                chosencol = k;
            end
        end
    end
    
    if chosencol > 0
        oldsyndrome = Tabmult(Hmat, decodebits);
        decodebits = savedecodebits;
        syndrome = Tabmult(Hmat, decodebits);
        checkpath = [checkpath; maxpos(chosencol), sum(syndrome(:))];
    end
end

decodebits = decodebits(:, 1:Zc*datacols(BG));
end

function perr = PPMperr(OSNR, Q, noisemodes)
cnq = importdata(sprintf('cnq%d.txt', noisemodes));
sum = 0;
for q = 1:Q-1
    currsum = 0;
    lagpp = 0;
    lagp = 0;
    lagn = 1;
    qn = 1 / (1 + q)^noisemodes;
    for n = 0:q*(noisemodes-1)
        tempsum = cnq(q+1,n+1) * qn * lagn;
        lagpp = lagp;
        lagp = lagn;
        lagn = (2*n + 1 + noisemodes - 1 + OSNR/(1+q)) / (n+1) * lagp - (n + noisemodes - 1) / (n+1) * lagpp;
        qn = qn / (1+q);
        currsum = currsum + tempsum;
    end
    newsum = currsum * nchoosek(Q-1, q) * exp(-OSNR*q/(1+q)) * (-1)^(q+1);
    sum = sum + newsum;
end
perr = Q / 2 / (Q - 1) * sum;
end

function span = FinddBspan(Q, noisemodes, BEPstart, BEPend)
i = 0;
tarperr = BEPstart;
while PPMperr(10^(i/10) * log2(Q), Q, noisemodes) > tarperr
    i = i + 1;
end
while PPMperr(10^(i/10) * log2(Q), Q, noisemodes) < tarperr
    i = i - 0.1;
end
dBstart = i;
tarperr = BEPend;
while PPMperr(10^(i/10) * log2(Q), Q, noisemodes) > tarperr
    i = i + 1;
end
while PPMperr(10^(i/10) * log2(Q), Q, noisemodes) < tarperr
    i = i - 0.1;
end
dBend = i;
span = [dBstart, dBend];
end




function sum = My0f1(z,f1coefflist)
sum = 1;
zpow = 1;
newterm = 1;
i = 1;

while newterm/sum > 1e-6
    zpow = zpow * z;
    newterm = f1coefflist(i) * zpow;
    sum = sum + newterm;
    i = i + 1;
end
end

function newvec = minsumfun(currvec)
absvec = abs(currvec);
pos = sort(absvec, 'ascend');
newvec = sign(currvec) * (1 - 2 * mod(sum(1 - double(currvec > 0)), 2)) * absvec(pos(1));
newvec(pos(1)) = newvec(pos(1)) * absvec(pos(2)) / absvec(pos(1));
end

function decodebits = LDPCsumprodDecoder(likelihoodvec,Hmat,Zc,coderows,codecols,datacols,BG,totalrowvec)

if length(likelihoodvec) ~= Zc*codecols
    error('Error: Wrong Dimensions in the Likelihood Vector!');
end

decodebits = likelihoodvec;
decodebits(decodebits > 0) = 0;
decodebits(decodebits < 0) = 1;

syndrome = Tabmult(Hmat, reshape(decodebits, [], Zc));
if sum(syndrome(:)) == 0
    decodebits = decodebits(1:Zc*datacols(BG));
    return;
end

Lmat = reshape(likelihoodvec(totalrowvec), Zc*coderows, []);
for tries = 1:10
    colsum = likelihoodvec;
    for i = 1:Zc*coderows
        currvec = Lmat(i, :);
        fvec = tanh(currvec / 2);
        fprod = prod(fvec);
        newvec = 2 * atanh(fprod ./ fvec);
        posvec = find(isnan(newvec));
        for idx = 1:length(posvec)
            pos = posvec(idx);
            newvec(pos) = min(abs(currvec([1:pos-1 pos+1:end]))) * prod(sign(currvec)) * sign(currvec(pos));
        end
        Lmat(i, :) = newvec;
        colsum(totalrowvec(i)) = colsum(totalrowvec(i)) + Lmat(i, :);
    end
    Lmat = reshape(colsum(totalrowvec)' - Lmat, Zc*coderows, []);
    
    decodebits = colsum;
    decodebits(decodebits > 0) = 0;
    decodebits(decodebits < 0) = 1;
    
    syndrome = Tabmult(Hmat, reshape(decodebits, [], Zc));
    if sum(syndrome(:)) == 0
        decodebits = decodebits(1:Zc*datacols(BG));
        return;
    end
end

decodebits = decodebits(1:Zc*datacols(BG));
end

function main_simulation(Kblock, CodeRate, conflvl, conflvlf, Q, noisemodes)
%     global Zc CoreCheckmat CorePartitymat ExtendedCheckmat checkcol coderows codecols datacols BG totalrowvec;

% Step 1: Find Base Graph
BG = 1;
if Kblock <= 292
    BG = 2;
end
if Kblock <= 3824 && CodeRate <= 0.67
    BG = 2;
end
if CodeRate <= 0.25
    BG = 2;
end

% Step 2: Find Kb
if BG == 1
    Kb = 22;
else
    if Kblock > 640
        Kb = 10;
    elseif Kblock > 560 && Kblock <= 640
        Kb = 9;
    elseif Kblock > 192 && Kblock <= 560
        Kb = 8;
    else
        Kb = 6;
    end
end

% Step 3: Find Zc
avals = [2, 3, 5, 7, 9, 11, 13, 15];
Zcvals = arrayfun(@(j) avals * 2^j, 0:7, 'UniformOutput', false);
Zcvals_flat = sort(cell2mat(Zcvals(:)));
Zc = Zcvals_flat(find(Zcvals_flat >= (Kblock/Kb), 1));
if Zc > 384 || Zc < 2
    error('ERROR: Zc Out of Bounds!');
end

% Step 4: Find iLS
[~, idx] = find(cellfun(@(x) ismember(Zc, x), Zcvals), 1);
ils = idx - 1;
fprintf('BG: %d Kb: %d Zc: %d iLS: %d\n', BG, Kb, Zc, ils);

% Step 5: Calculate Base Matrix
filename = sprintf('BG%d.csv', BG);
BGmatvals = readmatrix(filename);
len = size(BGmatvals, 1);

BGrowcol = [46, 68; 42, 52];
datacols = [22, 10];
rows = BGrowcol(BG, 1);
cols = BGrowcol(BG, 2);

BGmat = -ones(rows, cols);

for k = 1:len
    i = BGmatvals(k, 1);
    j = BGmatvals(k, 2);
    v = BGmatvals(k, ils + 3);
    BGmat(i + 1, j + 1) = mod(v, Zc);
end

% Shorten the Matrix According to Coderate
codecols = ceil(datacols(BG) / CodeRate);
coderows = codecols - datacols(BG);
if codecols > BGrowcol(BG, 2)
    error('Code Rate Requires Too Many Columns!');
end
if codecols < datacols(BG) + 4
    error('Code Rate Requires Too Few Columns!');
end
if coderows > BGrowcol(BG, 1)
    error('Code Rate Requires Too Many Rows!');
end
if coderows < 4
    error('Code Rate Requires Too Few Rows!');
end

% Check the Core Parity Matrix, Block Zero Matrix and Block Identity Matrix
Blockzeromat = BGmat(1:4, datacols(BG) + 5:codecols);
if nnz(Blockzeromat == -1) ~= 4 * (codecols - datacols(BG) - 4)
    error('ERROR: Block Zero Matrix has the Wrong Form!');
end

Identitymat = BGmat(5:coderows, datacols(BG) + 5:codecols);
if nnz(diag(Identitymat) == 0) ~= (codecols - datacols(BG) - 4)
    error('ERROR: Block Identity Matrix has the Wrong Form!');
end
if nnz(Identitymat == -1) ~= (coderows - 4) * (codecols - datacols(BG) - 5)
    error('ERROR: Block Identity Matrix has the Wrong Form!');
end

CorePartitymat = BGmat(1:4, datacols(BG) + 1:datacols(BG) + 4);
for j = 2:4
    checkcol = CorePartitymat(:, j);
    checkcol(checkcol == -1) = [];
    if j > 1 && (numel(checkcol) ~= 2 || nnz(checkcol == 0) ~= 2)
        error('ERROR: Core Parity Matrix has the Wrong Form!');
    end
end

checkcol = CorePartitymat(:, 1);
checkcol(checkcol == -1) = [];
remove_indx =[];
for i = 1:numel(checkcol)
    cnt = nnz(checkcol == checkcol(i));
    if mod(cnt, 2) == 0
        remove_indx = [remove_indx;i];
    end
end
checkcol(remove_indx)=[];
if numel(checkcol) ~= 1
    error('ERROR: Core Parity Matrix has the Wrong Form!');
end

checkcol = checkcol(1);
CoreCheckmat = BGmat(1:4, 1:datacols(BG));
ExtendedCheckmat = BGmat(5:coderows, 1:datacols(BG) + 4);
Hmat = BGmat(1:coderows, 1:codecols);

% Find the non-zero values that are required for row/column operations at the min-sum decoder
totalrowvec = zeros(coderows * Zc, 1);
index = 1;
for i = 1:coderows
    rowvec = [];
    for j = 1:codecols
        if Hmat(i, j) >= 0
            newrowvec = Zc * (j - 1) + circshift(1:Zc, [0, Hmat(i, j)]);
            rowvec = [rowvec, newrowvec(:)']; % Append newrowvec as rows to rowvec
        end
    end
    for k = 1:Zc
        totalrowvec(index) = rowvec(k);
        index = index + 1;
    end
end


totalcolvec = zeros(codecols * Zc, 1);
index = 1;

for j = 1:codecols
    colvec = [];
    for i = 1:coderows
        if Hmat(i, j) >= 0
            newcolvec = Zc * (i - 1) + circshift(1:Zc, [0, -Hmat(i, j)]);
            colvec = [colvec; newcolvec(:)];
        end
    end
    for k = 1:Zc
        totalcolvec(index) = colvec(k);
        index = index + 1;
    end
end


disp(CoreCheckmat);
disp(ExtendedCheckmat);
disp(BGmat);
disp(Hmat);
fprintf('%d\n', Zc * codecols);
fprintf('%f\n', datacols(BG) / codecols);

% Simulation with exact likelihood

% Confidence intervals for the simulation
alphaconf = 1 - conflvl;
alphaconff = 1 - conflvlf;

if Kblock < 4576
    dBstep = 0.1;
else
    dBstep = 0.05;
end


% For the calculation of the 0f1 hypergeometric - PPM symbol likelihood
f1coeff = 1;
f1coefflist = 1;
for i = 1:1000
    f1coeff = f1coeff / i / (i + noisemodes - 1);
    f1coefflist = [f1coefflist, f1coeff];
end

h = log2(Q);
codelen = Zc * codecols;
if mod(codelen, h) ~= 0
    error('Codelen is not divisible by h');
end
PPMsymbols = codelen / h;
BEPresults = {};
FEPresults = {};

% '0' and '1' sets of each bit
mat = de2bi(0:(Q - 1), 'left-msb', log2(Q));
zeropos = cell(1, h);
onepos = cell(1, h);
for i = 1:h
    tab = mat(:, i);
    zeropos{i} = find(tab == 0);
    onepos{i} = find(tab == 1);
end

fprintf('Running for: Q-%d and M-%d\n', Q, noisemodes);

dBspan = FinddBspan(Q, noisemodes, 0.1, 1e-5);

for EbNo = dBspan(1):dBstep:dBspan(2)
    fprintf('Running for Eb/No: %.1f\n', EbNo);
    
    lamda = 10^(EbNo / 10) * log2(Q) * CodeRate;
    totaldataerrors = 0;
    totalerrors = 0;
    frameerrors = 0;
    h1 = 1;
    h2 = 0;
    h1f = 1;
    h2f = 0;
    
    transmissions = 0;
    while (h1 - h2) * transmissions * Kblock >= 0.1 * totaldataerrors || ...
            (h1f - h2f) * transmissions >= 0.5 * frameerrors
        
        % Generate the LDPC code words - codewords are processed in parallel
        databits = zeros(1, Zc * datacols(BG));
        for ker = 1:1
            databits(ker, :) = [randi([0 1], 1, Kblock), zeros(1, Zc * datacols(BG) - Kblock)];
        end
        codebits = LDPCEncoder(databits,Zc,CoreCheckmat,CorePartitymat,ExtendedCheckmat,checkcol,coderows,codecols,datacols,BG);
        
        % Energy slot positions in symbols
        groupedcodebits = reshape(codebits, h, []);
        energyslotpositions = bi2de(groupedcodebits', 'left-msb') + 1;
        
        % Chi-square random variables
        PPMslotblock = reshape(chi2rnd(2 * noisemodes, Q, PPMsymbols * 1) / 2, Q, []);
        onepowers = ncx2rnd(2 * noisemodes, 2 * lamda, 1, PPMsymbols * 1) / 2;
        for i = 1:length(energyslotpositions)
            PPMslotblock(energyslotpositions(i), i) = onepowers(i);
        end
        
        % Calculation of PPM slot likelihoods and bit log-likelihoods
        PPMslotblock0f1 = exp(-lamda) * My0f1(lamda * PPMslotblock,f1coefflist);
        recloglikehood = zeros(PPMsymbols , h);
        for j = 1:h
            recloglikehood(:, j) = log(sum(PPMslotblock0f1(zeropos{j}, :))./sum(PPMslotblock0f1(onepos{j}, :)));
        end
        recloglikehood = reshape(recloglikehood', 1, []);
        
        decodedbits = LDPCsumprodDecoder(recloglikehood,Hmat,Zc,coderows,codecols,datacols,BG,totalrowvec);
        
        totaldataerrors = totaldataerrors + sum(sum(mod(databits + decodedbits, 2)));
        
        currenterrors = sum(mod(databits + decodedbits, 2), 2);
        frameerrors = frameerrors + sum(currenterrors > 0);
        
        receivedcodebits = reshape(recloglikehood, [], codelen);
        receivedcodebits(receivedcodebits > 0) = 0;
        receivedcodebits(receivedcodebits < 0) = 1;
        totalerrors = totalerrors + sum(sum(mod(codebits + receivedcodebits, 2)));
        
        transmissions = transmissions + 1;
        if totaldataerrors == Inf
            fprintf('ERROR! %d\n', totaldataerrors);
        end
        
        % Calculation of the confidence intervals
        if transmissions * Kblock > totaldataerrors
            h1 = betaincinv(1 - alphaconf / 2, totaldataerrors + 1, transmissions * Kblock - totaldataerrors, 'upper');
        end
        if totaldataerrors > 0
            h2 = betaincinv(alphaconf / 2, totaldataerrors, transmissions * Kblock - totaldataerrors + 1, 'upper');
        end
        
        if transmissions > frameerrors
            h1f = betaincinv(1 - alphaconff / 2, frameerrors + 1, transmissions - frameerrors, 'upper');
        end
        if frameerrors > 0
            h2f = betaincinv(alphaconff / 2, frameerrors, transmissions - frameerrors + 1, 'upper');
        end
    end
    
    fprintf('Regular sampling required: %d transmissions\n', transmissions);
    fprintf('Coded BEP: %f %f %f\n', h2, totaldataerrors / Kblock / transmissions, h1);
    fprintf('Uncoded BEP: %f\n', totalerrors / codelen / transmissions);
    fprintf('Coded FEP: %f %f %f\n', h2f, frameerrors / transmissions, h1f);
    
    BEPresults{end + 1} = [EbNo, totaldataerrors / Kblock / transmissions, h2, h1];
    FEPresults{end + 1} = [EbNo, frameerrors / transmissions, h2f, h1f];
    if totaldataerrors / Kblock / transmissions < 1e-5
        break;
    end
end

filename = sprintf('0f1-BEP-K-%d-Q%d-M%d.txt', Kblock, Q, noisemodes);
writematrix(cell2mat(BEPresults), filename, 'Delimiter', 'tab');

filename = sprintf('0f1-FEP-K-%d-Q%d-M%d.txt', Kblock, Q, noisemodes);
writematrix(cell2mat(FEPresults), filename, 'Delimiter', 'tab');
end


