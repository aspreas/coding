function [out,c] = calculate_ber_stable_terms(Q,M)
out= {};
for q=1:Q-1
    for n = 0:q*(M-1)
        tmp = [];
        for i=0:n
            s3 = factorial(n+M-1)/(factorial(i+M-1)*factorial(n-i)*factorial(i));
            tmp = [tmp;s3];
        end
        out{q,n+1} = tmp;
    end
end

pp = 'C:\Users\Antonis\Desktop\myMatlabFiles\mathematica\cnqFiles\cnq';
file = [pp num2str(M) num2str(Q) '.txt' ];
if ~exist(file,'file')==2
    error('No file %s', file)
    return;
end

fid = fopen(file);
line1 = fgetl(fid);
llfile=length(find(line1=='.'));
formating = repmat('%f',1,llfile);
fclose(fid);


fid = fopen( file);
data = textscan(fid, formating, 'Delimiter', '.	');
fclose(fid);

c = zeros(Q,llfile);
for mm =1:llfile
    c(:,mm) = data{mm};
    c(isinf(c(:,mm)),mm) = 0;
end





end