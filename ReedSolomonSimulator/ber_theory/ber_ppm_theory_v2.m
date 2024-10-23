function [Pe,Ser] = ber_ppm_theory_v2(EbNodB,Q,M,pp)

% EbNodB: the range of snrs in dB
% Q: modulation order
% M: noise modes

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


Pe = zeros(length(EbNodB),1);

parfor EbNo_db = 1:length(EbNodB)
    
    EbNo = 10^(EbNodB(EbNo_db)/10);
    ll = EbNo;
    tmp3= zeros(Q-1,1);
    for q=1:Q-1
        tmp2= zeros(q*(M-1),1);
        for n = 0:q*(M-1)
            tmp = zeros(n+1,1);
            for i=0:n
                s3 = factorial(n+M-1)/(factorial(i+M-1)*factorial(n-i)*factorial(i));
                tmp(i+1) = s3*(ll/(1+q))^i;   
            end
            ss1  = sum(abs(tmp));
            tmp2(n+1) = ss1*c(q+1,n+1)/((1+q)^(n+M));
        end
        ss2 = sum(abs(tmp2));
        tmp3(q) = ss2*( nchoosek(Q-1,q)* ((-1)^(q+1)) )*exp(-(ll*q)/(1+q)) ;
    end
    ss3 = sum((tmp3));
    Pe(EbNo_db) = (Q/(2*(Q-1)))*ss3;
    Ser(EbNo_db) =ss3;
end



end
