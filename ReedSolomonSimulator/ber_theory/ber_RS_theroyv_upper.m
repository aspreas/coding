function out =  ber_RS_theroyv_upper(s,m,k,Q)
%n=15
%m=4
%k=7
%
n=(2^m)-1;
q = 2^m;
dmin = n-k+1;
tt = floor(.5*(dmin-1));

distVec = tt+1:n;

% h=m/log2(Q);
% Ps = 1-(1-s).^h;
Ps = s;
out =zeros(length(Ps),1);
 ff= 2^(m-1)/n;
for kk =1:length(Ps)
    tmp = 0;
    for jj = 1:length(distVec)
        ii = distVec(jj);
        tmp =tmp+(ii) *nchoosek(n,ii) * (Ps(kk)^ii) *(1-Ps(kk))^(ii);
    end

    out(kk) =ff*tmp/(n);

end



end