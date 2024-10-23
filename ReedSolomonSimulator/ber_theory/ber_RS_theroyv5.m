function out =  ber_RS_theroyv5(s,m,k,Q)
%n=15
%m=4
%k=7
%

n=2^m-1;

t = floor((n-k)/2);
d=n-k+1;

rrvec = d+1:n;
% tt = t+1:d;
tt = t+1:n;

q = 2^m;
 h=ceil(m/log2(Q));
 Ps = 1-(1-s).^h;
% Ps = s;
out =zeros(length(Ps),1);
fact = q/(2*(q-1));
for kk =1:length(Ps)
    tmp1 = 0; tmp2 = 0;
%     for jj = 1:length(tt)
%         ii = tt(jj);
%         tmp1 =tmp1+(d/n) *nchoosek(n,ii) * (Ps(kk)^ii) *(1-Ps(kk))^(n-ii);
%     end
%     for jj = 1:length(rrvec)
%         ii = rrvec(jj);
%         tmp2 =tmp2+ nchoosek(n-1,ii-1) * (Ps(kk)^ii) *(1-Ps(kk))^(n-ii);
%     end
%   out(kk) = fact*(tmp1+tmp2);  
  
  for jj = 1:length(tt)
        ii = tt(jj);
        tmp1 =tmp1+nchoosek(n-1,ii-1) * (Ps(kk)^ii) *(1-Ps(kk))^(n-ii);
    end
  out(kk) = fact*(tmp1);

end



end