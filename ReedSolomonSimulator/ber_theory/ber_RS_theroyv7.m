function out =  ber_RS_theroyv7(s,m,k,Q)

n=2^m-1;

t = floor((n-k)/2);
tt = t+1:n;

q = 2^m;
h=ceil(m/log2(Q));
Ps = 1-(1-s).^h;
out =zeros(length(Ps),1);
fact = q/(2*(q-1));

for kk =1:length(Ps)
    tmp1 = 0;
    
    for jj = 1:length(tt)
        
        ii = tt(jj);
        %upper bound
        tmp1 =tmp1+ nchoosek(n,ii)* (Ps(kk)^ii) *(1-Ps(kk))^(n-ii);
        %prob        
%         tmp1 =tmp1+ ( (ii+t)/n)*nchoosek(n,ii)* (Ps(kk)^ii) *(1-Ps(kk))^(n-ii);

    end
    out(kk) = fact*(tmp1);
    
end



end