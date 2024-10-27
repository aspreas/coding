function coef = f1coeffs(M)
f1coeff = 1;
N=1E3;
coef = zeros(N,1);
for kk =1:N
    f1coeff=f1coeff/kk/(kk+M-1);
    coef(kk)=f1coeff;
end
end