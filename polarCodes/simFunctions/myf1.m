function out = myf1(input,f1coefs)
ss=1;zpow=1;newterm=1;
kk =1;
while newterm/ss > 10^-6
    zpow=zpow*input;
    newterm=f1coefs(kk)*zpow;
    if isinf(newterm) 
        newterm=1;
    end
    ss=ss+newterm;
    kk=kk+1;
end
out = ss;
end
