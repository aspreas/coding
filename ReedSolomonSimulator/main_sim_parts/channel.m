function chan = channel(cmode,L,n,severity)

switch cmode
    case 1
        chan = 1;
    case 2
        a = 4.04;
        b = 1.53;
        chan = gamrnd(a,1,L,n).*gamrnd(b,1,L,n)/(a*b*n^2);
    case 3
        switch severity
            case 'weak'
                a=50;
                b = 14;
                omega=1.0621;
                beta0 = 0.0216;
                rho = 0.86; 
                deltaphi=0;
            case'moderate'
                a=2.55;
                b = 22;
                omega=0.4618;
                beta0 = 0.6525;
                rho = 0.988; 
                deltaphi=pi/2;
            case'strong'
                a=2.2814;
                b = 33;
                omega=1.33;
                beta0 = 0.4231;
                rho = 0.84; 
                deltaphi=0;
        end
        gam = 2*beta0*(1-rho);
        omegaprime = omega+rho*2*beta0 +2*sqrt(2*beta0*rho)*cos(deltaphi);

        g = sqrt(gamrnd(b,1,L,n)/b);
        USprime = sqrt(2*beta0)*randn(L,n);
        phi = rand(L,n)*2*pi;
        r = g.*(sqrt(omega)+sqrt(2*beta0*rho).*exp(1i*deltaphi) +sqrt(1-rho).*USprime.*exp(1i*phi) );
        
        y = abs(r.^2);
        x = gamrnd(a,1,L,n)/a;
        irradiance = x.*y;
        f_mean=malaga_moment(1,a,b,gam,omegaprime);
        f_std=malaga_moment(2,a,b,gam,omegaprime);

        chan = (irradiance /f_mean );
end




end