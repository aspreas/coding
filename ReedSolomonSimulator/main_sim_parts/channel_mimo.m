function chan = channel_mimo(cmode,L,n,z,severity)

switch cmode
    case 1
        chan = 1;
    case 2
        switch severity
            case 'weak'
                %1-100m
                a = 16.5347;
                b = 14.9057;
            case 'strong'
                % 1-1000m
                a = 5.50966;
                b = 1.1138;
        end
        
        chan = gamrnd(a,1,L,n,z).*gamrnd(b,1,L,n,z)/(a*b*n^2);
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
        A = (2*(a^(a/2)))/ ((gam^(1+a/2))*gamma(a)) * ((gam*b)/(gam*b+omegaprime))^(b+a/2)
        g = sqrt(gamrnd(b,1,L,n,z)/b);
        USprime = sqrt(2*beta0)*randn(L,n,z);
        phi = rand(L,n,z)*2*pi;
        r = g.*(sqrt(omega)+sqrt(2*beta0*rho).*exp(1i*deltaphi) +sqrt(1-rho).*USprim e.*exp(1i*phi) );
        
        y = abs(r.^2);
        x = gamrnd(a,1,L,n,z)/a;
        irradiance = x.*y;
        f_mean=malaga_moment(1,a,b,gam,omegaprime);
        f_std=malaga_moment(2,a,b,gam,omegaprime);

        chan = (irradiance /f_mean );
        
end




end