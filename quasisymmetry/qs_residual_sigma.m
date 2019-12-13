function [f,Jacf] = qs_residual_sigma(D,sigma,tors,curv,sprime,etab,nNormal,iota,L,N,sigma0)
    
    F0 = 2*pi*(1/L)*(iota-nNormal);
    f = [D*sigma+sprime.*(F0.*(1+sigma.^2+etab^4./(curv.^4))+2*etab^2.*tors./(curv.^2));
        sigma(end)-sigma0];
     
    dfdsigma = D - diag(sprime.*F0*2.*sigma);
    dfdiota  = sprime.*2*pi*(1/L).*(1+sigma.^2+etab^4./(curv.^4));

    Jacf = [dfdsigma dfdiota;
            zeros(1,N-1) 1 0];

end

