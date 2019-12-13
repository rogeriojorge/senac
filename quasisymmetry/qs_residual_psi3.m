function [f,Jacf] = qs_residual_psi3(D,sigma,tors,curv,kpok,sprime,etab,nNormal,iota,L,curvpp,curvppp,torsp,torspp,psi31,psi32,psi33,psi34)
%Solve the psi3 differential equations with psi3 quasisymmetric
    %compute A and B matrices of psi3'=A.*psi3+B using mathematica
    
    B0 = 1;c0=1;
    F0 = 2.*pi.*(1./L).*(iota-nNormal);
    kpok = kpok./sprime;curvpp=curvpp./sprime./sprime;curvppp=curvppp./sprime./sprime./sprime;torsp=torsp./sprime;torspp=torspp./sprime./sprime;
    
%Input matrices here
    
    fApsi3131=2.*kpok - 3.*F0.*sigma;
    fApsi3132=(5.*(etab.^2).*F0)./(2.*(curv.^2)) - ((curv.^2).*F0.*(1 + (sigma.^2)))./(2.*(etab.^2)) + 3.*tors;
    fApsi3133=3.*kpok;
    fApsi3134=(3/2).*(((etab.^2).*F0)./(curv.^2) + ((curv.^2).*F0.*(1 + (sigma.^2)))./(etab.^2) + 2.*tors);
    
    fApsi3231=(F0.*(3.*(etab.^4) + (curv.^4).*(1 + (sigma.^2))))./(2.*(curv.^2).*(etab.^2)) + tors;
    fApsi3232=-2.*kpok + F0.*sigma;
    fApsi3233=(-3.*F0.*((etab.^4) + (curv.^4).*(1 + (sigma.^2))))./(2.*(curv.^2).*(etab.^2)) - 3.*tors;
    fApsi3234=3.*kpok;
    
    fApsi3331=kpok - F0.*sigma;
    fApsi3332=-(F0.*(3.*(etab.^4) + (curv.^4).*(1 + (sigma.^2))))./(2.*(curv.^2).*(etab.^2)) - tors;
    fApsi3333=0.*kpok;
    fApsi3334=3.*(((curv.^2).*F0.*(1 + (etab.^4)./(curv.^4) + (sigma.^2)))./(2.*(etab.^2)) + tors);
    
    fApsi3431=(F0.*(3.*(etab.^4) + (curv.^4).*(1 + (sigma.^2))))./(2.*(curv.^2).*(etab.^2)) + tors;
    fApsi3432=kpok - F0.*sigma;
    fApsi3433=(-3).*(((curv.^2).*F0.*(1 + (etab.^4)./(curv.^4) + (sigma.^2)))./(2.*(etab.^2)) + tors);
    fApsi3434=0.*kpok;
        
    fBpsi31=(B0.*c0.*(2.*(etab.^10).*(F0.^3).*sigma + (curv.^3).*(etab.^6).*(curvppp - 3.*curvpp.*kpok - 2.*curvpp.*F0.*sigma) - (curv.^7).*(etab.^2).*(curvppp - 3.*curvpp.*kpok - 2.*curvpp.*F0.*sigma).*(1 + (sigma.^2)) + ...
            8.*(curv.^5).*curvpp.*(etab.^4).*sigma.*tors + (curv.^2).*(etab.^8).*F0.*(6.*kpok.*tors + 14.*F0.*sigma.*tors - torsp) + ...
            2.*(curv.^4).*(etab.^6).*((F0.^2).*kpok + 2.*(kpok.^3) + (F0.^3).*sigma.*(3 + 4.*(sigma.^2)) - kpok.*(tors.^2) + F0.*sigma.*(4.*(kpok.^2) + 9.*(tors.^2)) - 3.*tors.*torsp) + ...
            (curv.^8).*(etab.^2).*(-4.*(kpok.^3).*(1 + (sigma.^2)) + 2.*(F0.^3).*sigma.*(2 + 5.*(sigma.^2) + 3.*(sigma.^4)) - 2.*(F0.^2).*kpok.*(1 + 5.*(sigma.^2) + 4.*(sigma.^4)) + 3.*(etab.^2).*sigma.*tors + ...
            2.*kpok.*(1 + (sigma.^2)).*(tors.^2) + 2.*F0.*(sigma + (sigma.^3)).*(2.*(kpok.^2) - 3.*(tors.^2)) + 6.*(1 + (sigma.^2)).*tors.*torsp) - ...
            (curv.^10).*(1 + (sigma.^2)).*(2.*(etab.^2).*(kpok - 2.*F0.*sigma) - F0.*(1 + (sigma.^2)).*(2.*kpok.*tors - 2.*F0.*sigma.*tors + torsp)) + ...
            (curv.^6).*(etab.^4).*((etab.^2).*(kpok + 4.*F0.*sigma) + 4.*(F0.^2).*sigma.*(4 + 5.*(sigma.^2)).*tors - 4.*F0.*(2.*kpok.*tors + 4.*kpok.*(sigma.^2).*tors + (sigma.^2).*torsp) + ...
            2.*sigma.*(-2.*(kpok.^2).*tors - 2.*(tors.^3) + 2.*kpok.*torsp + torspp))))./(2.*(curv.^7).*(etab.^4));
    
    fBpsi32=-(B0.*c0.*(2.*(curv.^3).*curvpp.*(etab.^8).*F0 - 4.*(etab.^12).*(F0.^3) + 2.*(curv.^11).*curvpp.*F0.*(1 + (sigma.^2)).^2 - 4.*(curv.^7).*(etab.^4).*(-(curvppp.*sigma) + curvpp.*(F0 + 3.*kpok.*sigma + F0.*(sigma.^2))) - ...
            8.*(curv.^5).*curvpp.*(etab.^6).*tors - 16.*(curv.^2).*(etab.^10).*(F0.^2).*tors + 8.*(curv.^9).*curvpp.*(etab.^2).*(1 + (sigma.^2)).*tors + ...
            2.*(curv.^12).*F0.*((1 + (sigma.^2)).^2).*(-2.*F0.*kpok.*sigma + (F0.^2).*(1 + 2.*(sigma.^2)) - (tors.^2)) - 2.*(curv.^4).*(etab.^8).*F0.*(6.*(kpok.^2) - 2.*F0.*kpok.*sigma + (F0.^2).*(1 + 6.*(sigma.^2)) + 7.*(tors.^2)) + ...
            2.*(curv.^8).*(etab.^4).*(8.*(kpok.^3).*sigma - 2.*F0.*(kpok.^2).*(-3 + (sigma.^2)) + 4.*(F0.^2).*kpok.*(sigma + 2.*(sigma.^3)) - 2.*(F0.^3).*(-1 + (sigma.^4)) + (etab.^2).*tors - 4.*kpok.*sigma.*(tors.^2) + ...
            8.*F0.*(1 + 2.*(sigma.^2)).*(tors.^2) - 12.*sigma.*tors.*torsp) + (curv.^6).*(etab.^6).*...
            (3.*(etab.^2).*F0 + 4.*(kpok.^2).*tors + 32.*F0.*kpok.*sigma.*tors + (F0.^2).*(4 - 8.*(sigma.^2)).*tors + 4.*(tors.^3) - 4.*kpok.*torsp - 2.*torspp) + ...
            (curv.^10).*((etab.^4).*(-2.*kpok.*sigma + F0.*(2 + 3.*(sigma.^2))) + 2.*(etab.^2).*(1 + (sigma.^2)).*...
            (-2.*(kpok.^2).*tors - 2.*(tors.^3) + 6.*(F0.^2).*(tors + 2.*(sigma.^2).*tors) + 2.*kpok.*torsp - 4.*F0.*sigma.*(2.*kpok.*tors + torsp) + torspp))))./(4.*(curv.^9).*(etab.^4));
    
    fBpsi33=(B0.*c0.*(-((curv.^3).*(etab.^6).*(curvppp - 3.*curvpp.*kpok)) + 2.*(etab.^10).*(F0.^3).*sigma - (curv.^7).*(etab.^2).*(curvppp - 3.*curvpp.*kpok).*(1 + (sigma.^2)) + ...
           (curv.^6).*(-((etab.^6).*(kpok - 2.*F0.*sigma)) - 2.*(etab.^4).*F0.*(1 + (sigma.^2)).*(2.*kpok.*tors - torsp)) + (curv.^2).*(etab.^8).*F0.*(-6.*kpok.*tors + 2.*F0.*sigma.*tors + torsp) + ...
           2.*(curv.^4).*(etab.^6).*(-2.*(kpok.^3) - (F0.^2).*(kpok + 2.*kpok.*(sigma.^2)) + 2.*(F0.^3).*(sigma + (sigma.^3)) + kpok.*(tors.^2) + 2.*F0.*sigma.*((kpok.^2) - (tors.^2)) + 3.*tors.*torsp) + ...
           (curv.^8).*(etab.^2).*(-4.*(kpok.^3).*(1 + (sigma.^2)) + 2.*(F0.^3).*sigma.*((1 + (sigma.^2)).^2) - 2.*(F0.^2).*kpok.*(1 + 3.*(sigma.^2) + 2.*(sigma.^4)) + (etab.^2).*sigma.*tors + 2.*kpok.*(1 + (sigma.^2)).*(tors.^2) + ...
           4.*F0.*sigma.*(1 + (sigma.^2)).*((kpok.^2) - (tors.^2)) + 6.*(1 + (sigma.^2)).*tors.*torsp) + (curv.^10).*F0.*(1 + (sigma.^2)).*((etab.^2).*sigma + (1 + (sigma.^2)).*(2.*kpok.*tors - 2.*F0.*sigma.*tors + torsp))))./...
           (2.*(curv.^7).*(etab.^4));
    
    fBpsi34=-(B0.*c0.*(-2.*(curv.^3).*curvpp.*(etab.^8).*F0 + 4.*(etab.^12).*(F0.^3) + 2.*(curv.^11).*curvpp.*F0.*((1 + (sigma.^2)).^2) + 8.*(curv.^5).*curvpp.*(etab.^6).*tors + 16.*(curv.^2).*(etab.^10).*(F0.^2).*tors + ...
             8.*(curv.^9).*curvpp.*(etab.^2).*(1 + (sigma.^2)).*tors + 2.*(curv.^12).*F0.*((1 + (sigma.^2)).^2).*(-2.*F0.*kpok.*sigma + (F0.^2).*(1 + 2.*(sigma.^2)) - (tors.^2)) + ...
             2.*(curv.^4).*(etab.^8).*F0.*(6.*(kpok.^2) - 2.*F0.*kpok.*sigma + (F0.^2).*(5 + 6.*(sigma.^2)) + 7.*(tors.^2)) + ...
             2.*(curv.^8).*(etab.^4).*(-4.*(F0.^2).*kpok.*sigma.*(1 + (sigma.^2)) + 2.*(F0.^3).*(2 + 5.*(sigma.^2) + 3.*(sigma.^4)) + (etab.^2).*tors + 6.*F0.*(1 + (sigma.^2)).*((kpok.^2) + (tors.^2))) + ...
             (curv.^6).*(etab.^6).*(3.*(etab.^2).*F0 - 4.*(kpok.^2).*tors + 4.*(F0.^2).*(7 + 8.*(sigma.^2)).*tors - 4.*(tors.^3) + 4.*kpok.*torsp - 4.*F0.*sigma.*(2.*kpok.*tors + torsp) + 2.*torspp) + ...
             (curv.^10).*((etab.^4).*(2.*kpok.*sigma + F0.*(2 + (sigma.^2))) + 2.*(etab.^2).*(1 + (sigma.^2)).*...
             (-2.*(kpok.^2).*tors + 2.*(F0.^2).*(3 + 4.*(sigma.^2)).*tors - 2.*(tors.^3) + 2.*kpok.*torsp - 2.*F0.*sigma.*(2.*kpok.*tors + torsp) + torspp))))./(4.*(curv.^9).*(etab.^4));
            
    fApsi3131=sprime.*fApsi3131;fApsi3132=sprime.*fApsi3132;fApsi3133=sprime.*fApsi3133;fApsi3134=sprime.*fApsi3134;
    fApsi3231=sprime.*fApsi3231;fApsi3232=sprime.*fApsi3232;fApsi3233=sprime.*fApsi3233;fApsi3234=sprime.*fApsi3234;
    fApsi3331=sprime.*fApsi3331;fApsi3332=sprime.*fApsi3332;fApsi3333=sprime.*fApsi3333;fApsi3334=sprime.*fApsi3334;
    fApsi3431=sprime.*fApsi3431;fApsi3432=sprime.*fApsi3432;fApsi3433=sprime.*fApsi3433;fApsi3434=sprime.*fApsi3434;
    fBpsi31=sprime.*fBpsi31;fBpsi32=sprime.*fBpsi32;fBpsi33=sprime.*fBpsi33;fBpsi34=sprime.*fBpsi34;
    fApsi31=psi31.*fApsi3131+psi32.*fApsi3132+psi33.*fApsi3133+psi34.*fApsi3134;
    fApsi32=psi31.*fApsi3231+psi32.*fApsi3232+psi33.*fApsi3233+psi34.*fApsi3234;
    fApsi33=psi31.*fApsi3331+psi32.*fApsi3332+psi33.*fApsi3333+psi34.*fApsi3334;
    fApsi34=psi31.*fApsi3431+psi32.*fApsi3432+psi33.*fApsi3433+psi34.*fApsi3434;
    
    f = [D*psi31-fApsi31-fBpsi31;
         D*psi32-fApsi32-fBpsi32;
         D*psi33-fApsi33-fBpsi33;
         D*psi34-fApsi34-fBpsi34;
        ];
    
    Jacf = [D-diag(fApsi3131) -diag(fApsi3132) -diag(fApsi3133) -diag(fApsi3134);
            -diag(fApsi3231) D-diag(fApsi3232) -diag(fApsi3233) -diag(fApsi3234);
            -diag(fApsi3331) -diag(fApsi3332) D-diag(fApsi3333) -diag(fApsi3334);
            -diag(fApsi3431) -diag(fApsi3432) -diag(fApsi3433) D-diag(fApsi3434)];
     
end

