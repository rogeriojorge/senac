function [psi31,psi32,psi33,psi34,psi32constraint,sigma,iota,curv,x,phi,nNormal] = QStriang(NFP,nphi,R0,epsilonR1,epsilonR2,epsilonZ1,epsilonZ2,etab,sigma0,iota0)

%Compute closed curve
phi=0:2*pi/(nphi-1):2*pi;
R = R0+epsilonR1.*cos(NFP.*phi)+epsilonR2.*cos(2.*NFP.*phi);
dR   = -epsilonR1.*NFP.*sin(NFP.*phi)-epsilonR2.*2.*NFP.*sin(2.*NFP.*phi);
ddR  = -epsilonR1.*NFP.*NFP.*cos(NFP.*phi)-epsilonR2.*4.*NFP.*NFP.*cos(2.*NFP.*phi);
dddR = epsilonR1.*NFP.*NFP.*NFP.*sin(NFP.*phi)+epsilonR2.*8.*NFP.*NFP.*NFP.*sin(2.*NFP.*phi);
z    = -epsilonZ1.*sin(NFP.*phi)-epsilonZ2.*sin(2.*NFP.*phi);
dz   = -epsilonZ1.*NFP.*cos(NFP.*phi)-epsilonZ2.*2.*NFP.*cos(2.*NFP.*phi);
ddz  = epsilonZ1.*NFP.*NFP.*sin(NFP.*phi)+epsilonZ2.*4.*NFP.*NFP.*sin(2*NFP.*phi);
dddz = epsilonZ1.*NFP.*NFP.*NFP.*cos(NFP.*phi)+epsilonZ2.*8.*NFP.*NFP.*NFP.*cos(2*NFP.*phi);
x = R.*cos(phi);
y = R.*sin(phi);

%Compute derivative matrix
N = nphi; h = 2*pi/N;
column = [0 .5*(-1).^(1:N-1).*cot((1:N-1)*h/2)]';
D = toeplitz(column,column([1 N:-1:2]));

%Derivatives of closed curve vector
r = [x; y; z];

%dr = gradient(r,2*pi/(nphi-1));
dr  = [-R.*sin(phi)+dR.*cos(phi);
       R.*cos(phi)+dR.*sin(phi); 
       dz];

%ddr = gradient(dr,2*pi/(nphi-1));
ddr = [-R.*cos(phi)-2*dR.*sin(phi)+ddR.*cos(phi);
       -R.*sin(phi)+2*dR.*cos(phi)+ddR.*sin(phi);
       ddz];

%dddr = gradient(ddr,2*pi/(nphi-1));
dddr = [R.*sin(phi)-3*dR.*cos(phi)-3*ddR.*sin(phi)+dddR.*cos(phi);
        -R.*cos(phi)-3*dR.*sin(phi)+3*ddR.*cos(phi)+dddR.*sin(phi);
        dddz];
   
sprime = sqrt(sum(dr.^2));
%sprimefit=fit(phi',sprime','fourier8');
%sprime=feval(sprimefit,phi)';

%Curvature
curv = sqrt(sum(cross(dr,ddr).^2))./sqrt(sum(dr.^2)).^3;
%curv(end-1)=curv(end-2);curv(end)=curv(end-1);curv(1)=curv(end);curv(2)=curv(1);
%curvfit=fit(phi',curv','fourier8');
%curv=feval(curvfit,phi)';
%dcurv1 = gradient(curv,2*pi/(nphi-1));
dcurv=D*curv';dcurv=dcurv';
%dcurv(end)=dcurv(end-1);dcurv(1)=dcurv(end);dcurv(2)=dcurv(1);
%dcurvfit=fit(phi',dcurv','fourier8');
%dcurv=feval(dcurvfit,phi)';
kpok=dcurv./curv;

%Torsion
tors = dot(dr,cross(ddr,dddr))./sum(cross(dr,ddr).^2);
%tors(end-2)=tors(end-3);tors(end-1)=tors(end-2);tors(end)=tors(end-1);tors(1)=tors(end);tors(2)=tors(1);tors(3)=tors(2);tors(4)=tors(3);
%torsfit=fit(phi',tors','fourier8');
%tors=feval(torsfit,phi)';

%Length of the axis
L = trapz(phi,sprime);

%Number of rotations of the normal vector
tangent = dr./sprime;
normal = gradient(tangent,2*pi/(nphi-1))./sprime./curv;
quadrant = 0; nNormal = 0; diffQuadrant = 0; qnew=0;
for i=1:nphi
    normalx=normal(1,i)*cos(2*pi*i/nphi)+normal(2,i)*sin(2*pi*i/nphi);
    normaly=normal(3,i);
    if norm(normalx) > norm(normaly)
        if normalx>0
            qnew=2;
        else
            qnew=4;
        end
    else
        if normaly > 0
            qnew=1;
        else
            qnew=3;
        end
    end
    if quadrant == 0
        quadrant = qnew;
    else
        diffQuadrant = qnew - quadrant;
    end
    if abs(diffQuadrant) == 3
        diffQuadrant = diffQuadrant/3;
    end
    nNormal = nNormal + diffQuadrant/2;
    quadrant = qnew;
end

%Find sigma
x = h*(1:N)';
sigma = 0.1.*cos(NFP.*x);iota=iota0;
sprime = sprime';tors=tors';curv=curv';
tol=1.e-4;totalIt=30;totalLineSearch=8;deltaold=0;
for i=1:totalIt
    [f,Jacf] = qs_residual_sigma(D,sigma,tors,curv,sprime,etab,nNormal,iota,L,N,sigma0);
    deltay = -Jacf\f;
    
    if i>5 && norm((deltay-deltaold))/sqrt(N)<tol
        break
    end
    
    for j=1:totalLineSearch
        [newf,~] = qs_residual_sigma(D,sigma+deltay(1:N),tors,curv,sprime,etab,nNormal,iota+deltay(N+1),L,N,sigma0);
        if abs(norm(newf))<abs(norm(f))
            break;
        else
            deltay=deltay/2;
        end
    end
    
    sigma = sigma+deltay(1:N);
    iota  = iota+deltay(N+1);
    deltaold = deltay;
end

%Derivatives of torsion/curvature
%dcurvfit=fit(phi',dcurv','fourier8');
%plot(dcurv)
%dcurv=feval(dcurvfit,phi)';
%hold on;plot(dcurv);

curvpp = gradient(dcurv,2*pi/(nphi-1));
%curvpp(end-1)=curvpp(end-2);curvpp(end)=curvpp(end-1);curvpp(1)=curvpp(end);curvpp(2)=curvpp(1);curvpp(3)=curvpp(2);
%curvppfit=fit(phi',curvpp','fourier8');
%curvpp=feval(curvppfit,phi)';

curvppp = gradient(curvpp,2*pi/(nphi-1));
%curvppp(end-1)=curvppp(end-2);curvppp(end)=curvppp(end-1);curvppp(1)=curvppp(end);curvppp(2)=curvppp(1);curvppp(3)=curvppp(2);
%curvpppfit=fit(phi',curvppp','fourier8');
%curvppp=feval(curvpppfit,phi)';

%torsp = gradient(tors',2*pi/(nphi-1));
torsp = D*tors;
%torspfit=fit(phi',torsp','fourier8');
%torsp=feval(torspfit,phi)';

torspp = gradient(torsp,2*pi/(nphi-1));
%torsppfit=fit(phi',torspp','fourier8');
%torspp=feval(torsppfit,phi)';

%Find psi3
psi31 = 0.2+0.*cos(NFP.*x);
psi32 = 0.01*sin(NFP.*x);
psi33 = 0.1+0.*cos(NFP.*x);
psi34 = 0.01*sin(NFP.*x);
kpok=kpok';curvpp=curvpp';curvppp=curvppp';%torsp=torsp';torspp=torspp';
tol=1.e-6;totalIt=30;totalLineSearch=6;deltaold=0;resd=0;
for i=1:totalIt
    [f,Jacf] = qs_residual_psi3(D,sigma,tors,curv,kpok,sprime,etab,nNormal,iota,L,curvpp,curvppp,torsp,torspp,psi31,psi32,psi33,psi34);
    deltay = -Jacf\f;
    resid=norm((deltay-deltaold))/sqrt(N);
    if i>5 && resid<tol
        break
    end
    
    for j=1:totalLineSearch
        [newf,~] = qs_residual_psi3(D,sigma,tors,curv,kpok,sprime,etab,nNormal,iota,L,curvpp,curvppp,torsp,torspp,psi31+deltay(1:N),psi32+deltay(N+1:2*N),psi33+deltay(2*N+1:3*N),psi34+deltay(3*N+1:4*N));
        if abs(norm(newf))<abs(norm(f))
            break;
        else
            deltay=deltay/2;
        end
    end
    
    psi31 = psi31+deltay(1:N);
    psi32 = psi32+deltay(N+1:2*N);
    psi33 = psi33+deltay(2*N+1:3*N);
    psi34 = psi34+deltay(3*N+1:4*N);
    deltaold = deltay;
end
%disp([num2str(i),' iterations with max psi31 = ',num2str(max(max(psi31)))]);

B0 = 1;c0 = 1;F0 = 2.*pi.*(1./L).*(iota-nNormal);
psi32constraint =(B0.*c0.*(4.*(curv.^3).*curvpp.*(etab.^2).*(-2.*kpok + F0.*sigma) + 4.*(etab.^6).*(F0.^2).*(kpok + F0.*sigma) + 2.*(curv.^2).*(etab.^4).*F0.*(2.*kpok.*tors + 6.*F0.*sigma.*tors - 3.*torsp) + ...
                 (curv.^6).*((etab.^2).*(-2.*kpok + F0.*sigma) - 2.*F0.*(1 + (sigma.^2)).*(2.*kpok.*tors - 2.*F0.*sigma.*tors + torsp)) + ...
                 4.*(curv.^4).*(etab.^2).*(-((F0.^2).*kpok.*(1 + (sigma.^2))) + (F0.^3).*(sigma + (sigma.^3)) - F0.*sigma.*((kpok.^2) - 2.*(tors.^2)) + 2.*((kpok.^3) - tors.*torsp))))./(4.*(curv.^5).*(etab.^2).*F0);

end