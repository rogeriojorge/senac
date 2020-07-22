tic;

inconds;

%NFP = 4;
%nphi = 500;
%R0 = 1;
% epsilonR = 0.04;
% epsilonZ = 0.04;
%etab = 1.51013;
%sigma0 = 0;
%iota0 = 3.1;

%Compute closed curve
phi=0:2*pi/(nphi-1):2*pi;
%R = R0+0.360407.*cos(1.*NFP.*phi)-0.0669714.*cos(2.*NFP.*phi)+0.00460894.*cos(3.*NFP.*phi)-0.000223451.*cos(4.*NFP.*phi);
R = R0+epsilonR.*cos(NFP.*phi);
x = R.*cos(phi);
y = R.*sin(phi);
%z = 0.4.*sin(1.*NFP.*phi)-0.0693047.*sin(2.*NFP.*phi)+0.00480132.*sin(3.*NFP.*phi)-0.000227014.*sin(4.*NFP.*phi);
z = epsilonZ.*sin(NFP.*phi);
%Derivatives of closed curve vector
r = [x; y; z];
dr = gradient(r,2*pi/(nphi-1));
sprime = sqrt(sum(dr.^2));
ddr = gradient(dr,2*pi/(nphi-1));
dddr = gradient(ddr,2*pi/(nphi-1));

%Curvature
curv = sqrt(sum(cross(dr,ddr).^2))./sqrt(sum(dr.^2)).^3;
curv(end-1)=curv(end-2);curv(end)=curv(end-1);curv(1)=curv(end);curv(2)=curv(1);
dcurv = gradient(curv,2*pi/(nphi-1));
%Torsion
tors = dot(dr,cross(ddr,dddr))./sum(cross(dr,ddr).^2);
tors(end)=tors(end-1);tors(1)=tors(end);tors(2)=tors(1);

curvfit=fit(phi',curv','fourier8');
torsfit=fit(phi',tors','fourier8');
sprimefit=fit(phi',sprime','fourier8');

curv=feval(curvfit,phi)';
tors=feval(torsfit,phi)';
sprime=feval(sprimefit,phi)';

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

%Compute derivative matrix
N = nphi; h = 2*pi/N; x = h*(1:N)';
column = [0 .5*(-1).^(1:N-1).*cot((1:N-1)*h/2)]';
D = toeplitz(column,column([1 N:-1:2]));

%Find sigma
sigma = 0.1*cos(NFP.*x);iota=iota0;
sprime = sprime';tors=tors';curv=curv';
tol=1.e-15;totalIt=50;totalLineSearch=15;deltaold=0;tic;
for i=1:totalIt
    [f,Jacf] = qs_residual_sigma(D,sigma,tors,curv,sprime,etab,nNormal,iota,L,N,sigma0);
    Jacf;
    deltay = -Jacf\f;
    if i>5 && norm((deltay-deltaold))/sqrt(N)<tol
        break
    end
    
    for j=1:totalLineSearch
        [newf,~] = qs_residual_sigma(D,sigma+deltay(1:N),tors,curv,sprime,etab,nNormal,iota+deltay(N+1),L,N,sigma0);
        if abs(norm(newf))<abs(norm(f))/1.2
            break;
        else
            deltay=deltay/2;
        end
    end
    
    sigma = sigma+deltay(1:N);
    iota  = iota+deltay(N+1);
    deltaold = deltay;
end
%plot(phi,sigma);
timeT=toc;
disp(['iota = ',num2str(iota),' after ',num2str(timeT),'s and ',num2str(i),' iterations. N = ',num2str(nNormal)]);

B0=1;
lambda=0;
pprime0=0;

Vpptemp=(1/32).*B0.^(-3).*etab.^(-2).*L.^(-2).*pi.^(-1).*sprime.*((-1).* ...
  B0.*curv.^(-2).*((-1).*L.^2.*lambda.^2+64.*((-1).*iota+nNormal) ...
  .^2.*pi.^2).*(etab.^4+curv.^4.*(1+sigma.^2))+L.^2.*((-16).*B0.* ...
  etab.^4+(-1).*B0.*curv.^(-2).*etab.^4.*lambda.^2+32.*etab.^2.*pi.* ...
  pprime0+B0.*curv.^2.*lambda.^2.*((-1)+(-1).*sigma.^2+16.*etab.^4.* ...
  (etab.^4+curv.^4.*(1+sigma.^2)).^(-1)))+64.*B0.*curv.*dcurv'.*L.*(( ...
  -1).*iota+nNormal).*pi.*sigma.*(1+curv.^(-4).*etab.^4+sigma.^2).* ...
  sprime.^(-1)+16.*B0.*curv.^2.*L.^2.*(1+curv.^(-4).*etab.^4+ ...
  sigma.^2).*tors.^2+16.*B0.*(1+curv.^(-4).*etab.^4+sigma.^2).*(4.* ...
  curv.^2.*((-1).*iota+nNormal).^2.*pi.^2.*(1+curv.^(-4).*etab.^4+ ...
  sigma.^2)+dcurv'.^2.*L.^2.*sprime.^(-2)+4.*etab.^2.*L.*(iota+(-1).* ...
  nNormal).*pi.*tors));

Vpp = 4*pi*pi*trapz(phi,Vpptemp)

%return;
%find mu and delta
mu=0.7+0.001.*cos(NFP.*x);
if nNormal~=0
    delta=0.05.*sin(NFP.*x)+0.1.*sin(2.*NFP.*x);
else
    delta=0.1.*sin(NFP.*x)+0.01.*sin(2.*NFP.*x)-NFP.*x/2;
end
tol=1.e-15;totalIt=70;totalLineSearch=12;deltaold=0;tic;
for i=1:totalIt
    [f,Jacf] = qs_residual_mercierGB(sigma,mu,delta,etab,curv);
    deltay = -Jacf\f;

    if i>5 && norm((deltay-deltaold))/sqrt(N)<tol
        break
    end
    
    for j=1:totalLineSearch
        [newf,~] = qs_residual_mercierGB(sigma,mu+deltay(1:N),delta+deltay(N+1:2*N),etab,curv);
        if abs(norm(newf))<abs(norm(f))/1.3
            break;
        else
            deltay=deltay/2;
        end
    end
    
    mu    = mu + deltay(1:N);
    delta = delta + deltay(N+1:2*N);
end
timeT=toc;
disp(['Found mu and delta after ',num2str(timeT),'s and ',num2str(i),' iterations.',]);
figure();hold on;plot(x,mu);

deltaold=delta;
if nNormal~=0
    plot(x,delta);
else
    plot(x,delta+NFP.*x/2);
    delta = delta+NFP.*x/2;
end

mufit=fit(phi',mu,'fourier8');
deltafit=fit(phi',delta,'fourier8');

muvalues = coeffvalues(mufit);
deltavalues = coeffvalues(deltafit);

%return;

%Write mu and delta to SENAC
fid = fopen('surf_input_QS.txt','w');
fprintf(fid,'NFP = %d\n', NFP);
fprintf(fid,'RAXIS = %5f %5f\n', R0, epsilonR);
fprintf(fid, 'ZAXIS = %5f', 0.0);
fprintf(fid, ' %5f\n', epsilonZ);
fprintf(fid, 'mu = %5f', muvalues(1));
for i=1:8
    fprintf(fid, ' %5f', muvalues(2*i));
end
if abs(deltaold(end))>abs(deltaold(1)+2*pi)
    fprintf(fid, '\ndelta = %5f', -NFP);
else
    fprintf(fid, '\ndelta = %5f', 0);
end
for i=1:8
    fprintf(fid, ' %5f', deltavalues(2*i+1));
end
fprintf(fid, '\niota = %5f', iota);
fprintf(fid, '\npsi31 = 0');
fprintf(fid, '\npsi32 = 0');
fprintf(fid, '\npsi33 = 0');
fprintf(fid, '\npsi34 = 0');
fclose(fid);

toc
