tic;
inconds;

%Compute closed curve
phi=0:2*pi/(nphi-1):2*pi;
R = R0+epsilonR.*cos(NFP.*phi);
x = R.*cos(phi);
y = R.*sin(phi);
z = -epsilonZ.*sin(NFP.*phi);


%Derivatives of closed curve vector
r = zeros(3,nphi);dr=r;ddr=r;dddr=r;
r = [x; y; z];
dr = gradient(r,2*pi/(nphi-1));
sprime = sqrt(sum(dr.^2));
ddr = gradient(dr,2*pi/(nphi-1));
dddr = gradient(ddr,2*pi/(nphi-1));

%Curvature
curv = sqrt(sum(cross(dr,ddr).^2))./sqrt(sum(dr.^2)).^3;
curv(end-1)=curv(end-2);curv(end)=curv(end-1);curv(1)=curv(end);curv(2)=curv(1);
dcurv = gradient(curv,2*pi/(nphi-1));
dcurv(end-2)=dcurv(end-3);dcurv(end-1)=dcurv(end-2);dcurv(end)=dcurv(end-1);dcurv(1)=dcurv(end);dcurv(2)=dcurv(1);dcurv(3)=dcurv(2);
kpok=dcurv./curv;

%Torsion
tors = dot(dr,cross(ddr,dddr))./sum(cross(dr,ddr).^2);
tors(end)=tors(end-1);tors(1)=tors(end);tors(2)=tors(1);

curvfit=fit(phi',curv','fourier8');
kpokfit=fit(phi',kpok','fourier8');
torsfit=fit(phi',tors','fourier8');
sprimefit=fit(phi',sprime','fourier8');

curv=feval(curvfit,phi)';
kpok=feval(kpokfit,phi)';
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

sprime = sprime';curv=curv';kpok=kpok';tors=tors';

y1 = 0.1.*sin(NFP.*x);
y2 = 0.5+0.01.*cos(NFP.*x);
F0=0.9;%kpok=0;tors=0;sprime=L/2/pi;
tol=1.e-5;totalIt=35;totalLineSearch=8;resid=zeros(1,totalIt);resid1=resid;deltaold=0;
for i=1:totalIt
    [f,Jacf] = qs_residual(D,y1,y2,curv,kpok,tors,sprime,F0,N,NFP,x,etab);
    deltay = -Jacf\f;
    deltay(1)=0;deltay(N)=0;
    
    for j=1:totalLineSearch
        [newf,~] = qs_residual(D,y1+deltay(1:N),y2+deltay(N+1:2*N),curv,kpok,tors,sprime,F0+deltay(2*N+1),N,NFP,x,etab);
        if abs(norm(newf))<abs(norm(f)) %&& min(y2+deltay(N+1:2*N))>-1 && max(y2+deltay(N+1:2*N))<1
            break;
        else
            deltay(2:N-1)=deltay(2:N-1)/2;
            deltay(N+1:2*N)=deltay(N+1:2*N)/2;
        end
    end

    resid(i)=abs((abs(deltay(2*N+1))-abs(deltaold)));
    if i>5 && resid(i)<tol
        break;
    end
    
    resid1(i)=deltay(2*N+1);
    y1 = y1 + deltay(1:N);
    y2 = y2 + deltay(N+1:2*N);
    F0 = F0 + deltay(2*N+1);
    deltaold = deltay(2*N+1);
end

figure();hold on;plot(x,y1);plot(x,tanh(y2));
figure();hold on;plot(x,curv.*sqrt(cosh(y2)-sinh(y2).*cos(2.*y1-NFP.*x)));
%figure();hold on;plot(resid);plot(resid1);

timeT=toc;
disp(['iota = ',num2str(F0*L/(2*pi)+nNormal),' after ',num2str(timeT),'s and ',num2str(i),' iterations.',]);