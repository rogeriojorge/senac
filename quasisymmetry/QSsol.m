tic;
inconds;

%Compute closed curve
phi=0:2*pi/(nphi-1):2*pi;
R = R0+epsilonR.*cos(NFP.*phi);
x = R.*cos(phi);
y = R.*sin(phi);
z = -epsilonZ.*sin(NFP.*phi);
%Derivatives of closed curve vector
r = [x; y; z];
dr = gradient(r,2*pi/(nphi-1));
sprime = sqrt(sum(dr.^2));
ddr = gradient(dr,2*pi/(nphi-1));
dddr = gradient(ddr,2*pi/(nphi-1));

%Curvature
curv = sqrt(sum(cross(dr,ddr).^2))./sqrt(sum(dr.^2)).^3;
dcurv = gradient(curv,2*pi/(nphi-1));
curv(end-1)=curv(end-2);curv(end)=curv(end-1);curv(1)=curv(end);curv(2)=curv(1);
dcurv(end-2)=dcurv(end-3);dcurv(end-1)=dcurv(end-2);dcurv(end)=dcurv(end-1);dcurv(1)=dcurv(end);dcurv(2)=dcurv(1);dcurv(3)=dcurv(2);
kpok=dcurv./curv;

%Torsion
tors = dot(dr,cross(ddr,dddr))./sum(cross(dr,ddr).^2);
tors(end-2)=tors(end-3);tors(end-1)=tors(end-2);tors(end)=tors(end-1);tors(1)=tors(end);tors(2)=tors(1);tors(3)=tors(2);

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

curvfit=fit(phi',curv','fourier8');
torsfit=fit(phi',tors','fourier8');
sprimefit=fit(phi',sprime','fourier8');

%Initial conditions for QS equations
solinit = bvpinit(phi, @sigmaGuess, iota0);

%Solution of QS equations
options = bvpset('RelTol',10^-3,'AbsTol',10^-5,'NMax',nphi);

sol = bvp4c(@sigmaEq, @sigmaBC, solinit, options, torsfit, curvfit, sprimefit, L, nNormal, etab);
sigma = sol.y(1,:);
iota=sol.parameters;

% %Scan in etabar
% nPoints=100;
% etabVec=linspace(-5,5,nPoints);
% iotaVec=zeros(nPoints,1);
% parfor i=1:nPoints
%     sol = bvp4c(@sigmaEq, @sigmaBC, solinit, options, torsfit, curvfit, sprimefit, L, nNormal, etabVec(i));
%     iotaVec(i)=sol.parameters;
% end
% figure();plot(etabVec,iotaVec);

disp(['iota on axis = ',num2str(iota)]);

nphiSigma = size(sol.x);
nphiSigma = nphiSigma(2);
options = optimset('Display','off','FinDiffType','central','MaxFunEvals',1000,'MaxIter',3000);
mercierParams = zeros(nphi,2);mutemp=0.5;deltatemp=0;
for i=1:nphi
    isigma = ceil(i*nphiSigma/nphi);
    fun = @(x) mercierFromGB(x,sigma(isigma),curv(i),etab);
    x0 = [mutemp,deltatemp-NFP*2*pi/2/nphi];
    mercierParams(i,:) = fsolve(fun,x0,options);
    mutemp=mercierParams(i,1);
    deltatemp=mercierParams(i,2);
end

mu = mercierParams(:,1);
if abs(mercierParams(end,2))>abs(mercierParams(1,2)+2*pi)
    delta = mercierParams(:,2)+NFP.*phi'/2;
else
    delta = mercierParams(:,2);
end

mufit=fit(phi',mu,'fourier8');
deltafit=fit(phi',delta,'fourier8');

muvalues = coeffvalues(mufit);
deltavalues = coeffvalues(deltafit);

%Write mu and delta to SENAC
fid = fopen('surf_input_QS.txt','w');
fprintf(fid,'NFP = %d\n', NFP);
fprintf(fid,'RAXIS = %5f %5f\n', R0, epsilonR);
fprintf(fid, 'ZAXIS = %5f', 0.0);
fprintf(fid, ' %5f\n', -epsilonZ);
fprintf(fid, 'mu = %5f', muvalues(1));
for i=1:8
    fprintf(fid, ' %5f', muvalues(2*i));
end
if abs(mercierParams(end,2))>abs(mercierParams(1,2)+2*pi)
    fprintf(fid, '\ndelta = %5f', -NFP);
else
    fprintf(fid, '\ndelta = %5f', 0);
end
for i=1:8
    fprintf(fid, ' %5f', deltavalues(2*i+1));
end
fprintf(fid, '\niota = %5f', iota);
fclose(fid);

toc

function yinit = sigmaGuess(s)
    inconds;
    yinit = sin(NFP*s);
end

function dYdt = sigmaEq(s,Y,iota,torsfit,curvfit,sprimefit,L,nNormal,etab)

    dYdt = -sprimefit(s)*(2*pi*(1/L)*(iota-nNormal)*(1+Y(1)^2+etab^4/curvfit(s)^4)+2*etab^2*torsfit(s)/curvfit(s)^2);
    
end

function res = sigmaBC(ya,yb,~,~,~,~,~,~,~)
    inconds;
    res = [ya(1) - yb(1)
           ya(1)-sigma0];
end

function F = mercierFromGB(x, sigma, k, etab)
    F(1) = sigma-x(1)*sin(2*x(2))/sqrt(1-x(1)^2);
    F(2) = etab^2/k^2-(1-x(1)*cos(2*x(2)))/sqrt(1-x(1)^2);
end