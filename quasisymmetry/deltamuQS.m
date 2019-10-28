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

kpokfit=fit(phi',kpok','fourier8');
torsfit=fit(phi',tors','fourier8');
sprimefit=fit(phi',sprime','fourier8');

%Initial conditions for QS equations
solinit = bvpinit(phi, @sigmaGuess, iota0);

%Solution of QS equations
options = bvpset('RelTol',10^-4,'AbsTol',10^-5,'NMax',nphi);

sol = bvp4c(@sigmaEq, @sigmaBC, solinit, options, torsfit, kpokfit, sprimefit);
%sigma = sol.y(1,:);
iota=sol.parameters
plot(sol.x,sol.y(1,:),'--o',sol.x,sol.y(2,:),'--o')

toc

function yinit = sigmaGuess(s)
    inconds;
    yinit = [0.67-0.03*cos(NFP*s)
             -NFP*s/2+0.44*sin(NFP*s)];
end

function dYdt = sigmaEq(s,Y,iota,torsfit,kpokfit,sprimefit)
    kpokfit1=kpokfit(s)/sprimefit(s);
    torsfit1=torsfit(s)/sprimefit(s);
    dYdt = [sprimefit(s)*(1-Y(1)^2)*2*(kpokfit1*cos(2*Y(2))-(iota/sqrt(1-Y(1)^2)+torsfit1)*sin(2*Y(2)))
            sprimefit(s)*(-(kpokfit1*sin(2*Y(2))+(iota/sqrt(1-Y(1)^2)+torsfit1)*cos(2*Y(2)))/Y(1)+iota/sqrt(1-Y(1)^2)+torsfit1)];
end

function res = sigmaBC(ya,yb,~,~,~,~)
    inconds;
    res = [ya(1) - yb(1)
           ya(2) - yb(2)
           ya(2)];
end