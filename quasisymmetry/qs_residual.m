function [f,Jacf] = qs_residual(D,y1,y2,curv,kpok,tors,sprime,F0,N,NFP,x,etab)
    
    F0 = F0.*sprime;
    tors = tors.*sprime;
     
    f = [D*y1 - NFP/2 + (tors + F0.*cosh(y2)).*(-1 + cos(2.*y1-NFP.*x).*coth(y2))+kpok.*coth(y2).*sin(2.*y1-NFP.*x);
         D*y2 - 2.*kpok.*cos(2.*y1-NFP.*x) + 2.*(tors + F0.*cosh(y2)).*sin(2.*y1-NFP.*x);
         cosh(y2)-sinh(y2).*cos(2.*y1-NFP.*x)-(etab^2)./(curv.^2)];
         %y1(1)];
     
     df1dy1 = D + diag(2.*coth(y2).*(kpok.*cos(2.*y1-NFP.*x) - (tors + F0.*cosh(y2)).*sin(2.*y1-NFP.*x)));
     df2dy1 = diag(4.*(cos(2.*y1-NFP.*x).*(tors + F0.*cosh(y2)) + kpok.*sin(2.*y1-NFP.*x)));
     df3dy1 = diag(sinh(y2).*2.*sin(2.*y1-NFP.*x));
     df1dy2 = diag(F0.*cos(2.*y1-NFP.*x).*cosh(y2) - (F0 + (csch(y2).^3).*(cos(2.*y1-NFP.*x).*(tors + F0.*cosh(y2)) + kpok.*sin(2.*y1-NFP.*x))).*sinh(y2));
     df2dy2 = D + diag(2.*F0.*sin(2.*y1-NFP.*x).*sinh(y2));
     df3dy2 = diag(sinh(y2)-cosh(y2).*cos(2.*y1-NFP.*x));
     df1dF0 = cosh(y2).*(-1 + cos(2.*y1-NFP.*x).*coth(y2));
     df2dF0 = 2.*cosh(y2).*sin(2.*y1-NFP.*x);
     df3dF0 = zeros(N,1);

%     f = [D*y1-NFP/2+(((tors.*(1-y2.^2)+F0.*sqrt(1-y2.^2)).*(y2-cos(2.*y1-NFP.*x)))./(-1+y2.^2)+kpok.*sin(2.*y1-NFP.*x))./y2;
%          D*y2+2.*kpok.*(-1+y2.^2).*cos(2.*y1-NFP.*x)+2.*(tors.*(1-y2.^2)+F0.*sqrt(1-y2.^2)).*sin(2.*y1-NFP.*x);
%          y1(1)];
%      
%     df1dy1 = D+diag((2.*(kpok.*cos(2.*y1-NFP.*x)+((tors.*(1-y2.^2)+F0.*sqrt(1 - y2.^2)).*sin(2.*y1-NFP.*x))./(-1 + y2.^2)))./y2);
%     df2dy1 = diag(4.*(tors.*(1-y2.^2)+F0.*sqrt(1-y2.^2)).*cos(2.*y1-NFP.*x)-4.*kpok.*(-1 + y2.^2).*sin(2.*y1-NFP.*x));
%     df1dy2 = diag(-((F0.*y2.^3+(F0-2.*F0.*y2.^2+tors.*(1 - y2.^2).^(3./2)).*cos(2.*y1-NFP.*x)+kpok.*(1 - y2.^2).^(3./2).*sin(2.*y1-NFP.*x))./(y2.^2.*(1-y2.^2).^(3./2))));
%     df2dy2 = D+diag(4.*kpok.*y2.*cos(2.*y1-NFP.*x)+(-4.*tors.*y2-(2.*F0.*y2)./sqrt(1-y2.^2)).*sin(2.*y1-NFP.*x));
%     df1dF0 = (-y2+cos(2.*y1-NFP.*x))./(y2.*sqrt(1-y2.^2));
%     df2dF0 = 2.*sqrt(1-y2.^2).*sin(2.*y1-NFP.*x);
% 
    Jacf = [df1dy1 df1dy2 df1dF0;
            df2dy1 df2dy2 df2dF0;
            df3dy1 df3dy2 df3dF0];
            %1 zeros(1,2*N)];

end

