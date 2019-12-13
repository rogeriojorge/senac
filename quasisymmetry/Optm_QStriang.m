tic;

inconds;

etabar=linspace(etabmin,etabmax,netab);
eR1=epsilonR1;
eR2=epsilonR2;
eZ1=epsilonZ1;
eZ2=epsilonZ2;
%QSconst=zeros(netab,nR,nZ);
QSconst=zeros(netab,1);

%h = waitbar(0,'Starting...');
%s = clock;
parfor i=1:netab
    eB=etabar(i);
    %parfor j=1:nR
        %eR1=epsilonR1;
        %for k=1:nZ
            %eZ=epsilonZ(k);
            [~,psi32n,~,~,psi32constraint] = QStriang(NFP,nphi,R0,eR1,eR2,eZ1,eZ2,eB,sigma0,iota0);
            QSconst(i)=norm(psi32n-psi32constraint)^2;
        %end
    %end
     %if i ==1
     % is = etime(clock,s);
     % esttime = is * netab;
     %end
     %h = waitbar(i/netab,h,['remaining time = ',secs2hms(esttime-etime(clock,s))]);
end
%close(h);

%[a,b]=min(QSconst,[],'all','linear');
%[a1,a2,a3]=ind2sub(size(QSconst),b);
[a,b]=min(QSconst,[],'all','linear');
[a1]=ind2sub(size(QSconst),b);

etab=etabar(a1);
%epsilonR=epsilonR(a2);
%epsilonZ=epsilonZ(a3);
nphi = 650;

[psi31n,psi32n,psi33n,psi34n,psi32constraint,sigma,iota,curv,x,phi,nNormal] = QStriang(NFP,nphi,R0,eR1,eR2,eZ1,eZ2,etab,sigma0,iota0);
figure();hold on;plot(psi32n);plot(psi32constraint);legend('\psi_{32}','constraint on \psi_{32}');
%timeT=toc;disp(['Found iota = ',num2str(iota),', etabar=',num2str(etab),', epsilonR = ',num2str(epsilonR),', epsilonZ = ',num2str(epsilonZ),' after ',secs2hms(timeT),'s']);
timeT=toc;disp(['Found iota = ',num2str(iota),', etabar=',num2str(etab),' after ',secs2hms(timeT),'s']);

%find mu and delta
N=length(x);
mu=0.8+0.1.*cos(NFP.*x);

if nNormal ~= 0
    delta=0.01.*sin(NFP.*x);
else
    delta=0.01.*sin(NFP.*x)-plusminus.*NFP.*x/2;
end

tol=1.e-4;totalIt=30;totalLineSearch=8;deltaold=0;tic;
for i=1:totalIt
    [f,Jacf] = qs_residual_mercierGB(sigma,mu,delta,etab,curv);
    deltay = -Jacf\f;

    if i>5 && norm((deltay-deltaold))/sqrt(N)<tol
        break
    end
    
    for j=1:totalLineSearch
        [newf,~] = qs_residual_mercierGB(sigma,mu+deltay(1:N),delta+deltay(N+1:2*N),etab,curv);
        if abs(norm(newf))<abs(norm(f)) & min(mu)+deltay(1:N)>-1 & max(mu)+deltay(1:N)<1
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
figure();hold on;plot(x,mu);plot(x,delta+plusminus.*NFP.*x/2);

deltaold=delta;
if abs(delta(end))>0.9*abs(delta(1)+2*pi)
    delta = delta+plusminus.*NFP.*x/2;
end

mufit=fit(phi',mu,'fourier8');
deltafit=fit(phi',delta,'fourier8');

muvalues = coeffvalues(mufit);
deltavalues = coeffvalues(deltafit);

%Write mu and delta to SENAC
fid = fopen('surf_input_QS.txt','w');
fprintf(fid,'NFP = %d\n', NFP);
fprintf(fid,'RAXIS = %5f %5f %5f\n', R0, epsilonR1, epsilonR2);
fprintf(fid, 'ZAXIS = %5f %5f %5f\n', 0.0, -epsilonZ1, -epsilonZ2);
fprintf(fid, 'mu = %5f', muvalues(1));
for i=1:8
    fprintf(fid, ' %5f', muvalues(2*i));
end
if abs(deltaold(end))>0.9*abs(deltaold(1)+2*pi)
    fprintf(fid, '\ndelta = %5f', -plusminus*NFP);
else
    fprintf(fid, '\ndelta = %5f', 0);
end
for i=1:8
    fprintf(fid, ' %5f', deltavalues(2*i+1));
end
fprintf(fid, '\niota = %5f', iota);
fprintf(fid, '\netabar = %5f', etab);

%Write psi3 to SENAC
psi31values=coeffvalues(fit(phi',psi31n,'fourier8'));
psi32values=coeffvalues(fit(phi',psi32n,'fourier8'));
psi33values=coeffvalues(fit(phi',psi33n,'fourier8'));
psi34values=coeffvalues(fit(phi',psi34n,'fourier8'));

fprintf(fid, '\npsi31 = ');
for i=1:17
    fprintf(fid, ' %5f', psi31values(i));
end
fprintf(fid, '\npsi32 = ');
for i=1:17
    fprintf(fid, ' %5f', psi32values(i));
end
fprintf(fid, '\npsi33 = ');
for i=1:17
    fprintf(fid, ' %5f', psi33values(i));
end
fprintf(fid, '\npsi34 = ');
for i=1:17
    fprintf(fid, ' %5f', psi34values(i));
end

fclose(fid);