etab=linspace(0.7,1.6,15);
NFP=3;
epsilonR=0.0322;
epsilonZ=0.0278;
nphi=250;

maxpsi31=zeros(1,length(etab));
for i = 1:length(etab)
    [psi31,psi32,psi33,psi34,psi32constraint,sigma,iota,curv,x,phi,nNormal] = QStriang(NFP,nphi,R0,epsilonR,epsilonZ,etab(i),sigma0,iota0);
    maxpsi31(i)=max(max(psi31));
end

semilogy(etab,maxpsi31)