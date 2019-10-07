function logL =Loglikelihood_fullsigome_nt(N,T1,k,m,Dy,Dy1,theta,theta_idx)
 

Phi1=theta(1:theta_idx(1),1);
Phi=[Phi1(1,1),Phi1(2,1) ;Phi1(3,1), Phi1(4,1) ];

omega11=theta(theta_idx(1)+1:theta_idx(2),1);
omega1=[omega11(1,1), omega11(3,1); omega11(3,1), omega11(2,1)];


Sigmae1=theta(theta_idx(2)+1:theta_idx(3),1);
Sigmae=[Sigmae1(1,1), Sigmae1(3,1); Sigmae1(3,1), Sigmae1(2,1)];

%Omega=2*eye(T1*(m+1))+diag(-1*ones((T1-1)*(m+1),1),m+1)+diag(-1*ones((T1-1)*(m+1),1),-(m+1)); % T(m+1) byT(m+1)
Omega=2*eye(T1+1)+diag(-1*ones((T1),1),1)+diag(-1*ones((T1),1),-1);

%L1=ones(T1+1); 
%L= tril(L1);
%e1=ones(T1+1,1);
psi=((Sigmae)^(-1))*omega1;
%invSigmaE=pinv((kron(pinv(L'*L),eye(m+1)) )+(kron((e1*e1'), (psi-eye(m+1)))))*kron(eye(T1+1), pinv(Sigmae));

SigmaE=kron(Omega,Sigmae);
BigOmega=kron(Omega,eye(m+1));
BigOmega(1:m+1,1:m+1)=psi;
SigmaE(1:m+1,1:m+1)=omega1;  % W
invSigmaE=(SigmaE)^(-1);
sqrSigmaE = invSigmaE^(0.5); %try inv



%%%%%% The problem may be happen here, I don't use Cholesky factorization to find incOvega^{-0.5},
%%%%%% because chol() function can not use in block tridiagonal matrix  %%%%%%%%%%%%%%%%%%%%%%%%   



phiDy=zeros(m+1,N,T1);
for tk=1:T1
    for ii=1:N
phiDy(:,ii,tk)= Phi*Dy1(:,ii,tk);  % m+1 by N by T1
    end
end

comphiDy=[];
for tkt=1:T1
comphiDy= cat(1,comphiDy,phiDy(:,:,tkt));   % T1*(m+1) by N   
end    

FirstDy=zeros(m+1,N); % first Dy m+1 by N


comphiDy=cat(1,FirstDy,comphiDy); % combine T0*(m+1) by N


invR_chi=Dy-comphiDy;   % T0*(m+1) by N   

C_N=(invR_chi*invR_chi')/N;  % T(m+1) by T(m+1)
D_N=(sqrSigmaE)*C_N*(sqrSigmaE);

[V,D]=eig(D_N); % Tm by k or P=pca(D_N) P=pca(Omega^(-0.5)*Q*A^(-0.5))

[Q,L1]=sortem(V,real(D));

lambda=diag(L1); 



logL=-((T1+1)/2)*(real(log(det(Sigmae))))-(1/2)*(real(log(det(BigOmega))))-(1/2)*sum(log(lambda(1:k*(m+1))))+(1/2)*(sum(lambda(1:k*(m+1)))-k*(m+1))-(1/2)*sum(lambda(1:(T1+1)*(m+1)));


logL = -logL;
% abs()imag() angle()

end