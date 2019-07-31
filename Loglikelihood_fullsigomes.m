function logL =Loglikelihood_fullsigomes(N,T1,k,m,Dy,DW1_N,DW2_N,theta,theta_idx)
 

Phi1=theta(1:theta_idx(1),1);
Phi2=[Phi1(1,1),Phi1(2,1) ;Phi1(3,1), Phi1(4,1) ];

omega11=theta(theta_idx(1)+1:theta_idx(2),1);
omega1=[omega11(1,1), omega11(3,1); omega11(3,1), omega11(2,1)];


Sigmae1=theta(theta_idx(2)+1:theta_idx(3),1);
Sigmae=[Sigmae1(1,1), Sigmae1(3,1); Sigmae1(3,1), Sigmae1(2,1)];


d1 =theta(theta_idx(3)+1:theta_idx(3)+1+m,1);
d2=   theta(theta_idx(3)+2+m:theta_idx(4),1);
d=[d1;d2]; % T
 

%Omega=2*eye(T1*(m+1))+diag(-1*ones((T1-1)*(m+1),1),m+1)+diag(-1*ones((T1-1)*(m+1),1),-(m+1)); % T(m+1) byT(m+1)
Omega=2*eye(T1)+diag(-1*ones((T1-1),1),1)+diag(-1*ones((T1-1),1),-1);

%L1=ones(T1); 
%L= tril(L1);
%e1=ones(T1,1);
psi=((Sigmae)^(-1))*omega1;
%invSigmaE=pinv((kron(pinv(L'*L),eye(m+1)) )+(kron((e1*e1'), (psi-eye(m+1)))))*kron(eye(T1), pinv(Sigmae));

SigmaE=kron(Omega,Sigmae);
BigOmega=kron(Omega,eye(m+1));
BigOmega(1:m+1,1:m+1)=psi;
SigmaE(1:m+1,1:m+1)=omega1;  % W
%invSigmaE=pinv(SigmaE);
sqrSigmaE = SigmaE^(-0.5); %try inv


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 Dd_hat_N = zeros(T1*(m+1),N);
    for j=1:size(d,1)
       Dd_hat_N = Dd_hat_N +  squeeze(DW1_N(:,j,:))*d(j);
    end
 
 Dy_hat_N = zeros(T1*(m+1),N);
    for jj=2:2:T1*(m+1)
       Dy_hat_N(jj-1:jj,:) = Phi2*DW2_N(jj-1:jj,:);
    end    
    
    
 
    
    
    xsi_N = Dy - Dy_hat_N- Dd_hat_N; % T0*(m+1) by N   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C_N=(xsi_N*xsi_N')/N;  % T(m+1) by T(m+1)
D_N=(sqrSigmaE)*C_N*(sqrSigmaE);
%B = balance(D_N);
[V,D]=eig(D_N); % Tm by k or P=pca(D_N) P=pca(Omega^(-0.5)*Q*A^(-0.5))

[Q,L1]=sortem(V,real(D));

lambda=diag(L1); 


%logL=-(1/2)*(log(det(SigmaE)))-(1/2)*sum(log(lambda(1:k*(m+1))))+(1/2)*(sum(lambda(1:k*(m+1)))-k*(m+1))-(1/2)*sum(lambda(1:(T1+1)*(m+1)));
%logL=-(1/2)*(abs(log(det(SigmaE))))-(1/2)*sum(log(lambda(1:k)))+(1/2)*(sum(lambda(1:k))-k)-(1/2)*sum(lambda(1:(T1+1)));

logL=-((T1)/2)*(real(log(det(Sigmae))))-(1/2)*(real(log(det(BigOmega))))-(1/2)*sum(log(lambda(1:k)))+(1/2)*(sum(lambda(1:k))-k)-(1/2)*sum(lambda(1:(T1+1)));

%logL=-(T1/2)*(log(det(Sigmae)))-(1/2)*(log(det(SigmaE)))-(1/2)*sum(log(lambda(1:k)))+(1/2)*(sum(lambda(1:k))-k)-(1/2)*sum(lambda(1:T1));

logL = -logL;
% abs()imag() angle()

end