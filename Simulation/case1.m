tic
% Case 1 (Covariance stationary and stable ;Assume factor is 2 )
rep = 50;         % number of replications
num_ini=6; % ***********************
list_T = [5 6 7];     % list  T
list_N = [100 300 500];   % list  N   
list_power = [0.00 0.10]; % list size and power

bias_mean_gam_11=zeros(size(list_T,2), size(list_N,2));  
bias_mean_gam_12=zeros(size(list_T,2), size(list_N,2));  
bias_mean_gam_21=zeros(size(list_T,2), size(list_N,2));  
bias_mean_gam_22=zeros(size(list_T,2), size(list_N,2));  

mean_specnorm= zeros (size(list_T,2), size(list_N,2)); 
mean_frobnorm=zeros (size(list_T,2), size(list_N,2)); 

MAE_gam11=zeros(size(list_T,2), size(list_N,2));    
MAE_gam21=zeros(size(list_T,2), size(list_N,2));   
MAE_gam12=zeros(size(list_T,2), size(list_N,2));   
MAE_gam22=zeros(size(list_T,2), size(list_N,2));   


std_gam11=zeros(size(list_T,2), size(list_N,2));        
std_gam12=zeros(size(list_T,2), size(list_N,2));        
std_gam21=zeros(size(list_T,2), size(list_N,2));        
std_gam22=zeros(size(list_T,2), size(list_N,2));        

rmse_gam11=zeros(size(list_T,2), size(list_N,2));      
rmse_gam21=zeros(size(list_T,2), size(list_N,2));      
rmse_gam12=zeros(size(list_T,2), size(list_N,2));       
rmse_gam22=zeros(size(list_T,2), size(list_N,2));       

bias_mean_sig_11=zeros (size(list_T,2), size(list_N,2));  
bias_mean_sig_21=zeros (size(list_T,2), size(list_N,2));
bias_mean_sig_22=zeros (size(list_T,2), size(list_N,2));  


se_gam_11=zeros(size(list_T,2), size(list_N,2));  
se_gam_12=zeros(size(list_T,2), size(list_N,2));  
se_gam_21=zeros(size(list_T,2), size(list_N,2));  
se_gam_22=zeros(size(list_T,2), size(list_N,2)); 


m=1;  % number of regressor
k=2; % number of factor
Dis_T=60;   % discard first time period
Phi =[0.2, 0.1 ; 0.1, 0.2]; % Panel VAR(1) coefficient    
Phi_11=Phi(1,1);
Phi_12=Phi(1,2);
Phi_21=Phi(2,1);
Phi_22=Phi(2,2);

mu=[0 0]; % individual mean  
sigmau=[1, 0 ;0, 1];      % individual Variance
sigma_idi=[1, 0.5; 0.5, 1];
sigma_idi11=sigma_idi(1,1);
sigma_idi12=sigma_idi(1,2);
sigma_idi22=sigma_idi(2,2);

sigmaul=[1/k, 0 ;0, 1/k];  
R1 = chol(sigmau);
R2= chol(sigma_idi);
Mueta=[0 0 0 0];   % factor loading
a=(1/k)*ones(1,(m+1)*(k));   
sigmaeta=diag(a);                  
R3=chol(sigmaeta);            

%se_gam=zeros(size(list_T,2), size(list_N,2));         % creat a space for saving bias in different setting
for idx_T=1:size(list_T,2)
T0 = list_T(idx_T);        % for loop of T     

TT= (T0+1)+Dis_T; % T* (because we want to discard first 50 units)
T1=T0-1;  

for idx_N=1:size(list_N,2)      
N = list_N(idx_N);        % for loop of N 

randn('state', 12345678) ;
rand('state', 1234567) ;
   RandStream.setGlobalStream (RandStream('mcg16807','seed',34));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=zeros(k,TT);  % gen factor
for t=2:TT
   f(:,t)= 0.5*f(:,t-1)+sqrt(1-0.5^2)*(repmat(mu,1,1) + randn(1,2)*R1)'; % k by TT  
end 

f1=f;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



est_gam_nsm11=zeros(rep,1); % creat the space for saving result 
est_gam_nsm12=zeros(rep,1); % creat the space for saving result 
est_gam_nsm21=zeros(rep,1); % creat the space for saving result 
est_gam_nsm22=zeros(rep,1); % creat the space for saving result 
 specnorm=zeros(rep,1);
 frobnorm=zeros(rep,1);
coef_Phi=zeros(m+1,m+1,rep);
D=zeros(m+1,m+1,rep);
est_sig_nsm11=zeros(rep,1); % creat the space for saving result 
est_sig_nsm22=zeros(rep,1); % creat the space for saving result 
est_sig_nsm12=zeros(rep,1); % creat the space for saving result 


se_gam_nsml1=zeros(rep,1);
se_gam_nsml2=zeros(rep,1);
se_gam_nsm21=zeros(rep,1);
se_gam_nsm22=zeros(rep,1);


size_power_gam_nsm11 = zeros(size(list_power,2),rep);
size_power_gam_nsm12 = zeros(size(list_power,2),rep);
size_power_gam_nsm21 = zeros(size(list_power,2),rep);
size_power_gam_nsm22 = zeros(size(list_power,2),rep);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

randn('state', 12345678) ;
rand('state', 1234567) ;
   RandStream.setGlobalStream (RandStream('mcg16807','seed',34));


skip=0;
sml=1;
while sml<=rep


u=repmat(mu,N,1) + randn(N,2)*R1;   % individual effect error (N by m+1) N(0,1)
idi_error=zeros(N,m+1,TT);
for id=1:TT
idi_error(:,:,id)=repmat(mu,N,1) + randn(N,2)*R2; 
end
mean_idie=zeros(N,m+1);
for mid=1:TT
mean_idie=mean_idie+idi_error(:,:,mid);
end
mean_idie=(mean_idie/TT);  %mean of idi error N by m+1


b0=1;
b1=1;
indi_u=(b0*mean_idie+b1*u)'; % gen idividual effect  m+1 by N

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eta=(repmat(Mueta,N,1) + randn(N,4)*R3);   % gen factor loading
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y=zeros(m+1,N,TT);  
y(:,:,1)=zeros(m+1,N);       

for tt=2:TT
    fac=zeros(m+1,N); 
  for ii=1:N     
    fac(:,ii)=(reshape(eta(ii,:),[m+1,k]))*f1(:,tt);  % gen factor loading*factor; m+1 by N  
  end
    y(:,:,tt)=Phi*y(:,:,tt-1)+indi_u+(fac)+(idi_error(:,:,tt))';  % gen data m by N by T
end
Y_NT=y(:,:,TT-T0:TT); % discard first 50 observation 
y = Y_NT(:,:,:); 


Dy=y(:,:,2:T0+1)-y( :, :,1:T0);  % gen data m by N by T-1  (Dy1,...,DyT)
Dy1=Dy(:,:,1:T1); % m by N by T1 (Dy1,...,DyT-1)


tranDy=[];
for kt=1:T0   
tranDy= cat(1,tranDy,Dy(:,:,kt));    % combine Dy (m+1)(T0) by N
end 
Dy=tranDy;  % (m+1)*(T0) by N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim_phi =((m+1)*(m+1)); 
dim_ome = ((m+2)*(m+1)/2);
dim_sig2 = ((m+2)*(m+1)/2);

theta_idx = cumsum([dim_phi; dim_ome; dim_sig2]);

dim_theta = theta_idx(end);



options_ML = optimoptions(@fmincon, 'Display','off','GradObj','off','Hessian','off','Algorithm','interior-point' );


 lb=[ -1*ones(2*(m+1),1)      ;[0.0001; 0.0001;-inf]       ;[0.0001; 0.0001; -inf]     ];
 ub=[ 1*ones(2*(m+1),1)    ; inf*ones(3,1)          ; inf*ones(3,1)    ];
 
 

logL_num_ini = zeros(2,num_ini); % *********
coef_num_ini = zeros(dim_theta,1,num_ini);

for j=1:num_ini   

    
gam_ini1=[Phi(1,1), Phi(1,2), Phi(2,1), Phi(2,2)];  % initial true value of phi
 omega_ini1=rand(1,2)+1;
 omega_ini2=-1+2*rand(1,1);
omega_ini3=[ omega_ini1, omega_ini2];  % symmertic variance of y_{i,1}

 
 sigma_ini3=[sigma_idi11,sigma_idi22,sigma_idi12]; %  initial true value of symmertic variance of e_{i,t}


theta_ini=[ gam_ini1,  omega_ini3, sigma_ini3]';



[theta_ML_j,fval_j,exitflag_j,~,~,~,~] = fmincon(@(theta)Loglikelihood_fullsigome_nt(N,T1,k,m,Dy,Dy1,theta,theta_idx),theta_ini,[],[],[],[],lb,ub,'constraint_fullsigome_nt',options_ML);


% theta_ML_j
% fval_j
% exitflag_j

if exitflag_j>0
  logL_num_ini(:,j)  = [fval_j; exitflag_j];
  coef_num_ini(:,:,j)  = theta_ML_j;
  
else
  logL_num_ini(:,j)  = NaN;
  coef_num_ini(:,:,j)  = NaN;
end  


end    

sml
logL_num_ini

[min_logL, indx_maxLL] = min(logL_num_ini(1,:));
max_logL = -min_logL;
indx=isnan(max_logL);

 if indx==0   
        
    coef_ML = coef_num_ini(:,:,indx_maxLL); 
    eps=1e-5;
       
    Hessian2 = compute_Hessian1(@(theta)Loglikelihood_fullsigome_nt(N,T1,k,m,Dy,Dy1,theta,theta_idx), coef_ML, eps);
  min_eig = min(eig((Hessian2)));
       
        if min_eig <=0   
        flag=1;
        end 
                
        score_N = TML_score1(N,T1,k,m,Dy,Dy1, coef_ML, theta_idx);
        Hessian3 =(score_N*score_N')*N;

        minHes3unres=min(eig(Hessian3));
     
        var5_ML = pinv(Hessian2)*Hessian3*pinv(Hessian2)/N;
        
        se5_ML = abs(sqrt(diag(var5_ML)));
    
 elseif indx==1;  
        ss=size(coef_num_ini,1); 
          coef_ML = NaN*ones(ss,m+1);     
 end
 



coef_ML1= coef_ML(1:theta_idx(1) ,1);
omega_ML= coef_ML( theta_idx(1)+1:theta_idx(2) ,1);
sigma_ML= coef_ML(theta_idx(2)+1:theta_idx(3) ,1);

est_gam_nsm11(sml,1)  =coef_ML1(1,1); 
est_gam_nsm12(sml,1)  =coef_ML1(2,1); 
est_gam_nsm21(sml,1)  =coef_ML1(3,1);  
est_gam_nsm22(sml,1)  =coef_ML1(4,1);  

coef_Phi(:,:,sml)=[coef_ML1(1,1), coef_ML1(2,1) ; coef_ML1(3,1) , coef_ML1(4,1)];
 

est_sig_nsm11(sml,1)  =sigma_ML(1,1);  
est_sig_nsm22(sml,1)  =sigma_ML(2,1);  
est_sig_nsm12(sml,1)  =sigma_ML(3,1);  

se_gam_nsml1(sml,1)  =se5_ML(1,1);  
se_gam_nsml2(sml,1)  =se5_ML(2,1);   
se_gam_nsm21(sml,1)  =se5_ML(3,1);   
se_gam_nsm22(sml,1)  =se5_ML(4,1);  

for idx_power=1:size(list_power,2)  % do test for power and size of test, power=P(reject H_0 | H_1 is true ), size=P(reject H_0 | H_0 is true) 
b = list_power(idx_power);

    if abs(est_gam_nsm11(sml,1) -Phi_11)/se_gam_nsml1(sml,1)> 1.96; 
        size_power_gam_nsm11(idx_power, sml)=1;  end
     
    if abs(est_gam_nsm12(sml,1) -Phi_12)/se_gam_nsml2(sml,1)> 1.96; 
        size_power_gam_nsm12(idx_power, sml)=1;  end
     
    if abs(est_gam_nsm21(sml,1) -Phi_21)/se_gam_nsm21(sml,1)> 1.96; 
        size_power_gam_nsm21(idx_power, sml)=1;  end
     
    if abs(est_gam_nsm22(sml,1) -Phi_22)/se_gam_nsm22(sml,1)> 1.96; 
        size_power_gam_nsm22(idx_power, sml)=1;  end
    
      
end    


sml=sml+1;

indxML=isnan(coef_ML);

  if sum(indxML)~=0 ;
    sml=sml-1;
    skip=skip+1;
    [skip sml ];
continue;
  end

end

 for sml=1:rep
  D(:,:,sml)= coef_Phi(:,:,sml)-Phi;
specnorm(sml,:) =sqrt(max(eig( D(:,:,sml)'* D(:,:,sml)))); 
frobnorm(sml,:) =sqrt(trace(D(:,:,sml)* D(:,:,sml)')); 
 end
 
 
 mean_specnorm(idx_T, idx_N)= nanmean(specnorm);
 mean_frobnorm(idx_T, idx_N)=nanmean(frobnorm);

mean_gam11= nanmean(est_gam_nsm11);       
mean_gam12= nanmean(est_gam_nsm12);       
mean_gam21= nanmean(est_gam_nsm21);       
mean_gam22= nanmean(est_gam_nsm22);       


mean_sig11= nanmean(est_sig_nsm11);      
mean_sig22= nanmean(est_sig_nsm22);       
mean_sig21= nanmean(est_sig_nsm12);       


bias_mean_gam_11(idx_T, idx_N) = mean_gam11 - Phi_11;     
bias_mean_gam_12(idx_T, idx_N) = mean_gam12 - Phi_12;     
bias_mean_gam_21(idx_T, idx_N) = mean_gam21 - Phi_21;     
bias_mean_gam_22(idx_T, idx_N) = mean_gam22 - Phi_22;    



bias_mean_sig_11(idx_T, idx_N)=mean_sig11-sigma_idi(1,1);
bias_mean_sig_22(idx_T, idx_N)=mean_sig11-sigma_idi(2,2);
bias_mean_sig_21(idx_T, idx_N)=mean_sig11-sigma_idi(1,2);



MAE_gam11(idx_T, idx_N) = sum(abs(est_gam_nsm11-Phi_11 ))/rep ;
MAE_gam12(idx_T, idx_N) = sum(abs(est_gam_nsm12-Phi_12 ))/rep ;
MAE_gam21(idx_T, idx_N) = sum(abs(est_gam_nsm21-Phi_21 ))/rep ;
MAE_gam22(idx_T, idx_N) = sum(abs(est_gam_nsm22-Phi_22 ))/rep ;


std_gam11(idx_T, idx_N) = nanstd(est_gam_nsm11);       
std_gam12(idx_T, idx_N) = nanstd(est_gam_nsm12);       
std_gam21(idx_T, idx_N) = nanstd(est_gam_nsm21);       
std_gam22(idx_T, idx_N) = nanstd(est_gam_nsm22);      


rmse_gam11(idx_T, idx_N) = sqrt( nanmean( (est_gam_nsm11-Phi_11).^2) ); 
rmse_gam12(idx_T, idx_N) = sqrt( nanmean( (est_gam_nsm12-Phi_12).^2) );  
rmse_gam21(idx_T, idx_N) = sqrt( nanmean( (est_gam_nsm21-Phi_21).^2) );  
rmse_gam22(idx_T, idx_N) = sqrt( nanmean( (est_gam_nsm22-Phi_22).^2) );  


size_power_gam11= nanmean(size_power_gam_nsm11,2);                            

size_power_gam12= nanmean(size_power_gam_nsm12,2);                            

size_power_gam21= nanmean(size_power_gam_nsm21,2);                            

size_power_gam22= nanmean(size_power_gam_nsm22,2);                           

end
end
%filename = '1.mat';
%save(filename)
toc



