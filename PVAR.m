tic

rep = 300;         % number of replications
num_ini=7; 
list_T = [10];     % list of time, T
list_N = [100 300 500];   % list of N   
 
 
bias_mean_gam_11=zeros (size(list_T,2), size(list_N,2));  
bias_mean_gam_12=zeros (size(list_T,2), size(list_N,2));  
bias_mean_gam_21=zeros (size(list_T,2), size(list_N,2));  
bias_mean_gam_22=zeros (size(list_T,2), size(list_N,2));  
 
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
bias_mean_sig_22=zeros (size(list_T,2), size(list_N,2));  
bias_mean_sig_12=zeros (size(list_T,2), size(list_N,2));  
 
 
 
 
 
 
m=1;  % number of endogenous variables=m+1
k=1; % number of factor
Dis_T=60;   % discard first time period
Phi =[0.4, 0.2 ; 0.2, 0.4]; % Panel VAR(1) coefficient    
Phi_11=Phi(1,1);
Phi_21=Phi(2,1);
Phi_12=Phi(1,2);
Phi_22=Phi(2,2);
 
mu=[0 0]; % individual mean  
sigmau=[1, 0 ;0, 1];      % individual Variance
sigma_idi=[0.07, 0.05; 0.05, 0.07];
R1 = chol(sigmau);
R2= chol(sigma_idi);
 
Mueta=[0 0 0 0];   % factor loading
a=(1/k)*ones(1,(m+1)*(k));   
sigmaeta=diag(a);
R3=chol(sigmaeta);
 

for idx_T=1:size(list_T,2)
T0 = list_T(idx_T);            
 
 
TT= (T0+1)+Dis_T; 
T1=T0-1;  
 
for idx_N=1:size(list_N,2)      
N = list_N(idx_N);       
 
randn('state', 12345678) ;
rand('state', 1234567) ;
    RandStream.setGlobalStream (RandStream('mcg16807','seed',34));
 
 
delta=zeros(TT,1);
for ttt=1:TT
delta(ttt,1)=((ttt^2-ttt)/2);  % normalize time effect
end

f=zeros(k,TT);  % creat a space for saving factor data
for t=2:TT
   f(:,t)= 0.9*f(:,t-1)+sqrt(1-0.5^2)*(normrnd(0,1)); % k by TT
end 


 
 
 
 
est_gam_nsm11=zeros(rep,1); 
est_gam_nsm12=zeros(rep,1); 
est_gam_nsm21=zeros(rep,1); 
est_gam_nsm22=zeros(rep,1);  
 
 
est_sig_nsm11=zeros(rep,1); 
est_sig_nsm22=zeros(rep,1); 
est_sig_nsm12=zeros(rep,1); 
 

randn('state', 12345678) ;
rand('state', 1234567) ;
   RandStream.setGlobalStream (RandStream('mcg16807','seed',34));
 
 
skip=0;
sml=1;
while sml<=rep
 

 

fac=(mvnrnd(mu,sigmau,N))';

 
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

 
 
 
 

 
y=zeros(m+1,N,TT);  
y(:,:,1)=zeros(m+1,N);         
for tt=2:TT
  for ii=1:N     
   y(:,ii,tt)=Phi*y(:,ii,tt-1)+indi_u(:,ii)+(delta(tt))*ones(m+1,1)+(fac(:,ii)*f(:,tt)) +(idi_error(ii,:,tt))';  % gen data m+1 by N by T
   end
end
Y_NT=y(:,:,TT-T0:TT); 
 
y = Y_NT(:,:,:); 
 

 
Dy=zeros(m+1,N,T0);
for dy= 1:T0
 Dy(:,:,dy)=y(:,:,dy+1)-y(:,:,dy);   
end

 
DW1_N = zeros(T0*(m+1),T0*(m+1),N); % T x dim_phi x N   for d
    for i=1:N
       DW1_N(:,:,i) =  eye(T0*(m+1)); 
    end
 
tranDy=[];
for kt=1:T0   
tranDy= cat(1,tranDy,Dy(:,:,kt));    % combine Dy (m+1)(T0) by N
end 
Dy=tranDy;  % (m+1)(T0) by N    
 
    
DW2_N = zeros(T0*(m+1),N); % T x dim_phi x N   for y_t-1
    for ii=1:N
       DW2_N(:,ii) =  [zeros(m+1,1);Dy(1:(T1)*(m+1),ii)]; 
    end    
    
    

dim_phi =(m+1)*(m+1);
dim_ome = ((m+2)*(m+1)/2);
dim_sig2 = ((m+2)*(m+1)/2);
dim_d = (T0)*(m+1); %   T0*(m+1)
 
theta_idx = cumsum([dim_phi; dim_ome; dim_sig2; dim_d]);
 
dim_theta = theta_idx(end);
 
 

options_ML = optimoptions(@fmincon, 'Display','off','GradObj','off','Hessian','off','Algorithm','interior-point' ,'TolX',10^(-6),'TolFun',10^(-6),'TolCon',10^(-5), 'MaxIter',600);
 
 
  lb=[ -1*ones((m+1)*(m+1),1)      ;[0.0001; 0.0001;-inf]       ;[0.0001; 0.0001; -inf]     ; -inf*ones((T0)*(m+1),1)];
 ub=[ 1*ones((m+1)*(m+1),1)    ; inf*ones(3,1)          ; inf*ones(3,1)   ; inf*ones((T0)*(m+1),1) ];
 
 
 
logL_num_ini = zeros(2,num_ini); % *********
coef_num_ini = zeros(dim_theta,1,num_ini);
 
for j=1:num_ini   
 
gam_ini1=-0.5+rand(1,(m+1)*(m+1));
 
omega_ini1=rand(1,2)+1;
omega_ini2=-1+2*rand(1,1); 
omega_ini3=[ omega_ini1, omega_ini2];
 
 
sigma_ini1=rand(1,2)+1;
sigma_ini2=-1+2*rand(1,1);
sigma_ini3=[sigma_ini1,sigma_ini2];
  
d_ini=rand(1,m+1)*2-1;
tdum_ini1=rand(1,(T1)*(m+1))*2-ones(1,(T1)*(m+1));
 
theta_ini=[ gam_ini1,  omega_ini3, sigma_ini3,d_ini,tdum_ini1]';
 
 
[theta_ML_j,fval_j,exitflag_j,~,~,~,~] = fmincon(@(theta)Loglikelihood_fullsigomes(N,T0,T1,k,m,Dy,DW1_N,DW2_N,theta,theta_idx),theta_ini,[],[],[],[],lb,ub,'constraints_fullsigome',options_ML);
 
 

 
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
 elseif indx==1;  
        ss=size(coef_num_ini,1); 
          coef_ML = NaN*ones(ss,m+1);     
 end
 
 

coef_ML1= coef_ML(1:theta_idx(1) ,1);
omega_ML= coef_ML( theta_idx(1)+1:theta_idx(2) ,1);
sigma_ML= coef_ML(theta_idx(2)+1:theta_idx(3) ,1);
d_ML = coef_ML(theta_idx(3)+1:theta_idx(4),1);
 
est_gam_nsm11(sml,1)  =coef_ML1(1,1);  % saving gamma in est_gam_nsml
est_gam_nsm12(sml,1)  =coef_ML1(2,1);  % saving gamma in est_gam_nsml
est_gam_nsm21(sml,1)  =coef_ML1(2,1);  % saving gamma in est_gam_nsml
est_gam_nsm22(sml,1)  =coef_ML1(4,1);  % saving gamma in est_gam_nsml
 
 
 
est_sig_nsm11(sml,1)  =sigma_ML(1,1);  % saving gamma in est_gam_nsml
est_sig_nsm22(sml,1)  =sigma_ML(2,1);  % saving gamma in est_gam_nsml
est_sig_nsm12(sml,1)  =sigma_ML(3,1);  % saving gamma in est_gam_nsml
 
 
sml=sml+1;
 
 
indxML=isnan(coef_ML);
 
  if sum(indxML)~=0 ;
    sml=sml-1;
    skip=skip+1;
    [skip sml ];
continue;
  end
 
end
 
 
 
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
bias_mean_sig_12(idx_T, idx_N)=mean_sig11-sigma_idi(1,2);
 
 
 
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
 
end
end
toc
