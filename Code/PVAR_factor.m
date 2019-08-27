function [  max_logL, coef_ML ] = PVAR_factor( DATA, num_fac, num_ini, time_effect  )


[T0, N, KK]  = size(DATA);
T1 = T0-1;
m = KK-1;

k = num_fac;  %number of factor


Dy  = DATA(2:T0,:) - DATA(1:T0-1,:) ; % T x N (Dy_1,...,Dy_T)       

DY_NT=[];
for t=1:T1
for mm=1:KK    
DY_NT=DY_NT+cat(1,Dy(t,:,mm)); % Tm by  N
end
end




DW1_N = zeros(T0*(m+1),T0*(m+1),N); % Tm by Tm by N   for d
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
    
    
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
dim_phi =(m+1)*(m+1);
dim_ome = ((m+2)*(m+1)/2);
dim_sig2 = ((m+2)*(m+1)/2);
dim_d = (T0)*(m+1); %   T0*(m+1)
 
theta_idx = cumsum([dim_phi; dim_ome; dim_sig2; dim_d]);
 
dim_theta = theta_idx(end);
 
 
options_ML = optimoptions(@fmincon, 'Display','off','GradObj','off','Hessian','off','Algorithm','interior-point' ,'TolX',10^(-6),'TolFun',10^(-6),'TolCon',10^(-5), 'MaxIter',600);
 
 
  lb=[ -1*ones((m+1)*(m+1),1)      ;[0.0001; 0.0001;-inf]       ;[0.0001; 0.0001; -inf]     ; -inf*ones((T0)*(m+1),1)];
 ub=[ 1*ones((m+1)*(m+1),1)    ; inf*ones(3,1)          ; inf*ones(3,1)   ; inf*ones((T0)*(m+1),1) ];
 
 
 
logL_num_ini = zeros(2,num_ini);
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
 

 
[min_logL, indx_maxLL] = min(logL_num_ini(1,:));
max_logL = -min_logL;
indx=isnan(max_logL);
 
 if indx==0   
        
    coef_ML = coef_num_ini(:,:,indx_maxLL); 
 elseif indx==1;  
        ss=size(coef_num_ini,1); 
          coef_ML = NaN*ones(ss,m+1);     
 end
 
 
%sum(isnan(logL_num_ini'))
 
coef_ML1= coef_ML(1:theta_idx(1) ,1);
omega_ML= coef_ML( theta_idx(1)+1:theta_idx(2) ,1);
sigma_ML= coef_ML(theta_idx(2)+1:theta_idx(3) ,1);
d_ML = coef_ML(theta_idx(3)+1:theta_idx(4),1);
    
    
    
    
    
    
    
    
end
    