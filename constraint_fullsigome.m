function [c,ceq]=constraints_fullsigome(theta_ML_j)


Phi=[theta_ML_j(1,1),theta_ML_j(2,1) ;theta_ML_j(3,1),theta_ML_j(4,1) ];
c1= max(abs(eig(Phi)))-0.999; 

 

omega1=[theta_ML_j(5,1), theta_ML_j(7,1); theta_ML_j(7,1), theta_ML_j(6,1)];
Sigmae1=[theta_ML_j(8,1), theta_ML_j(10,1); theta_ML_j(10,1), theta_ML_j(9,1)];

c2=0.0001-min(eig(Sigmae1));      % sig2
c3=0.0001-min(eig(omega1));   



c=[c1,c2,c3];

ceq=[];
