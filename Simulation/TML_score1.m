function score_N = TML_score1(N,T1,k,m,Dy,Dy1, coef_ML, theta_idx)

theta=coef_ML;
tiny = 0.00001;

dim_theta = size(theta,1);


score_N = zeros(dim_theta,N);

for i=1:N

logL_i = @(theta)Loglikelihood_fullsigome_nt(N,T1,k,m,Dy,Dy1,theta,theta_idx); 
for K=1:length(theta)
    theta1 = theta;  theta1(K) = theta(K) - tiny;
    theta2 = theta;  theta2(K) = theta(K) + tiny;
    if K==1 
        score_i = ( (logL_i(theta2) - logL_i(theta1))'/(2*tiny) )';
    else
        score_i = [score_i, ((logL_i(theta2) - logL_i(theta1))'/(2*tiny))'];
    end;
end;
score_N(:,i) = score_i;

end
end
