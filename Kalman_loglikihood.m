%Kalman loglikelihood
%Yt=dt+Zt*Xt+epsilon
%X(t)=ct+Qt*X(t-1)+eta

function [likelihoodsum]=Kalman_loglikihood(theta,data,T)
%data: m by n, each column of data is a measurement, T: m by 1
data=log(data);
[row,column]=size(data);

%parameters
kappa=theta(1);
mu=theta(2);
sigma=theta(3);
lamda=theta(4);
Xi=theta(5:4+row);

XiMatrix=diag(Xi.^2);

alpha=mu-sigma^2/2/kappa;
alphastar=alpha-lamda;

dt=(1-exp(-kappa*T))*alphastar+sigma^2/(4*kappa)*(1-exp(-2*kappa*T));
Zt=exp(-kappa*T);

dt=7/365;
eta=sigma^2*dt;

ct=kappa*alpha*dt;
Qt=1-kappa*dt;

%initialize
X0=mean((mean(data,2)-dt)./Zt);
P0=eta;

X_pos=X0;
P_pos=P0;

likelihoodsum=0;
All_x=zeros(1,column);
All_y=zeros(1,column);

for i=1:column
    All_x(i)=X_pos;
    %update time
    X_pri=ct+Qt*X_pos;
    P_pri=Qt*P_pos*Qt+eta;
    %measurement update
    cov=XiMatrix+Zt*P_pri*Zt';
    K=P_pri*Zt'*inv(cov);
    Yerror=data(:,i)-(dt+Zt*X_pri);
    X_pos=X_pri+K*Yerror;
    P_pos=(1-K*Zt)*P_pri;
    
    dumy=(dt+Zt*X_pri);
    All_y(i)=dumy(1);
    likelihoodsum=likelihoodsum+0.5*log(2*pi)+0.5*log(det(cov))+0.5*Yerror'*inv(cov)*Yerror;
    
end

plot(All_y,'.-')
hold on
plot(data(1,:),'r-')
plot(All_x,'k-')
likelihoodsum

end