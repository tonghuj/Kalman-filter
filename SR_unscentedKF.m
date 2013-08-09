%square root unscented kalman filter
% reference "THE SQUARE-ROOT UNSCENTED KALMAN FILTER
% FOR STATE AND PARAMETER-ESTIMATION "

function SR_unscentedKF(func,w0,pw0,pe,measurex,measurey)

[y_dim steps]=size(measurey);

Re=eye(length(y_dim))*pe;
root_Re=chol(Re);

Sw0=chol(eye(length(w0))*pw0);
Sw_pos=Sw0;

w_pos=w0;
w_all=zeros(length(w0),steps);


limda_RLS=0.995;
alpha=0.001;
kappa=0;
L=length(w0);
limda=alpha^2*(L+kappa)-L;
gamma=sqrt(L+limda);
beta=2;

%construct weight
weight=[limda/(L+limda);ones(2*L,1)*1/(2*(L+limda))];
weightc=limda/(L+limda)+(1-alpha^2+beta);


for i=1:steps
    w_pos;
    w_all(:,i)=w_pos;
    %time update
    w_pri=w_pos;
    Sw_pri=1/sqrt(limda_RLS)*Sw_pos;
    
    BigW=[w_pri w_pri*ones(1,L)+gamma*Sw_pri w_pri*ones(1,L)-gamma*Sw_pri];
    D=func(measurex(:,i),BigW);
    
    d=D*weight;
    
    %measure update
    
    Sd=qr([sqrt(weight(1))*(D(:,2:2*L+1)-d*ones(1,2*L)) root_Re]');
    cut=min(size(Sd));
    Sd=Sd(1:cut,1:cut);
    
    if(sign(weightc)>=0)
        Sd=cholupdate(Sd,(abs(weightc))^(1/4)*(D(:,1)-d),'+');
    else
        Sd=cholupdate(Sd,(abs(weightc))^(1/4)*(D(:,1)-d),'-');
    end
    
    Pwd=(BigW-w_pri*ones(1,2*L+1))*diag(weight)*(D-d*ones(1,2*L+1))';
    K=(Pwd/(Sd'))/Sd;
    measurey(:,i)-d
    w_pos=w_pri+K*(measurey(:,i)-d);
    U=K*Sd;
    A=Sw_pri*Sw_pri'-U*U'+eye(length(w0))*0.00001;
    Sw_pos=chol(A);
    %Sw_pos=cholupdate(Sw_pri,U,'-');
    
end

plot(w_all(1,:),'.-r')
hold on
plot(w_all(2,:),'.-k')

end