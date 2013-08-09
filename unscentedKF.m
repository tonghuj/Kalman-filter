%unscented Kalman filter  for reference " THE SQUARE-ROOT UNSCENTED KALMAN FILTER
%FOR STATE AND PARAMETER-ESTIMATION ", and " the unscented Kalman filter for nonlinear estimation"

% for nonlinear esitmation
%y=func(x,pv)
%y: m by 1
%x: n by 1
%pv: l by 1

%pn: m by 1
%measure: m by steps

%x,pv,pn are vertical vectors

function unscentedKF(func,x,px,pv,pn,measurex,measurey)

lenx=length(x);
lenv=length(pv);
lenn=length(pn);
[dummy,steps]=size(measurey);

Xa=[x;zeros(lenv,1);zeros(lenn,1)];
Pa=blkdiag(eye(length(x))*px,diag(pv),diag(pn));

X_pos=x;
P_pos=eye(length(x))*px;

X_all=zeros(lenx,steps);

alpha=0.001;
beta=2;
kappa=0;
L=length(Xa);

limda=alpha^2*(L+kappa)-lenx;

for i=1:steps
    
    X_all(:,i)=X_pos;
    
    %calculate sigma points
    root=chol(Pa*(L+limda));
    root
    weight_m0=limda/(L+limda);
    weight_c0=limda/(L+limda)+(1-alpha^2+beta);
    weight=ones(2*L,1)*1/(2*(L+limda));
    weight
    BigX=Xa*ones(1,2*L+1);
    BigX(:,2:L+1)=BigX(:,2:L+1)-root;
    BigX(:,L+2:2*L+1)=BigX(:,L+2:2*L+1)+root;
    BigX
    %time update
    X_pri=X_pos;
    P_pri=P_pos;
    ytemp=func(BigX(1:lenx,:),BigX(lenx+lenv+1:lenx+lenv+lenn,:),measurex(:,i));
    y_pri=weight_m0*ytemp(:,1)+(ytemp(:,2:2*L+1))*weight;
    ytemp
    y_pri
    %measurement update
    pyy=weight_c0*(ytemp(:,1)-y_pri)*(ytemp(:,1)-y_pri)';
    pyy=pyy+(ytemp(:,2:end)-y_pri*ones(1,2*L))*diag(weight)*(ytemp(:,2:end)-y_pri*ones(1,2*L))';
    pyy
    pxy=weight_c0*(BigX(1:lenx,1)-X_pri)*(ytemp(:,1)-y_pri)';
    pxy=pxy+(BigX(1:lenx,2:end)-X_pri*ones(1,2*L))*diag(weight)*(ytemp(:,2:end)-y_pri*ones(1,2*L))';
    pxy
    K=pxy*inv(pyy);
    X_pos=X_pri+K*(measurey(:,i)-y_pri);
    P_pos=P_pri-K*pyy*K';
    K
    X_pos
    P_pos
    %update Xa,Pa
    Xa(1:lenx)=X_pos;
    Pa(1:lenx,1:lenx)=P_pos;
    
end

for i=1:lenx
plot(X_all(i),'.-')
hold on
end

end