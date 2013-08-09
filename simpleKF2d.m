%simple 2d Kalman filter
% X*c=y
%[x1 x2][c1;c2]=[y]
% so H=X, 1*2 matrix


function simpleKF2d(steps)


measurex1=10:1/steps:11;
measurex2=2:1/steps:3;

measurex1=measurex1+normrnd(0,0.1,[1,steps+1]); %add noise
measurex2=measurex2+normrnd(0,0.1,[1,steps+1]);

measurey=measurex1*1.0+measurex2*2.0+normrnd(0,0.3,[1,length(measurex1)]);

measurex1_all=zeros(1,steps);
measurex2_all=zeros(1,steps);

X_pri=[0;0];
X_pos=[0;0];
P_pri=eye(2,2);
P_pos=eye(2,2);

Q=eye(2,2)*0.00001;
R=0.2;

for i=1:steps
    
    measurex1_all(i)=X_pos(1);
    measurex2_all(i)=X_pos(2);
    
    %update time
    X_pri=X_pos;
    P_pri=P_pos+Q;
    %update measure
    H=[measurex1(i) measurex2(i)];
    K=P_pri*H'*inv(H*P_pri*H'+R);
 
    X_pos=X_pri+K*(measurey(i)-H*X_pri);
    P_pos=(eye(2,2)-K*H)*P_pri;
    
end

plot(measurex1_all,'r.-')
hold on
plot(measurex2_all,'k.-')

figure
scatter3(measurex1,measurex2,measurey)

figure
hist(measurex1_all,100);
hold on
hist(measurex2_all,100);

end