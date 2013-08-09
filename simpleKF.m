% simple kalman filter to estimate a scale random variable

function simpleKF(steps)

% generate measurements: uniform random numbers

measure=normrnd(2,1,[1,steps]);
%measure=rand(1,steps);
v_all=zeros(1,steps);
%initialize

v_pri=0.1;
v_pos=0.1;

P_pri=1;
P_pos=0.5;

Q=0.001; % increase Q, measurement dominates
R=0.1;   % increase R, measurement ignored

for i=1:steps
    
    v_all(i)=v_pos;
    
    %time update
    v_pri=v_pos;
    P_pri=P_pos+Q;
    
    %measure update
    K=P_pri/(P_pri+R);
    v_pos=v_pri+K*(measure(i)-v_pri);
    P_pos=(1-K)*P_pri;
    
    
end


plot(measure,'r.-')
hold on
plot(v_all,'k.-')

end