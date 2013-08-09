

function y = nonlinearfunc(x,w)

%x: n by 1
%w: m by 2n+1
%y: l by 2n+1

x=[x(1)^2 x(2)^3];
y=x*w;

end