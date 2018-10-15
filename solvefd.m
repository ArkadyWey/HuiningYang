function solvefd(n,eps,w,c,xmin,xmax,left,right,symbol)

f = @(x)10*sin(20*x)+cos(x.^5);
h = (xmax-xmin)/n;
x = xmin:h:xmax;
b = f(x)'; 
b(1) = left;
b(end) = right;
gamma = 1/(h^2);alpha = -2/(h^2);beta = 1/(h^2);zeta = -1/(2*h);theta = 1/(2*h);
e = ones(n+1,1);
A = spdiags([gamma*e,alpha*e,beta*e],[-1 0 1],n+1,n+1);
N = spdiags([zeta*e,theta*e],[-1 1],n+1,n+1);
M = c*eye(n+1);
K = eps*A + w*N + M;
K(1,1) = 1;
K(end,end) = 1;
K(1,2:end) = 0;
K(end,1:end-1) = 0;
u = K\b;
plot(x,u,symbol), axis('square')
return
