function r = residual(u,f,N)
% Computing residual of Poisson equation
h = 1.0 / N ;
r = zeros(N+1,N+1) ;
for i=2:N
    for j=2:N
        r(i,j) = f(i,j) + ( (u(i+1,j) - 2*u(i,j) + u(i-1,j))/h*h...
            + (u(i,j-1) - 2*u(i,j) + u(i,j+1))/h*h ) ;
    end
end
end

