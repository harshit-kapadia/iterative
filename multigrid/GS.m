function u = GS(u,f,nu,N)
% Gauss-Seidel Smoother for Multigrid Algorithm (GS-LEX)
h = 1.0/N ;
for k=1:nu+1 % instead of 0 to nu, shifting the range by 1
    for j=2:N % as actual range is from 1 to N+1 in each space direction
        for i=2:N
            u(i,j) = (1/4) * ( h*h*f(i,j) + u(i-1,j) + u(i,j-1) + u(i,j+1) + u(i+1,j) ) ;
        end
    end
end
end

