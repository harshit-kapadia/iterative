function u = GS_test(u,f,N)
% Gauss-Seidel Smoother for Multigrid Algorithm (GS-LEX)
h = 1.0/N ;

nu = 0 ; % iteration counter
diff_norm = 1 ; % dummy value for 1st iteration
while diff_norm > 1e-10
    nu = nu + 1 ;
    u_old = u ;
    for j=2:N % as actual range is from 1 to N+1 in each space direction
        for i=2:N
            u(i,j) = (1/4) * ( h*h*f(i,j) + u(i-1,j) + u(i,j-1) +...
                                              u(i,j+1) + u(i+1,j) ) ;
        end
    end
    diff = u - u_old ;
    diff_norm = norm_inf(diff,N) ;
end
fprintf('Maximum Iteration for N=%d\n',N) ;
display(nu) ;
end

