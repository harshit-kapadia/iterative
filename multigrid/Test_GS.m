% Script to test Gauss-Seidel Smoother
% NOTE: GS_test() is used for the test as iterations, nu, are not fixed
% Author: Harshit J. Kapadia

clc ;
clear all ;
close all ;

N_values = [10 100] ;

for N=N_values
    
    h = 1.0/N ;
    
    u = zeros(N+1,N+1) ; % as index varies from 1 to N+1 instead of 0 to N
    u_exact = zeros(N+1,N+1) ;
    f = zeros(N+1,N+1) ;
    
    for i=1:N+1
        for j=1:N+1
            u_exact(i,j) = sin(2*pi*(i-1)*h) * sin(2*pi*(j-1)*h) ;
            f(i,j) = 8*pi*pi * sin(2*pi*(i-1)*h) * sin(2*pi*(j-1)*h) ;
        end
    end
    
    u = GS_test(u,f,N) ;
    
    error = u - u_exact ;
    error_norm = norm_inf(error,N) ;
    
    fprintf('Converged maximum error for N=%d: %f \n\n',N,error_norm) ;
    
end
