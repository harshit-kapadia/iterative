% Script to test Prolongation Operator
% NOTE: PROLONG() is tested
% Author: Harshit J. Kapadia

clc ;
clear all ;
close all ;

n_values = [4 7] ;

for n=n_values
    
    N = 2^n ;
    
    h = 1.0/N ;
    h_c = 2.0/N ;
    
    u_exact = zeros(N+1,N+1) ;
    
    u_c_exact = zeros(N/2+1,N/2+1) ;
    
    for i=1:N+1
        for j=1:N+1
            u_exact(i,j) = sin(2*pi*(i-1)*h) * sin(2*pi*(j-1)*h) ;
        end
    end
    for i=1:N/2+1
        for j=1:N/2+1
            u_c_exact(i,j) = sin(2*pi*(i-1)*h_c) * sin(2*pi*(j-1)*h_c) ;
        end
    end
    
    u_c = u_c_exact ;
    u = PROLONG(u_c,N/2) ;
    
    error = u - u_exact ;
    error_norm = norm_inf(error,N) ;
    
    fprintf('Error norm for PROLONGATION from N=%d to Nc=%d: %f \n\n',N/2,N,error_norm) ;
    
end
