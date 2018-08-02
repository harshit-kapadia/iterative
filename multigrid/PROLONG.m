function [outputArg1,outputArg2] = PROLONG(u,Nc)
%UNTITLED2 Summary of this function goes here
g = zeros(2*Nc+1,2*Nc+1) ;
for i=2:Nc
    ii = 2*i ;
    for j=2:Nc
        jj = 2*j ;
        g(ii-1,jj-1) = g(ii-1,jj-1) + (1/4) * u(i,j) ;
        g(ii,jj-1) = g(ii,jj-1) + (1/2) * u(i,j) ;
        g(ii+1,jj-1) = g(ii+1,jj-1) + (1/4) * u(i,j) ;
        
        g(ii-1,jj) = g(ii-1,jj) + (1/2) * u(i,j) ;
        g(ii,jj) = g(ii,jj) + u(i,j) ;
        g(ii+1,jj) = g(ii+1,jj) + (1/2) * u(i,j) ;
        
        g(ii-1,jj+1) = g(ii-1,jj+1) + (1/4) * u(i,j) ;
        g(ii,jj+1) = g(ii,jj+1) + (1/2) * u(i,j) ;
        g(ii+1,jj+1) = g(ii+1,jj+1) + (1/4) * u(i,j) ;
    end
end
end

