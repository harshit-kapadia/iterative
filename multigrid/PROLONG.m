function g = PROLONG(u,Nc)
%UNTITLED2 Summary of this function goes here
g = zeros(2*Nc+1,2*Nc+1) ;
for i=2:Nc
    ii = 2*i - 1 ; % can't be 2*i as index varies from 1 to N+1 this Matlab code
    for j=2:Nc
        jj = 2*j - 1 ;
        
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

