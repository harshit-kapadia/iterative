function g = RESTR(u,Nc)
% Restriction Operator for Multigrid Algorithm
g = zeros(Nc+1,Nc+1) ;
for i=2:Nc % as actual range is from 1 to Nc+1 in each space direction
    ii = 2*i - 1 ; % can't be 2*i as index varies from 1 to N+1 this Matlab code
    for j=2:Nc
        jj = 2*j - 1 ;
        g(i,j) = (1/16) * ( u(ii-1,jj-1) + 2*u(ii,jj-1) + u(ii+1,jj-1)...
            + 2*u(ii-1,jj) + 4*u(ii,jj) + 2*u(ii+1,jj)...
            + u(ii-1,jj+1) + 2*u(ii,jj+1) + u(ii+1,jj+1)) ;
    end
end
end

