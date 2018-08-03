function u = MG(l,u,f,gamma,nu_1,nu_2,N)
% Multigrid Algorithm
u = GS(u,f,nu_1,N) ; % smoothing and overwriting u
r = residual(u,f,N) ;
r_c = RESTR(r,N/2) ;

if l==1
    e_c = zeros(N/2+1,N/2+1) ;
    e_c = GS(e_c,-r_c,1,N/2) ; % domain has 9 grid nodes, only one interior
    % so nu=1 sufficient to get exact value
else
    e_c = zeros(N/2+1,N/2+1) ;
    for j=1:gamma
        e_c = MG(l-1,e_c,-r_c,gamma,nu_1,nu_2,N/2) ;
    end
end

e = PROLONG(e_c,N/2) ;

for i=1:N+1
    for j=1:N+1
        u(i,j) = u(i,j) - e(i,j) ;
    end
end

u = GS(u,f,nu_2,N) ;
end

