function e_l = MG(l,u,f,gamma,nu_1,nu_2,N)
% Multigrid Algorithm
u = GS(u,f,nu_1,N) ; % smoothing and overwriting u_l
r = residual(u,f,N) ;
r_c = RESTR(r,N/2) ;

if l==1
    e_c = zeros(N/2+1,N/2+1) ;
    e_c = GS(e_c,-r_c,1,N/2) ;
else
    e_c = zeros(N/2+1,N/2+1) ;
    for j=1:gamma
        e_c = MG(l-1,e_c,-r_c,gamma,nu_1,nu_2,N/2) ;
    end
end

e = PROLONG(e_c,N/2) ;
u = u - e ;
u = GS(u,f,nu_2,N) ;
end

