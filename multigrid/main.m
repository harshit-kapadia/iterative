% Main script of the Multigrid Project
% Author: Harshit J. Kapadia

clc ;
clear all ;
close all ;

n = 4 ; % n is as defined in Project 2 sheet, i.e. N=2^n
% n = 7 ;
N = 2^n ; % instead of multiplying 2 n times in a loop ^ is used here, as
% pow() function already present in math.h for C
l = n-1 ; % using as many grid levels as possible, i.e. in such a way so that
% last level has only one interior node so that GS becomes exact in
% one step
gamma = 2 ; % for W cycles
% gamma = 1 ; % for V cycles
nu_1 = 1 ; nu_2 = 1 ;
% nu_1 = 2 ; nu_2 = 1 ;

h = 1.0/N ;
m_max = 50 ; % maximum multigrid iterations

u = zeros(N+1,N+1) ; % as index varies from 1 to N+1 instead of 0 to N
u_exact = zeros(N+1,N+1) ;
f = zeros(N+1,N+1) ;
r = zeros(N+1,N+1) ; % for residual
r_0 = zeros(N+1,N+1) ; % for initial residual

for i=1:N+1
    for j=1:N+1
        u_exact(i,j) = sin(2*pi*(i-1)*h) * sin(2*pi*(j-1)*h) ;
        f(i,j) = 8*pi*pi * sin(2*pi*(i-1)*h) * sin(2*pi*(j-1)*h) ;
    end
end
r_0 = residual(u,f,N) ;
r_0_norm = norm_inf(r_0,N) ;

r_norm = zeros(m_max,1) ;
for m=1:m_max
    u = MG(l,u,f,gamma,nu_1,nu_2,N) ;
    r = residual(u,f,N) ;
    r_norm(m,1) = norm_inf(r,N) ;
end

ratio = zeros(m_max,1) ;
for m=1:m_max
    ratio(m,1) = r_norm(m,1) / r_0_norm ;
end

figure
semilogy([1:1:m_max],ratio, '-o') 


