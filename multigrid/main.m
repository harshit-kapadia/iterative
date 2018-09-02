%% Main script of the Multigrid Project
% Author: Harshit J. Kapadia

clc ;
clear all ;
close all ;

%% Selection of the CASE

% % CASE 1
% nu_1 = 1 ; nu_2 = 1 ; 

% CASE 2
nu_1 = 2 ; nu_2 = 1 ; 

%% Assigning other parameters

m_max = 60 ; % maximum multigrid iterations

gamma = 2 ; % for W cycles
% gamma = 1 ; % for V cycles

n_values = [4 7] ;

%% Loop over all required grids

for n=n_values % n is as defined in Project 2 sheet, i.e. N=2^n
    
    N = 2^n ; % instead of multiplying 2 n times in a loop ^ is used here, 
    % as pow() function already present in math.h for C
    
    l = n-1 ; % using as many grid levels as possible, i.e. in such a way 
    % so that last level has only one interior node so that GS becomes 
    % exact in one step. For n=4; l=3,2,1,0 are the levels. 
    
    h = 1.0/N ;
    
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
    
    semilogy([1:1:m_max],ratio, '-o')
    hold on
    
end

title(['Convergence plot with $\gamma$, $\nu_1$, $\nu_2$ as ',...
    num2str(gamma), ', ', num2str(nu_1), ', ', num2str(nu_2),...
    ' (semilog)'],'FontSize',15, 'Interpreter','latex')
leg = legend('N=16','N=128') ;
leg.FontSize = 12 ;
xlabel('No. of multigrid iteration', 'Interpreter', 'latex', 'FontSize',13)
ylabel('$||r^{(m)}||_{inf} / ||r^{(0)}||_{inf}$', 'Interpreter', 'latex',...
    'FontSize', 14)

