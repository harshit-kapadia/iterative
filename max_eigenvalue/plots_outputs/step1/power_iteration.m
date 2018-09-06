clc ;
close all ;
clear all ;

load('power-iteration.mat') ;

semilogy(iteration, residual, '-', 'LineWidth', 1.7) ;

title('Convergence plot (Power iteration on nos6.mtx)' ,'FontSize',15, 'Interpreter','latex')
xlabel('Iteration index, k', 'Interpreter', 'latex', 'FontSize',14)
ylabel('$| \lambda^{(k)} - \lambda^{(k-1)} |$', 'Interpreter', 'latex',...
    'FontSize', 15)