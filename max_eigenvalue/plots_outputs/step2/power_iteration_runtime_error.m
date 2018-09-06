clc ;
close all ;
clear all ;

runtime = [3.726800e-02 1.829910e-01 3.840750e-01 1.849739] ;
error = [6.498889e+01 2.948186e-01 2.797706e-03 4.365575e-11] ;

semilogy(runtime, error, 'o-', 'LineWidth', 1.7) ;

labels = {'k=100', 'k=500', 'k=1000', 'k=5000'} ;
text(runtime, error, labels, 'horizontal','left', 'vertical','bottom', 'FontSize', 12)

title('Error vs Execution time (Power iteration on s3rmt3m3.mtx)' ,'FontSize',15, 'Interpreter','latex')
xlabel('Runtime (sec)', 'Interpreter', 'latex', 'FontSize',14)
ylabel('error ($ \lambda^{(true)} - \lambda^{(computed)} $)', 'Interpreter', 'latex',...
    'FontSize', 15)