clc ;
close all ;
clear all ;

runtime = [1.209100e-02 2.097200e-02 3.272000e-02 4.664700e-02] ;
error = [1.568834e+01 7.686811e-02 1.520357e-04 1.666558e-08] ;

semilogy(runtime, error, 'o-', 'LineWidth', 1.7) ;

labels = {'m=30', 'm=50', 'm=75', 'm=100'} ;
text(runtime, error, labels, 'horizontal','left', 'vertical','bottom', 'FontSize', 12)

title('Error vs Execution time (Lanczos method on s3rmt3m3.mtx)' ,'FontSize',15, 'Interpreter','latex')
xlabel('Runtime (sec)', 'Interpreter', 'latex', 'FontSize',14)
ylabel('error ($ \lambda^{(true)} - \lambda^{(computed)} $)', 'Interpreter', 'latex',...
    'FontSize', 15)