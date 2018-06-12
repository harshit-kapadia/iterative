y1 = log(cg.e_normA) ;
y2 = log(cg.r_norm2) ;
x = cg.iteration ;

figure
plot(x, y1, x, y2)
xlabel('# Iteration'), ylabel('Norm value')
legend('Error(A-norm)','Residual(2-norm)')
title('Conjugate Gradient: Error/ Residual Plot (semi-log scale)')