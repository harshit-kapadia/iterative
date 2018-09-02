function result = norm_inf(matrix,N)
% To compute infinity norm as defined in Project 2 sheet
result = 0.0 ;
for i=1:N+1
    for j=1:N+1
        if( abs(matrix(i,j)) > result )
            result = abs(matrix(i,j)) ;
        end
    end
end
end

