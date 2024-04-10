function meanX = nanOutMean3(X)

    sizex = size(X,2); % По условию задачи
    
    meanX = zeros(1,sizex);
    
    meanX = sum(X,'omitnan') ./ (sizex - sum(isnan(X)) ) ;

end

