function A = getstat(V)
    
    %A = zeros(4,5); 
            A(:,1)    = mean(V,3);
            A(:,2)    = std(V,0,3);
            A(:,3)    = max(V,[],3);
            A(:,4)    = min(V,[],3);
            A(:,5)    = median(V,3);
end
      