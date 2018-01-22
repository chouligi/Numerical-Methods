function [l_collection,c_collection,res_collection] = efficiency_assignment3_exercise1(data)

% this function checks the residuals for vector of initial guesses 
%with variable length from 1 to 20 and stores the results of each iteration 

l_collection = [];
c_collection =[];
res_collection = [];

for i = 1:20
    
    g = linspace(-4,0,i);
    [l_best,c_best,residue] = Chouliaras_assignment3_exercise1(data, g);
   
    l_collection{i} = l_best;
     c_collection{i} = c_best;
    res_collection(i) = residue;
    
     end
    
end