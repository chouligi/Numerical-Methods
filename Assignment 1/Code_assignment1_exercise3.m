function Chouliaras_assignment1_exercise3(lmin,lmax,N)

% G.C. Chouliaras
%This function produces a picture of the "attractor" of the logistic map.
%The inputs are the minimum and maximum values of lambda lmin, lmax and
%optionally the number of values of lamda between lmin and lmax. The
%default value for N is 1000.
%Appropriate values are 0<= lmin <= lmax <= 4, in any other case the
%function returns an error
%The function has not output, it just produces the picture of the
%attractor.
%
% Example:
%Chouliaras_assignment1_exercise3(0,4,10000)
%

%%
% make input N optional, default value N = 1000
if ~exist('N')
    N=1000;
end

%determine n
n = 1000;

%check that input is valid
if (lmin >= 0 && lmin <= lmax  && lmax <= 4)
    %determine the step for the loop so as to have N values between lmin
    %and lmax
    step = (lmax - lmin)/(N + 1); % without +1 it gives N-1 values
      
        l = (lmin:step:lmax);
        itemsplot = 50;
        x = zeros(itemsplot);
        x(1) = 0.01;
        store = ones(N+2,n);% create array to store the values ,+2 to have the same size as l sequence
        store(:,1) = store(:,1).*x(1); %fill with x0 values
        
      
        plotarray = zeros(N+2, itemsplot);%rows = number of l values, columns = number of x values that will be plotted 
       
    
        counter = 1; %initialize counter to loop through columns of plotarray
        for index = 1:n 
          
            %store the values into store
          store(:,index + 1) = l'.*(store(:,index).*(ones(N+2,1)-store(:,index)));
            
            %pass to plotarray only the values that will be plotted
            if index > n-itemsplot
                plotarray(:,counter) = store(:,index + 1); 
                counter = counter + 1;
            end
        end
        
        %plot
        plot(l,plotarray,'.','MarkerSize',3,'MarkerEdgeColor','k')
        title('Bifurcation Diagram of the Logistic Map','fontsize',13)
        xlabel('Bifurcation parameter lambda','fontsize',13)
        ylabel('x','fontsize',13)
        set(gca,'fontsize',10)
      

else
   
    %in case that input is not appropriate return error
    error('Check input: 0 <= lmin <= lmax <=4')
    
end


