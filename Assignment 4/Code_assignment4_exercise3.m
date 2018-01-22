function [rankings] = Chouliaras_assignment4_exercise3(A,websites)
%G.C. Chouliaras
%This function makes a page ranking of webpages. Given a sparse nxn connection
%matrix A for which A_ij = 1 if page j contains a link to page i and 0
%otherwise, we consider a markov process where we follow an arbitrary link
%with probability p=0.85 and we jump to an arbitrary link with probability
%1-p. This markov process is described by a matrix B_ij = p/c_j A_ij +
%(1-p)/n if c_j is not zero and B_ij = 1/n if c_j is zero. C_j is the
%number of 1's in the jth column of matrix A. The eigenvector corresponding
%to the largest eigenvalue of matrix B is the limiting long term
%distribution of the markov process. This eigenvector contains the page
%ranks of the pages. The largest eigenvalue of B is found using the power
%method.
%The inputs of the function are the nxn sparse connection matrix A and
%optionally a list of associated websites.
%The output of the function is a vector (rankings) which contain the page
%rankings. Moreover, the function produces a bar chart with the page
%rankings. If the list of websites is given as optional input, the function
%also produces a list with the top 5 websites and their page ranks.
%
%Example syntax:
%To find the page ranks of the sparse connection matrix sp100 from
%wwwdata.mat with the associated list of websites url100 one should type:
%[rankings] = Chouliaras_assignment4_exercise3clean(sp100,url100)

n = length(A);

p = 0.85;

%initialize y0
y0 = zeros(n,1);

%initialize y1
y1 = zeros(n,1);



% create a vector c of size n which contains the number of links for each
% page

c = sum(A);

%create auxiliary diagonal sparse matrix to store 1/c if c ~= 0 otherwise 0

%spdiags creates a nxn sparse matrix from the columns of 1./c' and places
%them along the diagonals specified by 0.
H = spdiags(1./c',0,n,n);

%since there are infinite values in the matrix we set them to zero

H(isinf(H)) =0;


%implement the power method
%difference = 1 as initial value
difference = 1;
while difference > eps
    
    
    y0 = y1;
    
    x1 = p*A*H*y0 + ((1-p)/n)* ones(n,1); 
    
    %normalize the vector by dividing with the sum of all elements to be 1
    y1 = x1/sum(x1);
    difference = norm(y0-y1);
    
end

%return the eigenvector as output
rankings = y1;

%plot a barchart with the page rankings
bar(rankings,'k');
title('Barchart with the website rankings');
xlabel('Websites');
ylabel('Page Rank');
set(gca,'fontsize',14);

%check if the variable websites is given as input
if ~exist('websites','var')
    
else
    %some url's contain row vectors so we turn them to column vectors
    if isrow(websites)
        websites = websites';
    end
    
    if n>5
        %we add one more row to put the headers
        top5 = cell(6,2);
    else
        top5 = cell(n+1,2);
    end
    %put the headers
    top5{1,1} = 'Website';
    top5{1,2} = 'Page Rank';
    
    %we start from 2 because the 1st row contains the headers
    for i = 2:length(top5)
        [sranks,index] = sort(rankings,'descend');
        %store the names of the websites
        top5{i,1} = websites{index(i-1),1};
        %store the page ranks
        top5{i,2} = sranks(i-1);
    end
    %present the top 5 websites
    top5
end

end
