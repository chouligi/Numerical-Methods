function FFT = Chouliaras_assignment3_exercise3(v,type,count)
%G.C. Chouliaras
%This function calculates both the discrete Fourier transform and its
%inverse Fourier transform.
%The input is a vector of length 2^n and a string variable which indicates
%whethere the normal ('nor') or the inverse ('inv') transform needs to be calculated.
%The output is a vector of the Fourier transform.
%
%
%
%Example syntax:
%Calculate the Fourier transform of the vector v = linspace(1,8,8):
% fft = Chouliaras_assignment3_exercise3(v,'nor')
%Calculate its inverse Fourier transform:
% fft = Chouliaras_assignment3_exercise3(v,'inv')



%store the length of the vector v
N = length(v);

%check if count exists and if it doesn't give value equal to N
if ~exist('count')
    count = length(v);
end

%chech if the type is normal
if type == 'nor'
    
    omega = exp(-2*pi*1i/N);
    
    %when N has length 1, finish the algorithm and store to fft the latest
    %value of v
    if N==1
        FFT = v;
        return;
        
    else
        %calculate the FFT
        FFT = [Chouliaras_assignment3_exercise3(v(1:2:N),'nor') + omega.^(0:(N/2-1)).*Chouliaras_assignment3_exercise3(v(2:2:N),'nor'),Chouliaras_assignment3_exercise3(v(1:2:N),'nor') -  omega.^(0:(N/2-1)).*Chouliaras_assignment3_exercise3(v(2:2:N),'nor')];
    end
    
else
    %in inverse, omega has positive exponential
    omega = exp(2*pi*1i/N);
    
    %if N cannnot be divided further, find the inverse FFT by dividing with initial
    %length of vector
    if N==1
        FFT = v/count;
        return;
    else
        FFT = [Chouliaras_assignment3_exercise3(v(1:2:N),'inv',count) + omega.^(0:(N/2-1)).*Chouliaras_assignment3_exercise3(v(2:2:N),'inv',count),Chouliaras_assignment3_exercise3(v(1:2:N),'inv',count) -  omega.^(0:(N/2-1)).*Chouliaras_assignment3_exercise3(v(2:2:N),'inv',count)];
        
    end
    
end