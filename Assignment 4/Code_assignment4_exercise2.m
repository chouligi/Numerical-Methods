function [digits] = Chouliaras_assignment4_exercise2(thesignal,rate)
%G.C. Chouliaras
%This function distills a phone number from a sound signal. Given the sound
%signal and its sampling rate, the function identifies each digit in the
%signal and finds the low and high frequency associated with each digit.
%Then, based on the Dual Tone Multi-Frequency system it recognizes the
%digits using the low and high frequencies.
%The inputs of the function are the vector containing the sound signal (thesignal) and the
%sampling rate (rate).
%The output is a row vector (digits) containing the digits of the phone
%number.
%
%Example syntax:
%In order to find the digits of the signal "number1" from signal.mat with
%sampling rate 4096, one should type:
%digits = Chouliaras_assignment4_exercise2clean(number1,4096)

%store the length of the vector

l = length(thesignal);

%put a threshold over which the signal is not noise
%we choose 0.6 since after this there is not noise, found by looking at
%the plots of the signals
thres = 0.6;
%counter denotes how many numbers are present in the signal
counter = 0;
index = 1;



%to determine the start of a number we check the absolute value of the
%signal to be larger than the threshold.
%to determine the end of a number we check when 50 subsequent points have
%maximum less than the threshold

while index < (l - 50)
    index = index + 1;
    
    
    if  (abs(thesignal(index))) >= thres
        counter = counter + 1;
        %store the starting points 
        thestart(counter) = index;
        while max(abs(thesignal(index: index + 50))) >= thres
            index = index + 1;
        end
        %store the ending points
        theend(counter) = index;
    end
    
    
end


%store the number of digits
numberdigits = counter;


%put the start and end points in a matrix of size counter x 2
cleansignal = [thestart; theend]';


%allocate space to store the low and high frequencies
low_f = zeros(counter,1);
high_f = zeros(counter,1);

for i = 1:counter
    %the digits
    d = thesignal(cleansignal(i,1): cleansignal(i,2));
    
    %length of digit
    N = length(d);
    
    
    %magnitude of FFT for digit d
    m = abs(fft(d))';
    
  
    %find the peaks of m, set the minimum peak height to 30
    
    [~,peaks] = findpeaks(m,'MinPeakHeight',30);
    
    %we use only the first 2 peaks since the next 2 peaks contain the same
    %information
    
    low_f(i) = peaks(1)*rate/N;
    high_f(i) = peaks(2)*rate/N;
    
    
end

%store the frequencies of each digit in a matrix from 0 to 9 and also * and
%#
numbers = [941,1336; 697,1209; 697,1336; 697, 1477; 770,1209 ; 770,1336; 770,1477; 852,1209; 852,1336; 852,1477; 941,1209; 941,1477];


%store the low and high frequencies into a (counter x 2) matrix
frequencies = [low_f,high_f];

for i = 1:counter
    
    %replicate the matrix in order to take the difference with the numbers
    fr = repmat(frequencies(i,:),12,1);
    
    %take the differences with the numbers
    diff = numbers - fr;
    
    %we find the row for which the 2-norm of the differences is the minimum, the row
    %corresponds to the number
    [~, row] = min(diff(:,1).^2 + diff(:,2).^2);
    
    %we stored the numbers from 0 to 9 so we take -1 to find the
    %appropriate digit
    if row == 11
        digits(i) = '*';
    elseif row == 12
        digits(i) = '#';
    else
        digits(i) = row - 1;
    end
    
    
end
end
