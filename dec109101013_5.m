% Parameters;
fc = 100; % carrier frequency
fs = 1000; % sampling rate

% Read samples of x5(t) into x
fileID = fopen('5.txt', 'r');
x = fscanf(fileID, '%f');
fclose(fileID);

% An array from t = 0 to t = 32, with sampling rate = 1000Hz
t = 0:0.001:32-0.001;
t = transpose(t);

% Differentiator
dt = 0.001;
% returns [x2-x1, x3-x2, ... ]
dx = diff(x);
dxdt = dx/dt;
% Padding, since dxdt is of length 31999
e(1) = 0;
e(2:32000) = dxdt;

% Squarer
y = e.^2;

% Low-pass filter, cutoff @ 20Hz
z = lowpass(y, 20, fs);

% Square rooter
w = sqrt(z);

% DC Blocker: Since both m(t) and the noise are zero-mean,
% we can remove the DC offset by subtracting the mean of w(t) from w(t)
m = transpose(w - mean(w));

% Recover its magnitude
m = m * sqrt(2);
m = m / (2 * pi * 70);

% dlmwrite('109101013_results.txt', m, '-append');

% Plot m(t)
for i = 1:4
    subplot(4, 1, i);
    plot(t(1+8000*(i-1):8000*i), m(1+8000*(i-1):8000*i));
    title(sprintf('rom t = %d to t = %d',  8 * (i-1), 8 * i));
end
set(gcf, 'position', [10, 10, 1280, 720]);
saveas(gcf, 'm5', 'epsc');


% Recover the bit sequence
A = sin(2*pi*t); % A represents how a sequence of 1's is sent
bin_seq = zeros(32, 1); % sequence with all elements set to zero
% Check the entire signal duration of 32 seconds
for i = 0:31
    pos = 0;
    neg = 0;
    % For each second, check if this part of m(t) 
    % is in-phase with sin(2*pi*t)
    for j = 1:fs
        % For every sample of m(t), 
        % if it has the same sign as the
        % corresponding point in A, then pos += 1;
        eval = m(i*fs+j) / A(i*fs+j);
        if eval > 0
            pos = pos + 1;
        else
            neg = neg + 1;
        end
    end
    % If pos > neg, then this part of m(t) 
    % is "probably" in-phase with sin(2*pi*t)
    % meaning that this bit should be set to 1
    % Otherwise, it should be left unchanged (0)
    if pos > neg
        bin_seq(i+1) = 1;
    end
end

% Join the binary sequence into a string
bitString = strjoin(string(int8(bin_seq)));
% Convert to decimal
dec = bin2dec(bitString);
% Display results
disp(bitString);
disp(dec);