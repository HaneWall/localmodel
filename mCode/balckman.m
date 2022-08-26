function [win] = balckman(wlen, option)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if option == 'periodic'
    L = wlen + 1;
    win = zeros(L, 1);
    for n = 1:1:L
        win(n) = 0.42 - 0.5*cos(2*pi*n / (L - 1)) + 0.08*cos(4*pi*n / ( L - 1));
    end
else
    L = wlen;
    win = zeros(L, 1);
    for n = 1:1:L
        win(n) = 0.42 - 0.5*cos(2*pi*n / (L - 1)) + 0.08*cos(4*pi*n / ( L - 1));
    end
end
end