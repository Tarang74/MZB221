clear; close all; clc;

f = @(t) t.^2;
T = 1;
T2 = 2*T;

t = linspace(-3, 3, 2001); t(end) = [];
n = (1:25).';

% Fourier series
fF = 1/3 + sum(1./(n.^2 * pi^2) .* cos(2 * pi / T * n * t) - 1./(n * pi) .* sin(2 * pi / T * n * t), 1);
% Fourier cosine series
fC = 1/3 + sum(4 * (-1).^n ./ (n.^2 * pi^2) .* cos(2 * pi / T2 * n * t), 1);
% Fourier sine series
fS = sum(2 * ((-1).^n .* (2 - n.^2 * pi^2) - 2) ./ (n.^3 * pi^3) .* sin(2 * pi / T2 * n * t), 1);

ttrunc = t(0 <= t & t < T);
ftrunc = f(ttrunc);

figure
plot(ttrunc, ftrunc)
hold on
plot(t, fF)
plot(t, fC)
plot(t, fS)

% Export data
writetable(array2table([ttrunc.' ftrunc.'], 'VariableNames', {'t', 'f'}), 'half_range_trunc.csv') 
writetable(array2table([t.' fF.', fC.', fS.'], 'VariableNames', {'t', 'fF', 'fC', 'fS'}), 'half_range.csv') 