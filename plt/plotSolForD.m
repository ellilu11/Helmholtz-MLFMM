clear; clc;

dir = "C:\Users\ellil\Documents\WORK\MLFMA\MLFMA\out\build\x64-debug\out\sol\";

solDir = readmatrix(dir+"solDir.txt");
nvec = 1:length(solDir);

%% 
digits = 2*(2:6);
ndigits = length(digits);

meanRelErr = zeros(ndigits,2);

for id = 1:ndigits
    sol = readmatrix(dir+"sol_d" + digits(id)+".txt");
    
    for comp = 1:2
        meanRelErr(id,comp) = mean(abs(sol(:,comp)-solDir(:,comp))./abs(solDir(:,comp)));
    end
end

close all;
figure(1);
semilogy(digits, meanRelErr)