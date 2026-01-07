clear; clc;

digits = 6;

dir = "C:\Users\ellil\Documents\WORK\MLFMA\MLFMA\out\build\x64-debug\out\sol\";

sol1 = readmatrix(dir+"rvec_nq7.txt");
sol2 = readmatrix(dir+"rvec_nq7_anterp.txt");
solDir = readmatrix(dir+"rvecDir_nq7.txt");
% sol1 = readmatrix(dir+"curr_nq7.txt");
% sol2 = readmatrix(dir+"currDir_nq7.txt");

nvec = 1:length(sol1);

%% 
close all;
for i=1:2
    sol1Sort = sortrows(sol1,i);
    sol2Sort = sortrows(sol2,i);
    solDirSort = sortrows(solDir,i);
    figure(2*i-1);
    plot(nvec, sol1Sort(:,i), nvec, solDirSort(:,i));
    % semilogy(nvec, abs(sol2Sort(:,i)), nvec, abs(solDirSort(:,i)));

    % relErr0 = abs(sol1Sort(:,i)-sol2Sort(:,i)) ./ abs(sol1Sort(:,i));
    relErr1 = abs(sol1Sort(:,i)-solDirSort(:,i)) ./ abs(solDirSort(:,i));
    relErr2 = abs(sol2Sort(:,i)-solDirSort(:,i)) ./ abs(solDirSort(:,i));

    figure(2*i);
    semilogy(nvec, relErr1, nvec, relErr2);
end

% mean(relErr0)
mean(relErr1)
mean(relErr2)
