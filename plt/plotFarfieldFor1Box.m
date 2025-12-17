nth1 = 26;
nth2 = 39;
nph = 13;

dir = "C:\Users\ellil\Documents\WORK\MLFMA\MLFMA\out\build\x64-debug\";
obss = readmatrix(strcat(dir,"config\obss.txt"));

fld1 = readmatrix(dir+"out\ff_nl1_nth"+string(nth1)+"_nph"+string(nph)+".txt");
fld2 = readmatrix(dir+"out\ff_nl1_nth"+string(nth2)+"_nph"+string(nph)+".txt");

% fldAnl = readmatrix(strcat(dir,"out\ffDir.txt"));

obsvec = 1:size(obss,1);

close all;
%%
figure(1);
plot(obsvec, fld1(:,1), obsvec, fld2(:,1));
% plot(thvec, fld1(1,:), thvec, fld2(1,:))

% hold on;
% figure(2);
% relErr = abs(phiSort-phiAnlSort)./abs(phiAnlSort);
% semilogy(nvec, relErr, '-o');
