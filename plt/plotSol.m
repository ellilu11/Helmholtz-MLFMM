nleaves = [1, 8, 64];

dir = "C:\Users\ellil\Documents\WORK\MLFMA\MLFMA\out\build\x64-debug\";
obss = readmatrix(dir+"config\obss.txt");

fld1 = readmatrix(dir+"out\ff_nl"+string(nleaves(1))+".txt");
fld2 = readmatrix(dir+"out\ff_nl"+string(nleaves(2))+".txt");
fld3 = readmatrix(dir+"out\ff_nl"+string(nleaves(3))+".txt");

fldDir = readmatrix(dir+"out\ffDir.txt");

obsvec = 1:size(obss,1);

%%
close all;
for i=1:3
    figure(i);
    hold on;
    % plot(obsvec, fld1(:,i));
    plot(obsvec, fld2(:,i));
    plot(obsvec, fld3(:,i));
    hold off;
end

figure(4)
plot(obsvec, -fldDir(:,1));

