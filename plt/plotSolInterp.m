clear; clc;

maxLvls = 0:2;
nLvls = size(maxLvls,2);

flds = cell(nLvls,1);
thetas = cell(nLvls,1);
phis = cell(nLvls,1);

nth = zeros(nLvls,1);
nph = zeros(nLvls,1);

dir = "C:\Users\ellil\Documents\WORK\MLFMA\MLFMA\out\build\x64-debug\";

for i=1:size(maxLvls,2)
    flds{i} = readmatrix(dir+"out\ff_maxlvl"+string(maxLvls(i))+".txt");
    thetas{i} = readmatrix(dir+"out\thetas_lvl"+string(maxLvls(1))+".txt");
    phis{i} = readmatrix(dir+"out\phis_lvl"+string(maxLvls(1))+".txt");

    nth(i) = length(thetas{i});
    nph(i) = length(phis{i});
end

% Select levels and field components to plot
lvl1 = 1; lvl2 = 2;
comp = 1; % E_x = 1, E_y = 2, E_z = 3

%%
% fldDir = readmatrix(dir+"out\ffDir.txt");
% nangles = size(fldDir,1);
% nvec = 1:nangles;

% close all;
% lvlpp = 1;
% 
% figure(1);
% for i=1:3
%     plot(nvec, flds{lvlpp}(:,i), nvec, fldDir(:,i));
% end
% 
% relErrFld = abs(flds{lvlpp}(:,1) - fldDir(:,1))./abs(fldDir(:,1));
% figure(2);
% semilogy(nvec, relErrFld)
%% Reshape fld matrices
fldArr1 = reshape(flds{lvl1}, nth(lvl1), nph(lvl1), 3);
fldArr2 = reshape(flds{lvl2}, nth(lvl2), nph(lvl2), 3);

%% Interpolate sampled functions on to regular grid
% [th_1, ph_1] = meshgrid(thetas{lvl1}, phis{lvl1});
[ph_1, th_1] = meshgrid(phis{lvl1}, thetas{lvl1});

% [th_2, ph_2] = meshgrid(thetas{lvl2}, phis{lvl2});
[ph_2, th_2] = meshgrid(phis{lvl2}, thetas{lvl2});

nthQ = 50;
nphQ = 100;
thetasQ = linspace(0,pi,nthQ);
phisQ = linspace(0,2*pi,nphQ);
% [thQ, phQ] = meshgrid(thetasQ, phisQ);
[phQ, thQ] = meshgrid(phisQ, thetasQ);

fldInterp1 = interp2(th_1, ph_1, fldArr1(:,:,comp), thQ, phQ, 'linear');
fldInterp2 = interp2(th_2, ph_2, fldArr2(:,:,comp), thQ, phQ, 'linear');

ymax = max([max(fldInterp1), max(fldInterp2)]);

%% Plot as function of theta
close all;
figure(1);
for iph = 1:nphQ
    plot(thetasQ, fldInterp1(:,iph), thetasQ, fldInterp2(:,iph));
    xlim([0 pi]);
    ylim([0 ymax]);
    xlabel('th');
    title('ph = '+string((iph-1)/nphQ*2*pi))
    pause(0.5);
end

%% Plot as function of phi
close all;
figure(2);
for ith = 1:nthQ
    plot(phisQ, fldInterp1(ith,:), phisQ, fldInterp2(ith,:));
    xlim([0 2*pi]);
    ylim([0 ymax]);
    xlabel('ph');
    title('th = '+string((ith-1)/nthQ*pi))
    pause(0.5);
end


