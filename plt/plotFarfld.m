clear; clc;

maxLvls = 0:1;
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

%% Interpolate sampled functions on to regular grid
lvl1 = 1; lvl2 = 2;
comp = 1; % E_x = 1, E_y = 2, E_z = 3

[th_1, ph_1] = meshgrid(thetas{lvl1}, phis{lvl1});
fldArr1 = reshape(flds{lvl1}, nph(lvl1), nth(lvl1), 3);

[th_2, ph_2] = meshgrid(thetas{lvl2}, phis{lvl2});
fldArr2 = reshape(flds{lvl2}, nph(lvl2), nth(lvl2), 3);

nthQ = 50;
nphQ = 100;
thetasQ = linspace(0,pi,nthQ);
phisQ = linspace(0,2*pi,nphQ);
[thQ, phQ] = meshgrid(thetasQ, phisQ);

fldInterp1 = interp2(th_1, ph_1, fldArr1(:,:,comp), thQ, phQ, 'linear');
fldInterp2 = interp2(th_2, ph_2, fldArr2(:,:,comp), thQ, phQ, 'linear');

%% Plot 
% as function of theta
close all;
figure(1);
for iph = 1:nphQ
    plot(thetasQ, fldInterp1(iph,:), thetasQ, fldInterp2(iph,:));
    xlim([0 pi]);
    % ylim([0 1.0]);
    xlabel('th');
    title('ph = '+string((iph-1)/nphQ*2*pi))
    pause(0.05);
end

% as function of phi
close all;
figure(2);
for ith = 1:nthQ
    plot(phisQ, fldInterp1(:,ith), phisQ, fldInterp2(:,ith));
    xlim([0 2*pi]);
    % ylim([0 1.0]);
    xlabel('ph');
    title('th = '+string((ith-1)/nthQ*pi))
    pause(0.05);
end


