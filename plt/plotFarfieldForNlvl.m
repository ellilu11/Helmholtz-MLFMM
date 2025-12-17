clear; clc;

nth = 47;

maxLvls = 0:3;
nLvls = size(maxLvls,2);

flds = cell(nLvls,1);
thetas = cell(nLvls,1);
phis = cell(nLvls,1);

dir = "C:\Users\ellil\Documents\WORK\MLFMA\MLFMA\out\build\x64-debug\out\ff\";

for i=1:size(maxLvls,2)
    flds{i} = readmatrix(dir+"ff_maxlvl"+string(maxLvls(i))+"_nth"+string(nth)+".txt");
    thetas{i} = readmatrix(dir+"thetas_nth"+string(nth)+".txt");
    phis{i} = readmatrix(dir+"phis_nth"+string(nth)+".txt");
end

% fldDir = readmatrix(dir+"out\ffDir.txt");
% nangles = size(fldDir,1);
% nvec = 1:nangles;

% Select levels and field components to plot
lvl1 = 1; lvl2 = 3;
comp = 2; % E_x = 1, E_y = 2, E_z = 3

ymax = max([max(flds{lvl1}), max(flds{lvl2})]);

%% Reshape fld matrices
fldArr1 = reshape(flds{lvl1}(:,comp), nph(lvl1), nth(lvl1));
fldArr2 = reshape(flds{lvl2}(:,comp), nph(lvl2), nth(lvl2));

thetas1 = thetas{lvl1};
thetas2 = thetas{lvl2};

phis1 = phis{lvl1};
phis2 = phis{lvl2};

%% Plot as function of theta
close all;

figure(1);
for iph = 1:nph(lvl1)
    plot(thetas1, fldArr1(iph,:), thetas2, fldArr2(iph,:));
    xlim([0 pi]);
    ylim([0 ymax]);
    xlabel('th');
    title('ph = '+string(phis1(iph)))
    pause(0.05);
end

%% Plot as function of phi
close all;
figure(2);
for ith = 1:nth(lvl1)
    plot(phis1, fldArr1(:,ith), phis2, fldArr2(:,ith));
    xlim([0 2*pi]);
    ylim([0 ymax]);
    xlabel('ph');
    title('th = '+string(thetas1(ith)))
    pause(0.05);
end
