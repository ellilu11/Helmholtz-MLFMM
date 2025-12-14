clear; clc;

nths = [30, 31];
nLvls = length(nths);

coeffs = cell(nLvls,1);
thetas = cell(nLvls,1);
phis = cell(nLvls,1);

nth = zeros(nLvls,1);
nph = zeros(nLvls,1);

dir = "C:\Users\ellil\Documents\WORK\MLFMA\MLFMA\out\build\x64-debug\";

for i=1:size(nths,2)
    coeffs{i} = readmatrix(dir+"out\lcoeffs_nth"+string(nths(i))+".txt");
    thetas{i} = readmatrix(dir+"out\thetas_nth"+string(nths(i))+".txt");
    phis{i} = readmatrix(dir+"out\phis_nth"+string(nths(i))+".txt");

    nth(i) = length(thetas{i});
    nph(i) = length(phis{i});
end

lvl1 = 1; lvl2 = 2;
comp = 1;

%% Reshape fld matrices
coeffArr1 = reshape(coeffs{lvl1}(:,comp), nph(lvl1), nth(lvl1));
coeffArr2 = reshape(coeffs{lvl2}(:,comp), nph(lvl2), nth(lvl2));

thetas1 = thetas{lvl1};
thetas2 = thetas{lvl2};

phis1 = phis{lvl1};
phis2 = phis{lvl2};

%% Plot as function of theta
close all;
ymax = max([max(coeffArr1), max(coeffArr2)]);

figure(1);
for iph = 1:nph(lvl1)
    plot(thetas1, coeffArr1(iph,:), thetas2, coeffArr2(iph,:));
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
    plot(phis1, coeffArr1(:,ith), phis2, coeffArr2(:,ith));
    xlim([0 2*pi]);
    ylim([0 ymax]);
    xlabel('ph');
    title('th = '+string(thetas1(ith)))
    pause(0.05);
end

%%
% figure(1);
% for ith = 1:nth
%     plot(1:nph, coeffs1(:,ith,1), 1:nph, coeffs2(:,ith,1))
%     title('th = '+string(pi*(ith-1)/nth))
%     pause(0.2);
% end
% 
% figure(2);
% for iph = 1:nph
%     plot(1:nth, coeffs1(iph,:,1), 1:nth, coeffs2(iph,:,1))
%     title('ph = '+string(2*pi*(iph-1)/nph))
%     pause(0.2);
% end