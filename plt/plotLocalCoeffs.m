clear; clc;

nth = 18;
nph = 2*nth;

dir = "C:\Users\ellil\Documents\WORK\MLFMA\MLFMA\out\build\x64-debug\";
coeffs1 = readmatrix(dir+"out\lcoeffs.txt");
coeffs2 = readmatrix(dir+"out\lcoeffs_p6.txt");

coeffs1 = reshape(coeffs1, nph, nth, 6);
coeffs2 = reshape(coeffs2, nph, nth, 6);

%%
figure(1);
for ith = 1:nth
    plot(1:nph, coeffs1(:,ith,1), 1:nph, coeffs2(:,ith,1))
    title('th = '+string(pi*(ith-1)/nth))
    pause(0.2);
end

figure(2);
for iph = 1:nph
    plot(1:nth, coeffs1(iph,:,1), 1:nth, coeffs2(iph,:,1))
    title('ph = '+string(2*pi*(iph-1)/nph))
    pause(0.2);
end