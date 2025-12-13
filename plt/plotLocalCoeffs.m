clear; clc;

nth = 18;
nph = 2*nth;

dir = "C:\Users\ellil\Documents\WORK\MLFMA\MLFMA\out\build\x64-debug\";
coeffs1 = readmatrix(dir+"out\lcoeffs_nolut.txt");
coeffs2 = readmatrix(dir+"out\lcoeffs_lut.txt");

coeffs1 = reshape(coeffs1, nth, nph, 6);
coeffs2 = reshape(coeffs2, nth, nph, 6);

%%
figure(1);
for ith = 1:nth
    plot(1:nph, coeffs1(ith,:,1), 1:nph, coeffs2(ith,:,1))
    title('th = '+string(pi*(ith-1)/nth))
    pause(0.1);
end

figure(2);
for iph = 1:nph
    plot(1:nth, coeffs1(:,iph,1), 1:nth, coeffs2(:,iph,1))
    title('ph = '+string(2*pi*(iph-1)/nph))
    pause(0.1);
end