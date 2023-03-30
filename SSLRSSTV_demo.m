%%  Name: SSLRSSTV
%   [1]F. Yang, X. Chen, and L. Chai, ¡°Hyperspectral image
%    destriping and denoising using stripe and spectral low-
%   rank matrix recovery and global spatial-spectral total
%   variation,¡± Remote Sensing, vol. 13, no. 4, p. 827, 2021.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all;
load Pavia_80.mat;
load oriData3_noise.mat;
[M N p] = size(oriData3_noise);
[noipsnr, noissim, noimsam] = MSIQA3(OriData3*255, oriData3_noise*255);
tic;
[opts] = ParSetSSTV;
[SSLRSSTVoutput_image B S J W] =SSLR_SSTV(oriData3_noise,opts);
for i=1:p
SSLRSSTVoutput_image(:,:,i)=SSLRSSTVoutput_image(:,:,i)+mean(mean(S(:,:,i)));
end
time(1) = toc;
[SSLRSSTVpsnr, SSLRSSTVssim, SSLRSSTVmsam] = MSIQA(OriData3*255, SSLRSSTVoutput_image*255);
