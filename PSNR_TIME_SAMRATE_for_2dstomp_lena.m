%% Fig about: PSNR/SSIM/MSE & running time for diff samrate 
% for 2d reconstruct algorithms
% composed by Rinabell
% version1.0 @18-04-24
% version2.0 @18-05-03 to change 2DOMP & add 2DROMP
% version3.0 @18-06-21 find error in pramenter psn and err
% and add reconstruction probability 
% version3.5 to fix changes in k
% version4.0 @191217 to add SSIM/MSE

%% ini
clc;clear;close all;
addpath('otheralg'); load('lena.mat'); 
% kind = [ '2DSTOMP';'2DROMP';'2DOMP'; 'SL0_2D'];
kind = 4;

img = lena;
% img = imread('uu.png');
% img = im2double(img);
[~,n] = size(img);
THR = 1e2;
NOL = 100;

samrate = 0.1:0.1:0.9;
rate = length(samrate);  %from 0.1 to 1;
num = 0;
CHART_TIME_SAMRATE = zeros(rate,kind);
CHART_PSNR_SAMRATE = zeros(rate,kind);
CHART_SSIM_SAMRATE = zeros(rate,kind);
CHART_MSE_SAMRATE = zeros(rate,kind);

%% script
for samrate = 0.1:0.1:0.9
    m = floor(n * samrate);
    k = floor(sqrt(m));
    Psi = (2/n)^0.5 * cos(pi * ((0:(n-1))+0.5)' * (0:(n-1)) / n);
    Psi(:,1) = Psi(:,1) / (2^0.5);
    num = num + 1;
    % NOTICE: here to put them to 0;
    err = zeros(kind,1); psn = zeros(kind,1); 
    pro = zeros(kind,1); ssi = zeros(kind,1); 
    ms = zeros(kind,1); 
    
    
    for time = 1:NOL
        fprintf('samrate = %1.2f, time = %d;\n', samrate, time);
        Phi = randn(m,n) / (m^0.5);
        A = Phi * Psi;
        A_t = A';
        C = A_t * A;
        N = zeros(n);
        for i=1:n
            for j=1:n
                N(i,j) = sqrt(C(i,i)*C(j,j));
            end
        end
        y = Phi * double(img) * Phi';
        
        % 2dstomp        
        tic
        out1 = stomp2d(y, A, A_t, C, N, k);
        toc
        out1 = Psi * out1 *Psi';
%         figure;imshow(out1);
        err(1) = err(1) + toc;
        psn(1) = psn(1) + psnr(out1,img);
        ssi(1) = ssi(1) + ssim(out1,img);
        ms(1) = ms(1) + mse(out1,img);
       
        % 2DROMP
        tic
        out2 = romp2d_v2_5(y, A, A_t, C, N, k, 8);
        toc
        out2 = Psi * out2 *Psi';
        err(2) = err(2) + toc;
        psn(2) = psn(2) + psnr(out2,img);
        ssi(2) = ssi(2) + ssim(out2,img);
        ms(2) = ms(2) + mse(out2,img);
        
        % 2DOMP        
        tic
        out3 = fomp22(y, A, A_t, C, N, k);
        toc
        out3 = Psi * out3 *Psi';
%         figure;imshow(out3);
        err(3) = err(3) + toc;
        psn(3) = psn(3) + psnr(out3,img);
        ssi(3) = ssi(3) + ssim(out3,img);
        ms(3) = ms(3) + mse(out3,img);
        
        %2d sl0
%         h1 = randn(m,n); h2 = randn(m,n);
        y4 = Phi * im2double(img) * Phi';
        tic
        out4 = SL0_2D(y4,Phi,Phi,n);
        toc
%         out4 = im2uint8(out4);
        err(4) = err(4) + toc;
        psn(4) = psn(4) + psnr(out4,img);
        ssi(4) = ssi(4) + ssim(out4,img);
        ms(4) = ms(4) + mse(out4,img);
                
    end
    
    for temp = 1:kind
        CHART_TIME_SAMRATE(num,temp) = err(temp) / NOL;
        CHART_PSNR_SAMRATE(num,temp) = psn(temp) / NOL;
        CHART_SSIM_SAMRATE(num,temp) = ssi(temp) / NOL;
        CHART_MSE_SAMRATE(num,temp) = ms(temp) / NOL;
    end
end

%% print
save OUTPUT_ALL_TIME_SAMRATE_100

s = [' -rs';'--ko';'--kd';'--kx'];
ff1 = figure(1);
samrate = 0.1:0.1:0.9;
rate = length(samrate);
for temp = 1:kind
    plot(samrate,CHART_PSNR_SAMRATE(1:rate,temp),s(temp,:),'LineSmoothing','on');
    hold on;
end
hold off; legend('2D-STOMP','2D-SP','2D-OMP','2D-SL_0');
xlabel('Sampling Rate','Fontname', 'Times New Roman');
ylabel('PSNR Values (db)','Fontname', 'Times New Roman');

ff2 = figure(2);
for temp = 1:kind
    plot(samrate,CHART_TIME_SAMRATE(1:rate,temp),s(temp,:),'LineSmoothing','on');
    set(gca,'YScale','log')
    hold on;
end
hold off; legend('2D-STOMP','2D-SP','2D-OMP','2D-SL_0');
xlabel('Sampling Rate','Fontname', 'Times New Roman');
ylabel('Running Time (s)','Fontname', 'Times New Roman');

ff4 = figure(3);
for temp = 1:kind
    plot(samrate,CHART_SSIM_SAMRATE(1:rate,temp),s(temp,:),'LineSmoothing','on');
    hold on;
end
hold off; legend('2D-STOMP','2D-SP','2D-OMP','2D-SL_0');
xlabel('Sampling Rate','Fontname', 'Times New Roman');
ylabel('SSIM Values (1)','Fontname', 'Times New Roman');

ff5 = figure(4);
for temp = 1:kind
    plot(samrate,CHART_MSE_SAMRATE(1:rate,temp),s(temp,:),'LineSmoothing','on');
    hold on;
end
hold off; legend('2D-STOMP','2D-SP','2D-OMP','2D-SL_0');
xlabel('Sampling Rate','Fontname', 'Times New Roman');
ylabel('MSE Values (1)','Fontname', 'Times New Roman');

pack;



