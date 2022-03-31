% NMF Algorithm to Separate 2 overlapped AIS Signals
clc;
clear all;
close all;

%% Read Dataset and Extract 2 AIVDM message from Dataset and get AIS Signals
filename_dataset = 'D:\0-MSc_Course\MSc Project\AIS_Python\ais_msgs_sample_V1.0.txt';
dataset = importdata(filename_dataset);
index = randperm(size(dataset,1),2);
index = [616,109];
aismsg1=char(dataset(index(1)));
aismsg2=char(dataset(index(2)));

[sigI1,sigQ1,sigIQ1]=ais_gmsk_signal_gen(aismsg1);
[sigI2,sigQ2,sigIQ2]=ais_gmsk_signal_gen(aismsg2);

figure();
subplot(2,1,1)
plot(sigI1)
title('In-phase');xlabel('Samples');ylabel('Amplitude')
set(gca,'XLimSpec', 'Tight','fontname','Times New Roman','fontsize',9)

subplot(2,1,2)
plot(sigQ1)
title('Quadrature');xlabel('Samples');ylabel('Amplitude')
set(gca,'XLimSpec', 'Tight','fontname','Times New Roman','fontsize',9)  

figure();
subplot(2,1,1)
plot(sigI2)
title('In-Phase');xlabel('Samples');ylabel('Amplitude')
set(gca,'XLimSpec', 'Tight','fontname','Times New Roman','fontsize',9)  

subplot(2,1,2)
plot(sigQ2)
title('Quadrature');xlabel('Samples');ylabel('Amplitude')
set(gca,'XLimSpec', 'Tight','fontname','Times New Roman','fontsize',9)  


same_len = max(size(sigI1,2),size(sigI2,2));
delay = randi([1,5000]);
delay = 2109;

s1=[sigI1,zeros(1,20000)];
s2=[zeros(1,delay),sigI2,zeros(1,20000)];
s1 = s1(1:same_len);
s2 = s2(1:same_len);

figure();
subplot(2,1,1)
plot(s1)
title('Signal 1 (In-Phase)');xlabel('Samples');ylabel('Amplitude')
set(gca,'XLimSpec', 'Tight','fontname','Times New Roman','fontsize',9)  

subplot(2,1,2)
plot(s2)
title('Signal 2 (In-Phase)');xlabel('Samples');ylabel('Amplitude')
set(gca,'XLimSpec', 'Tight','fontname','Times New Roman','fontsize',9)  

%% Generation Hybrid channel and mixedsignal with 2x2 mixing matrix
ch = [0.75,0.44;0.30,0.79];
x = ch*[s1;s2];

figure();
subplot(2,1,1)
plot(x(1,:))
title('Antenna 1');xlabel('Samples');ylabel('Amplitude')
set(gca,'XLimSpec', 'Tight','fontname','Times New Roman','fontsize',9)  
hold on; plot(x(1,1:delay))

subplot(2,1,2)
plot(x(2,:))
title('Antenna 2');xlabel('Samples');ylabel('Amplitude')
set(gca,'XLimSpec', 'Tight','fontname','Times New Roman','fontsize',9)  
hold on; plot(x(2,1:delay))

n1=0;
n2=0;
%% Noise
snr = 12; % [dB]
Ps_Pn = 10^(snr/10);
n1=(1/sqrt(Ps_Pn))*randn(1,size(x,2))*sqrt(mean(x(1,:).^2));
n2=(1/sqrt(Ps_Pn))*randn(1,size(x,2))*sqrt(mean(x(2,:).^2));

x_noisy(1,:) = x(1,:) + n1;
x_noisy(2,:) = x(2,:) + n2;


%%
x_noisy(1,:) = x(1,:) + n1;
x_noisy(2,:) = x(2,:) + n2;

W = [0.75,rand();0.3,rand()];
H = rand(2,size(x,2))+5;
V=x_noisy+4;
while 1
    if norm(V-W*H,'fro')<10
        break;
    else
        
    
        wv=W'*V;
        wwh=(W)'*W*H;
        for i=1:2
            for j=1:size(V,2)
                H(i,j) = H(i,j)*(wv(i,j)/wwh(i,j));
%                 if H(i,j)>4
%                     H(i,j) = 4;
%                 end
            end
        end

        vh_n1 = V*H';
        whh_n1 = W*H*(H)';
        for i=1:2
            for j=1:2
                if (i == 1 && j == 1) || (i == 2 && j == 1)
                    continue;
                else
                    W(i,j) = W(i,j)*(vh_n1(i,j)/whh_n1(i,j));
            
                end
%                 W(i,j) = W(i,j)*(vh_n1(i,j)/whh_n1(i,j));
            end
        end
    end
end


%%
B =H;

figure();plot(s2);hold on;
plot(B(2,:)-1.5)
title('Component 1');xlabel('Samples');ylabel('Amplitude')
set(gca,'fontname','Times New Roman','fontsize',9)  
axis tight

figure();plot(s1);hold on;
plot(B(1,:)-1)
title('Component 2');xlabel('Samples');ylabel('Amplitude')
set(gca,'fontname','Times New Roman','fontsize',9)  
axis tight

