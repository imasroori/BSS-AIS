% FastICA Algorithm to Separate 2 overlapped AIS Signals
clc;
clear all;
close all;

%% Read Dataset and Extract 2 AIVDM message from Dataset and get AIS Signals
filename_dataset = 'D:\0-MSc_Course\MSc Project\AIS_Python\ais_msgs_sample_V1.0.txt';
dataset = importdata(filename_dataset);
index = randperm(size(dataset,1),2);
%index = [747,736];
aismsg1=char(dataset(index(1)));
aismsg2=char(dataset(index(2)));

[sigI1,sigQ1,sigIQ1,stored_msgb1]=ais_gmsk_signal_gen(aismsg1);
[sigI2,sigQ2,sigIQ2,stored_msgb2]=ais_gmsk_signal_gen(aismsg2);

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
%delay=2427;

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

subplot(2,1,2)
plot(x(2,:))
title('Antenna 2');xlabel('Samples');ylabel('Amplitude')
set(gca,'XLimSpec', 'Tight','fontname','Times New Roman','fontsize',9)  

%% Noise
n1 = 0.1*randn(1,size(x,2));
n2 = 0.1*randn(1,size(x,2));
x(1,:) = x(1,:) + n1;
x(2,:) = x(2,:)+ n2;

%% centering data and whitening
m1 = mean(x(1,:));
m2 = mean(x(2,:));
x(1,:) = x(1,:)-m1;
x(2,:) = x(2,:)-m2;

[EigenVectors,eigenValues] = eig(x*transpose(x));
D_12 = eigenValues^(-1/2);
E_T=transpose(EigenVectors);
z = D_12*E_T*x;

z(1,:) = z(1,:)/std(z(1,:));
z(2,:) = z(2,:)/std(z(2,:));

%% Extract components
C = size(z,1); % number of components
M = size(z,2); % number of samples

w = zeros(C,C);
for p=1:C
    stepp=0;
    w_p=rand(2,1);
    w_p_old = w_p;
    
    %w_p_new = (1/M)*(z*(tanh(w_p_old'*z))')-(1/M)*(1-(tanh(w_p_old'*z)).^2)*ones(M,1)*w_p_old;
    w_p_new = z*(tanh(w_p_old'*z))'-mean(1-(tanh(w_p_old'*z)).^2)*w_p_old;
    summ = [0;0];
    for j=1:p-1
        summ = summ + (w_p_new'*w(:,j))*w(:,j);
    end
    w_p_new = w_p_new - summ;

    w_p_new = w_p_new/norm(w_p_new);
    
    while norm(w_p_new-w_p_old) > 1e-12
       disp(w_p_old)
       disp(w_p_new)
       
       stepp = stepp + 1;
       w_p_old = w_p_new;
       %w_p_new = (1/M)*(z*(tanh(w_p_old'*z))')-(1/M)*(1-(tanh(w_p_old'*z)).^2)*ones(M,1)*w_p_old;
       w_p_new = z*(tanh(w_p_old'*z))'-mean(1-(tanh(w_p_old'*z)).^2)*w_p_old;
       summ = [0;0];
       for j=1:p-1
           summ = summ + (w_p_new'*w(:,j))*w(:,j);
       end
       w_p_new = w_p_new - summ;
       
       w_p_new = w_p_new/norm(w_p_new);

    end
    disp(stepp)
    w(:,p)=w_p_new;
end

%%
B =w*z;

figure();
plot(-s1);hold on;
plot(B(2,:))
title('Component 1');xlabel('Samples');ylabel('Amplitude')
set(gca,'fontname','Times New Roman','fontsize',9)  
axis tight

figure();plot(-s2);hold on;
plot(B(1,:))
title('Component 2');xlabel('Samples');ylabel('Amplitude')
set(gca,'fontname','Times New Roman','fontsize',9)  
axis tight

