% AIS-ICA Algorithm to Separate 2 overlapped AIS Signals
clc;
clear all;
close all;
%% Read Dataset and Extract 2 AIVDM message from Dataset and get AIS Signals
filename_dataset = 'D:\0-MSc_Course\MSc Project\AIS_Python\ais_msgs_sample_V1.0.txt';
dataset = importdata(filename_dataset);
index = randperm(size(dataset,1),2);
% index = [76,214];
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
% delay=2732;
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
% ch=[0.2,0.7;0.5,0.2];
% ch = [0.7483,0.4432;0.3005,0.7891];
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
n1 = 0.5*randn(1,size(x,2));
n2 = 0.5*randn(1,size(x,2));
x(1,:) = x(1,:) + n1;
x(2,:) = x(2,:) + n2;
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

%%
% p=1;
C = size(z,1); % number of components
M = size(z,2); % number of samples
h=C;
N=10;
beta =100;
mutation_rate=0.08;

n = 2;

s_t = 100;
w = zeros(h,h);
Abs = rand(N,2);
for p=1:C
    step = 0;
    % Condition to stopping 
    while step<=s_t
        o = randn(1,size(z,2));
        % Affinity
        for i=1:N
           b_til_k = Abs(i,1:2)*z;
           G = 3*(mean(b_til_k.^2))^2;
           
           
%            b_til = Abs(i,1:2)*z;
%            G = (mean(b_til.^4) -  mean(o.^4))^2;
           
           Abs(i,3) = G;
        end

        Abs_sorted = sortrows(Abs,3);
        b_A = Abs_sorted(end-n+1:end,:);
        
        C_A=[];
        for i=1:size(b_A,1)
            clone_num = floor(beta*N/(n-i+1));
            clone_tmp = repmat(b_A(i,:),clone_num,1);
            C_A=cat(1,C_A,clone_tmp);
        end
        
        
        for pi=1:size(C_A,1)
            
            if rand() <= C_A(pi,3)/(mutation_rate * 100)
                for col=1:h
                    if rand() <= C_A(pi,3)/(mutation_rate * 100)
                        C_A(pi,col) = rand();
                    else
                        ;
                    end
                end
                                
            else
                ;
            end

            
            
            b_til_2_k = C_A(pi,1:2)*z;
            G_aff = 3*(mean(b_til_2_k.^2))^2;
           
%             o = randn(1,size(z,2));
%             b_til_2 = C_A(pi,1:2)*z;
%             G_aff = (mean(b_til_2.^4) -  mean(o.^4))^2;
            C_A(pi,3) = G_aff;
          
        end
        
        b_A = cat(1,b_A,C_A);
        b_A = sortrows(b_A,3);
        b_A = b_A(end-n+1:end,:);
        
        Abs_rand = [rand(N-n,2),zeros(N-n,1)];
        b_A = cat(1,b_A,Abs_rand);


        step = step+1;
    end
    
    
    b_A = sortrows(b_A,3);
    w_p = b_A(end,1:2);
    
    summ = [0,0];
    for j=1:p-1
       summ = summ + (w_p*w(j,:)')*w(j,:);
    end
    w_p = w_p - summ;

    w_p = w_p/norm(w_p);
    w(p,:) = w_p;

end

%%

B=w*z;

figure();plot(-s1);hold on;
plot(B(2,:))
title('Component 1');xlabel('Samples');ylabel('Amplitude')
set(gca,'fontname','Times New Roman','fontsize',9)  
axis tight

figure();plot(-s2);hold on;
plot(B(1,:))
title('Component 2');xlabel('Samples');ylabel('Amplitude')
set(gca,'fontname','Times New Roman','fontsize',9)  
axis tight

%%
recovered_sig1 = movmean(B(2,:),64);
recovered_sig2 = movmean(B(1,:),64);

figure();plot(recovered_sig1);hold on; plot(-s1)
figure();plot(recovered_sig2);hold on; plot(-s2)
    
