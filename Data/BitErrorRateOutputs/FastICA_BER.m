% FastICA Algorithm to Separate 2 overlapped AIS Signals
clc;
clear all;
close all;
%%
SNRs = 1:2:27;
ber_11 = zeros(250,length(SNRs));
ber_22 = zeros(250,length(SNRs));
for hhh = 1:size(ber_11,1)
for kkk=1:size(ber_11,2)
%% Read Dataset and Extract 2 AIVDM message from Dataset and get AIS Signals
filename_dataset = 'D:\0-MSc_Course\MSc Project\AIS_Python\ais_msgs_sample_V1.0.txt';
dataset = importdata(filename_dataset);
index = randperm(size(dataset,1),2);
%index = [747,736];
aismsg1=char(dataset(index(1)));
aismsg2=char(dataset(index(2)));

[sigI1,sigQ1,sigIQ1,stored_msgb1]=ais_gmsk_signal_gen(aismsg1);
[sigI2,sigQ2,sigIQ2,stored_msgb2]=ais_gmsk_signal_gen(aismsg2);

% figure();
% subplot(2,1,1)
% plot(sigI1)
% title('In-phase');xlabel('Samples');ylabel('Amplitude')
% set(gca,'XLimSpec', 'Tight','fontname','Times New Roman','fontsize',9)
% 
% subplot(2,1,2)
% plot(sigQ1)
% title('Quadrature');xlabel('Samples');ylabel('Amplitude')
% set(gca,'XLimSpec', 'Tight','fontname','Times New Roman','fontsize',9)  
% 
% figure();
% subplot(2,1,1)
% plot(sigI2)
% title('In-Phase');xlabel('Samples');ylabel('Amplitude')
% set(gca,'XLimSpec', 'Tight','fontname','Times New Roman','fontsize',9)  
% 
% subplot(2,1,2)
% plot(sigQ2)
% title('Quadrature');xlabel('Samples');ylabel('Amplitude')
% set(gca,'XLimSpec', 'Tight','fontname','Times New Roman','fontsize',9)  


same_len = max(size(sigI1,2),size(sigI2,2));
min_len = min(size(sigI1,2),size(sigI2,2));
delay = randi([60,3000]);
%delay=2368;

s1=[sigI1,zeros(1,20000)];
sq1=[sigQ1,zeros(1,20000)];
s2=[zeros(1,delay),sigI2,zeros(1,20000)];
sq2=[zeros(1,delay),sigQ2,zeros(1,20000)];
s1 = s1(1:same_len);
sq1 = sq1(1:same_len);
s2 = s2(1:same_len);
sq2 = sq2(1:same_len);

% figure();
% subplot(2,1,1)
% plot(s1)
% title('Signal 1 (In-Phase)');xlabel('Samples');ylabel('Amplitude')
% set(gca,'XLimSpec', 'Tight','fontname','Times New Roman','fontsize',9)  
% 
% subplot(2,1,2)
% plot(s2)
% title('Signal 2 (In-Phase)');xlabel('Samples');ylabel('Amplitude')
% set(gca,'XLimSpec', 'Tight','fontname','Times New Roman','fontsize',9)  
% 
% figure();
% subplot(2,1,1)
% plot(sq1)
% title('Signal 1 (Quadrature)');xlabel('Samples');ylabel('Amplitude')
% set(gca,'XLimSpec', 'Tight','fontname','Times New Roman','fontsize',9)  
% 
% subplot(2,1,2)
% plot(sq2)
% title('Signal 2 (Quadrature)');xlabel('Samples');ylabel('Amplitude')
% set(gca,'XLimSpec', 'Tight','fontname','Times New Roman','fontsize',9)  

%% Generation Hybrid channel and mixedsignal with 2x2 mixing matrix
ch = [0.75,0.44;0.30,0.79];
x = ch*[s1;s2];
xq=ch*[sq1;sq2];

% figure();
% subplot(2,1,1)
% plot(x(1,:))
% title('Antenna 1 (I)');xlabel('Samples');ylabel('Amplitude')
% set(gca,'XLimSpec', 'Tight','fontname','Times New Roman','fontsize',9)  
% 
% subplot(2,1,2)
% plot(x(2,:))
% title('Antenna 2 (I)');xlabel('Samples');ylabel('Amplitude')
% set(gca,'XLimSpec', 'Tight','fontname','Times New Roman','fontsize',9)  
% 
% figure();
% subplot(2,1,1)
% plot(xq(1,:))
% title('Antenna 1 (Q)');xlabel('Samples');ylabel('Amplitude')
% set(gca,'XLimSpec', 'Tight','fontname','Times New Roman','fontsize',9)  
% 
% subplot(2,1,2)
% plot(xq(2,:))
% title('Antenna 2 (Q)');xlabel('Samples');ylabel('Amplitude')
% set(gca,'XLimSpec', 'Tight','fontname','Times New Roman','fontsize',9) 

n1=0;
nq1=0;
n2=0;
nq2=0;

%% Noise
snr = SNRs(kkk); % [dB]
Ps_Pn = 10^(snr/10);
n1=(1/sqrt(Ps_Pn))*randn(1,size(x,2))*sqrt(mean(x(1,:).^2));
nq1=(1/sqrt(Ps_Pn))*randn(1,size(xq,2))*sqrt(mean(xq(1,:).^2));
n2=(1/sqrt(Ps_Pn))*randn(1,size(x,2))*sqrt(mean(x(2,:).^2));
nq2=(1/sqrt(Ps_Pn))*randn(1,size(xq,2))*sqrt(mean(xq(2,:).^2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

x_noisy(1,1:min_len) = x(1,1:min_len) + n1(1,1:1:min_len);
x_noisy(2,1:min_len) = x(2,1:min_len) + n2(1,1:1:min_len);

xq_noisy(1,1:min_len) = xq(1,1:min_len) + nq1(1,1:1:min_len);
xq_noisy(2,1:min_len) = xq(2,1:min_len) + nq2(1,1:1:min_len);

%%
x_noisy(1,1:min_len) = x(1,1:min_len) + n1(1,1:1:min_len);
x_noisy(2,1:min_len) = x(2,1:min_len) + n2(1,1:1:min_len);
xq_noisy(1,1:min_len) = xq(1,1:min_len) + nq1(1,1:1:min_len);
xq_noisy(2,1:min_len) = xq(2,1:min_len) + nq2(1,1:1:min_len);

%% centering data and whitening
m1 = mean(x_noisy(1,:));
m2 = mean(x_noisy(2,:));
mq1 = mean(xq_noisy(1,:));
mq2 = mean(xq_noisy(2,:));

x_noisy(1,:) = x_noisy(1,:)-m1;
x_noisy(2,:) = x_noisy(2,:)-m2;
xq_noisy(1,:) = xq_noisy(1,:)-mq1;
xq_noisy(2,:) = xq_noisy(2,:)-mq2;

[EigenVectors,eigenValues] = eig(x_noisy*transpose(x_noisy));
[EigenVectors_q,eigenValues_q] = eig(xq_noisy*transpose(xq_noisy));

D_12 = eigenValues^(-1/2);
D_12_q = eigenValues_q^(-1/2);

E_T=transpose(EigenVectors);
E_T_q=transpose(EigenVectors_q);

z = D_12*E_T*x_noisy;
zq = D_12_q*E_T_q*xq_noisy;

z(1,:) = z(1,:)/std(z(1,:));
z(2,:) = z(2,:)/std(z(2,:));

zq(1,:) = zq(1,:)/std(zq(1,:));
zq(2,:) = zq(2,:)/std(zq(2,:));

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
%        disp(w_p_old)
%        disp(w_p_new)
       
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
%     disp(stepp)
    w(:,p)=w_p_new;
end

%%
B =w*z;

% figure();
% plot(s1);hold on;
% plot(B(2,:))
% title('Component 1 (I)');xlabel('Samples');ylabel('Amplitude')
% set(gca,'fontname','Times New Roman','fontsize',9)  
% axis tight
% 
% figure();plot(-s2);hold on;
% plot(B(1,:))
% title('Component 2 (I)');xlabel('Samples');ylabel('Amplitude')
% set(gca,'fontname','Times New Roman','fontsize',9)  
% axis tight
%%
Bq =w*zq;

% figure();
% plot(sq1);hold on;
% plot(Bq(2,:))
% title('Component 1 (Q)');xlabel('Samples');ylabel('Amplitude')
% set(gca,'fontname','Times New Roman','fontsize',9)  
% axis tight
% 
% figure();plot(-sq2);hold on;
% plot(Bq(1,:))
% title('Component 2 (Q)');xlabel('Samples');ylabel('Amplitude')
% set(gca,'fontname','Times New Roman','fontsize',9)  
% axis tight

%%
%Demodolator


%if  ((var(B(2,1:delay))) < (var(B(1,1:delay))) && (var(Bq(2,1:delay))) < (var(Bq(1,1:delay)))) || ((abs(mean(B(2,1:delay)))) < (abs(mean(B(1,1:delay)))) && (abs(mean(Bq(2,1:delay)))) < (abs(mean((Bq(1,1:delay))))))
%if  (abs(mean(B(2,1:delay))) < abs(mean(B(1,1:delay)))) && (var(B(2,1:delay)) < var(B(1,1:delay)))
% if  (abs(mean(B(2,1:delay))) < abs(mean(B(1,1:delay)))) && (abs(mean(Bq(2,1:delay))) < abs(mean(Bq(1,1:delay))))
if  (norm(B(2,1:delay)) < norm(B(1,1:delay))) 
    %
    test1 = [mean(Bq(1,1:delay)),mean(Bq(2,1:delay)),mean(B(1,1:delay)),mean(B(2,1:delay))];
    test2 = [var(Bq(1,1:delay)),var(Bq(2,1:delay)),var(B(1,1:delay)),var(B(2,1:delay))];
    test3 = [std(Bq(1,1:delay)),std(Bq(2,1:delay)),std(B(1,1:delay)),std(B(2,1:delay))];
    disp("flag 1");disp([hhh,kkk]);disp(test1);disp(test2);disp(test3);
q_part_sig1 = Bq(1,1:min_len);
i_part_sig1 = B(1,1:min_len);

q_part_sig2 = Bq(2,1:min_len);
i_part_sig2 = B(2,1:min_len);

else 
        test11 = [mean(Bq(1,1:delay)),mean(Bq(2,1:delay)),mean(B(1,1:delay)),mean(B(2,1:delay))];
    test22 = [var(Bq(1,1:delay)),var(Bq(2,1:delay)),var(B(1,1:delay)),var(B(2,1:delay))];
    test33 = [std(Bq(1,1:delay)),std(Bq(2,1:delay)),std(B(1,1:delay)),std(B(2,1:delay))];
    disp("flag 2");disp([hhh,kkk]);disp(test11);disp(test22);disp(test33)
q_part_sig1 = Bq(2,1:min_len);
i_part_sig1 = B(2,1:min_len);

q_part_sig2 = Bq(1,1:min_len);
i_part_sig2 = B(1,1:min_len);  

end

q_part_mov_sig1 = movmean(q_part_sig1,4);
i_part_mov_sig1 = movmean(i_part_sig1,4);
    
q_part_mov_sig2 = movmean(q_part_sig2,4);
i_part_mov_sig2 = movmean(i_part_sig2,4);

Arctan1 = zeros(1,min_len);
Arctan2 = zeros(1,min_len);
for jj=1:min_len

    Arctan1(1,jj) = atan(q_part_mov_sig1(jj)/i_part_mov_sig1(jj));
    
end
for jj=delay:min_len
    Arctan2(1,jj) = atan(q_part_mov_sig2(jj)/i_part_mov_sig2(jj));
end
diffrentiate1 = diff(Arctan1);
diffrentiate2 = diff(Arctan2);

% figure();plot(diffrentiate1);title('diferentiate............')
% figure();plot(diffrentiate2)
%%
extracted_bits1 = zeros(1,256);
sign_dif1 = sign(diffrentiate1);
sign_dif2 = sign(diffrentiate2);
for i=1:min_len/64 -1 
extracted_bits1(1,i) =sum((sign_dif1(32+(64*(i-1)+1:64*i))));
end
extracted_bits2 = zeros(1,256);
for i = 1:floor((16448-delay)/64) -1 
extracted_bits2(1,i) =sum((sign_dif2(delay+32+(64*(i-1)+1:64*i))));
end

bits_sign1 = sign(extracted_bits1(1:256));
bits_sign2 = sign(extracted_bits2(1:floor((16448-delay)/64)));
% bits_sign2 = sign(extracted_bits2(1:256));
%% %BER
first_nrz_bits1 =(-1).^(stored_msgb1+1);
first_nrz_bits2 =(-1).^(stored_msgb2+1);

checker1 = (first_nrz_bits1(1:255) == bits_sign1(1:255));
checker2 = (first_nrz_bits2(1:floor((16448-delay)/64)-1) == bits_sign2(1:floor((16448-delay)/64)-1));
% checker2 = (first_nrz_bits2(1:floor((16448)/64)-1) == bits_sign2(1:floor((16448)/64)-1));

ber1 = (numel(checker1) - nnz(checker1))/numel(checker1);
ber2 = (numel(checker2) - nnz(checker2))/numel(checker2);
% ber_10(hhh,kkk) = ber;

ber_11(hhh,kkk) = ber1;
ber_22(hhh,kkk) = ber2;
end
end

% disp('(Complete Signal) number error bits:');disp(numel(checker1) - nnz(checker1))
% disp('ber rate 1:');disp(ber1)
% 
% disp('(Delayed Signal) number error bits:');disp(numel(checker2) - nnz(checker2))
% disp('ber rate 2:');disp(ber2)
