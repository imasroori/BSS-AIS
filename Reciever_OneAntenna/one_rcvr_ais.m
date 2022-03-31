% Oone Algorithm
clc;
clear all;
close all;

%%


filename_dataset = 'D:\0-MSc_Course\MSc Project\AIS_Python\ais_msgs_sample_V1.0.txt';
dataset = importdata(filename_dataset);
index = randperm(size(dataset,1),2);
index = [76,214];
aismsg1=char(dataset(index(1)));
aismsg2=char(dataset(index(2)));

[sigI1,sigQ1,sigIQ1,stored_msgb1]=ais_gmsk_signal_gen(aismsg1);
[sigI2,sigQ2,sigIQ2,stored_msgb2]=ais_gmsk_signal_gen(aismsg2);

figure();
subplot(2,1,1)
plot(sigI1)
title('In-Phase');xlabel('Samples');ylabel('Amplitude')
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
delay=2732;
s1=[sigI1,zeros(1,20000)];
sq1=[sigQ1,zeros(1,20000)];
s2=[zeros(1,delay),sigI2,zeros(1,20000)];
sq2=[zeros(1,delay),sigQ2,zeros(1,20000)];

s1 = s1(1:same_len);
sq1 = sq1(1:same_len);
s2 = s2(1:same_len);
sq2 = sq2(1:same_len);

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
ch = [0.8,0.65];
% ch = [0.7483,0.4432;0.3005,0.7891];
% ch=[0.2,0.7;0.5,0.2];
x = ch*[s1;s2];
xq = ch*[sq1;sq2];
n=0;
nq=0;
%% Noise
snr = 15; % [dB]
Ps_Pn = 10^(snr/10);
n=(1/sqrt(Ps_Pn))*randn(1,size(x,2))*sqrt(mean(x.^2));
nq=(1/sqrt(Ps_Pn))*randn(1,size(xq,2))*sqrt(mean(xq.^2));
x_noisy = x+n;
xq_noisy = xq+nq;

%%
x_noisy = x+n;
xq_noisy = xq+nq;
figure();
plot(x_noisy(1,:))
title('Mixed-in Signals');xlabel('Samples');ylabel('Amplitude')
set(gca,'XLimSpec', 'Tight','fontname','Times New Roman','fontsize',9)


coef_mat = zeros(1,80);
for i=1:80
    coef_mat(i) = estmate_coef(pre_shape(1:32*i),x_noisy(1:32*i));
end
coef = mean(coef_mat);
xx = [(1/coef)*x_noisy,zeros(1,1000)];
xx_q = [(1/coef)*xq_noisy,zeros(1,1000)];

recovered_sig = zeros(1,20000);
recovered_sig_q = zeros(1,20000);
recovered_sig(1:2592) = pre_shape(1:2592);
recovered_sig_q(1:2592) = pre_shape_q(1:2592);
%%

stepp=1;
for i=2593:128:size(x_noisy,2)
    if stepp == 1
        idx1 = norm(xx(i:i+127)-sh1);
        idx2 = norm(xx(i:i+127)-sh2);
        idx3 = norm(xx(i:i+127)-sh3);
        idx4 = norm(xx(i:i+127)-sh4);
        
        shapes = [sh1,idx1;
                  sh2,idx2;
                  sh3,idx3;
                  sh4,idx4];
              
        tmp = sortrows(shapes,129);
        if tmp(1,1:128) == sh1
            rec_num(1,stepp) = 1;
        elseif tmp(1,1:128) == sh2
            rec_num(1,stepp) = 2;
        elseif tmp(1,1:128) == sh3
            rec_num(1,stepp) = 3;
        else
            rec_num(1,stepp) = 4;
        end
          
        recovered_sig(i:i+127) = tmp(1,1:128);
    else
        if (rec_num(1,stepp-1) == 1) ||(rec_num(1,stepp-1) == 4)
            idx1 = norm(xx(i:i+127)-sh1);
            idx2 = norm(xx(i:i+127)-sh2);
            shapes = [sh1,idx1;
                      sh2,idx2];
            tmp = sortrows(shapes,129);
        elseif (rec_num(1,stepp-1) == 2) || (rec_num(1,stepp-1) == 3)
            idx3 = norm(xx(i:i+127)-sh3);
            idx4 = norm(xx(i:i+127)-sh4);
            shapes = [sh3,idx3;
                      sh4,idx4];
            tmp = sortrows(shapes,129);
            
        end        
        
        
        if tmp(1,1:128) == sh1
            rec_num(1,stepp) = 1;
        elseif tmp(1,1:128) == sh2
            rec_num(1,stepp) = 2;
        elseif tmp(1,1:128) == sh3
            rec_num(1,stepp) = 3;
        else
            rec_num(1,stepp) = 4;
        end
        
                
        recovered_sig(i:i+127) = tmp(1,1:128);
    end
    stepp = stepp + 1;
end


%%
figure();
subplot(3,1,1)
plot([x_noisy(1:16448)])
title('Primary Interfered Signal');xlabel('Samples');ylabel('Amplitude')
set(gca,'XLimSpec', 'Tight','fontname','Times New Roman','fontsize',9)  

subplot(3,1,2)
plot(recovered_sig(1:16448))
title('Extracted Component 1');xlabel('Samples');ylabel('Amplitude')
set(gca,'XLimSpec', 'Tight','fontname','Times New Roman','fontsize',9) 
hold on; plot(s1)

subplot(3,1,3)
plot(xx(1:16448)-recovered_sig(1:16448))
title('Reminder of Interfered Signal (Component 2)');xlabel('Samples');ylabel('Amplitude')
set(gca,'XLimSpec', 'Tight','fontname','Times New Roman','fontsize',9) 
hold on; plot(s2)

%%
%Quadrature

stepp_q=1;
for i=2529:128:size(xq_noisy,2)
    if stepp_q == 1
        idxq1 = norm(xx_q(i:i+127)-sh1);
        idxq2 = norm(xx_q(i:i+127)-sh2);
        idxq3 = norm(xx_q(i:i+127)-sh3);
        idxq4 = norm(xx_q(i:i+127)-sh4);
        
        shapes_q = [sh1,idxq1;
                    sh2,idxq2;
                    sh3,idxq3;
                    sh4,idxq4];      
        tmp_q = sortrows(shapes_q,129);
        
        if tmp_q(1,1:128) == sh1
            rec_num_q(1,stepp_q) = 1;
        elseif tmp_q(1,1:128) == sh2
            rec_num_q(1,stepp_q) = 2;
        elseif tmp_q(1,1:128) == sh3
            rec_num_q(1,stepp_q) = 3;
        else
            rec_num_q(1,stepp_q) = 4;
        end
                
        recovered_sig_q(i:i+127) = tmp_q(1,1:128);
    else
         if (rec_num_q(1,stepp_q-1) == 1) ||(rec_num_q(1,stepp_q-1) == 4)
            idxq1 = norm(xx_q(i:i+127)-sh1);
            idxq2 = norm(xx_q(i:i+127)-sh2);
            shapes_q = [sh1,idxq1;
                      sh2,idxq2];
            tmp_q = sortrows(shapes_q,129);
        elseif (rec_num_q(1,stepp_q-1) == 2) || (rec_num_q(1,stepp_q-1) == 3)
            idxq3 = norm(xx_q(i:i+127)-sh3);
            idxq4 = norm(xx_q(i:i+127)-sh4);
            shapes_q = [sh3,idxq3;
                        sh4,idxq4];
            tmp_q = sortrows(shapes_q,129);
            
        end       
        
        
               
        
        if tmp_q(1,1:128) == sh1
            rec_num_q(1,stepp_q) = 1;
        elseif tmp_q(1,1:128) == sh2
            rec_num_q(1,stepp_q) = 2;
        elseif tmp_q(1,1:128) == sh3
            rec_num_q(1,stepp_q) = 3;
        else
            rec_num_q(1,stepp_q) = 4;
        end
                
        
        recovered_sig_q(i:i+127) = tmp_q(1,1:128);
    end
    stepp_q = stepp_q + 1;
end

%%
%Quadrature
figure();
subplot(3,1,1)
plot([xq_noisy(1:16448)])
title('(Quadrature)-Primary Interfered Signal');xlabel('Samples');ylabel('Amplitude')
set(gca,'XLimSpec', 'Tight','fontname','Times New Roman','fontsize',9)  

subplot(3,1,2)
plot(recovered_sig_q(1:16448))
title('(Quadrature)-Extracted Component 1');xlabel('Samples');ylabel('Amplitude')
set(gca,'XLimSpec', 'Tight','fontname','Times New Roman','fontsize',9) 
hold on; plot(sq1)

subplot(3,1,3)
plot(xx_q(1:16448)-recovered_sig_q(1:16448))
title('(Quadrature)-Reminder of Interfered Signal (Component 2)');xlabel('Samples');ylabel('Amplitude')
set(gca,'XLimSpec', 'Tight','fontname','Times New Roman','fontsize',9) 
hold on; plot(sq2)



%% NOOOO NOOOOO  for step to and run second algorithm
x_new = xx-recovered_sig(1:size(xx,2));
coef_mat2 = zeros(1,floor((3488-delay)/32));
for i=1:size(coef_mat2,2)
    coef_mat2(i) = estmate_coef(pre_shape(1:32*i),x_new(delay+(1:32*i)));
end
coef2 = mean(coef_mat2);

xx_new = [(1/coef2)*x_new,zeros(1,delay+1)];
recovered_sig2 = zeros(1,20000);
recovered_sig2(delay+(1:2592)) = pre_shape(1:2592);

%% NOOOO NOOOOO

stepp2=1;
for i=delay+(2593:128:16512)
    if stepp2 == 1
        idx1 = norm(xx_new(i:i+127)-sh1);
        idx2 = norm(xx_new(i:i+127)-sh2);
        idx3 = norm(xx_new(i:i+127)-sh3);
        idx4 = norm(xx_new(i:i+127)-sh4);
        shapes = [sh1,idx1;
                  sh2,idx2;
                  sh3,idx3;
                  sh4,idx4];
        tmp = sortrows(shapes,129);
        if tmp(1,1:128) == sh1
            rec_num(1,stepp2) = 1;
        elseif tmp(1,1:128) == sh2
            rec_num(1,stepp2) = 2;
        elseif tmp(1,1:128) == sh3
            rec_num(1,stepp2) = 3;
        else
            rec_num(1,stepp2) = 4;
        end
                
        recovered_sig2(i:i+127) = tmp(1,1:128);
    else
        if (rec_num(1,stepp2-1) == 1) ||(rec_num(1,stepp2-1) == 4)
            idx1 = norm(xx_new(i:i+127)-sh1);
            idx2 = norm(xx_new(i:i+127)-sh2);
            shapes = [sh1,idx1;
                      sh2,idx2];
            tmp = sortrows(shapes,129);
        elseif (rec_num(1,stepp2-1) == 2) || (rec_num(1,stepp2-1) == 3)
            idx3 = norm(xx_new(i:i+127)-sh3);
            idx4 = norm(xx_new(i:i+127)-sh4);
            shapes = [sh3,idx3;
                      sh4,idx4];
            tmp = sortrows(shapes,129);
            
        end
        
        if tmp(1,1:128) == sh1
            rec_num(1,stepp2) = 1;
        elseif tmp(1,1:128) == sh2
            rec_num(1,stepp2) = 2;
        elseif tmp(1,1:128) == sh3
            rec_num(1,stepp2) = 3;
        else
            rec_num(1,stepp2) = 4;
        end
                
        recovered_sig2(i:i+127) = tmp(1,1:128);
    end
    stepp2 = stepp2 + 1;
end

%%  NOOOO NOOOOO
figure();
subplot(3,1,1)
plot(x_new(1:16448))
title('Secondary Interfered Signal');xlabel('Samples');ylabel('Amplitude')
set(gca,'XLimSpec', 'Tight','fontname','Times New Roman','fontsize',9)  

subplot(3,1,2)
plot(recovered_sig2(1:16448))
title('Extracted Component 2');xlabel('Samples');ylabel('Amplitude')
set(gca,'XLimSpec', 'Tight','fontname','Times New Roman','fontsize',9) 
hold on; plot(s2)

subplot(3,1,3)
plot(xx_new(1:16448)-recovered_sig2(1:16448))
title('Reminder of Signal (Pure Noise)');xlabel('Samples');ylabel('Amplitude')
set(gca,'XLimSpec', 'Tight','fontname','Times New Roman','fontsize',9) 
hold on; plot(s3)

%%
%Demodolator
q_part_sig2 = xx_q(1:16448)-recovered_sig_q(1:16448);
i_part_sig2 = xx(1:16448)-recovered_sig(1:16448);

figure();
plot(q_part_sig2);
hold on; plot(sq2)
figure();
plot(i_part_sig2);
hold on; plot(s2)

q_part_mov = movmean(q_part_sig2,4);
i_part_mov = movmean(i_part_sig2,4);
figure();
plot(q_part_mov);hold on; plot(q_part_sig2)
figure();
plot(i_part_mov);hold on; plot(i_part_sig2)

Arctan = zeros(1,16448);
Arctan_1 = zeros(1,16448);
for jj=delay:16448

    Arctan(1,jj) = atan(q_part_mov(jj)/i_part_mov(jj));
    
end
for jj=1:16448
    Arctan_1(1,jj) = atan(recovered_sig_q(jj)/recovered_sig(jj));
end
diffrentiate = diff(Arctan);
diffrentiate_1 = diff(Arctan_1);

figure();plot(diffrentiate)
figure();plot(diffrentiate_1)
%%

extracted_bits = zeros(1,256);
sign_dif = sign(diffrentiate);
sign_dif_1 = sign(diffrentiate_1);
for i=1:floor((16448-delay)/64)
extracted_bits(1,i) =sum((sign_dif(2752+(64*(i-1)+1:64*i-1))));
end
extracted_bits_1 = zeros(1,256);
for i = 1:16448/64 -1
    extracted_bits_1(1,i) =sum((sign_dif_1(32+(64*(i-1)+1:64*i-1))));
end

bits_sign = sign(extracted_bits(1:214));
bits_sign_1 = sign(extracted_bits_1(1:256));
%%
%BER

nrz_bits =(-1).^(stored_msgb2+1);
checker = (nrz_bits(1:214) == bits_sign);
ber = (numel(checker) - nnz(checker))/numel(checker);

disp('number error bits:');disp(numel(checker) - nnz(checker))
disp('ber rate:');disp(ber)

%%
SNRs = [1,3,5,7,9,11,13,15,17,19,21,23,25];
BERs = [0.3925,0.2897,0.2336,0.1589,0.1168,0.0561,0.0327,0.00001,0.000001,0,0,0,0];
figure();semilogy(SNRs,BERs)
grid on;
title('BER-SNR (one reciever) and MA-Window= 4');xlabel('SNR[dB]');ylabel('BER')
set(gca,'YLimSpec', 'Tight','XLimSpec', 'Tight','fontname','Times New Roman','fontsize',9) 



