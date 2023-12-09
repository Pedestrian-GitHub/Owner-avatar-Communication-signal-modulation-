clc,clear,close all;

% 参数设置
M = 4; % M进制PAM
% total_t = 0.0001;
d = 2; % 
Rs = 100000; % 符号速率
Rb = Rs*log2(M); % 比特速率
fc = 1e6; % 载频
Fs = 10*fc; % 采样频率
T = 1/Rs; % 每符号脉冲持续时间
Tn = round(1/Rs*Fs); % 每符号数据点数
% N = round(Fs*total_t); % 数据点数
% N_symbol = total_t*Rs;
SNR = 10; % 信噪比，以分贝（dB）为单位

% % 生成随机数据
% signal_raw = randi([0 1], 1, N_symbol*2);

I = imread('test_picture1.jpg');
figure;
imshow(I); % 显示图像
title('发送的图片');
disp('图片读取完成.../n');

[A B C] = size(I);
c111 = dec2bin(I(10,1,1));
size(c111)
binary_image = zeros(1,8,A,B,C);
% binary_image = rgb2gray(I);
% figure;
% imshow(binary_image);
temp = zeros(1,8);
for c = 1:C
    for b = 1:B
        for a = 1:A
            temp = cellstr(dec2bin(I(a,b,c)));
            temp = char(temp);
            s = size(temp);
            for i = 1:s(2)
                binary_image(1,8-s(2)+i,a,b,c) = str2num(temp(i));
            end
        end
    end
end
disp('图片转换为二进制.../n');
% fileID = fopen('BINARY_IMAGE.mat','r');
% binary_image = fread(fileID);
% fclose(fileID);

N_symbol = 1*8*A*B*C/2;
total_t = N_symbol/Rs;
N = round(Fs*total_t); % 数据点数
signal_raw = binary_image;

% % 生成M进制PAM信号
% pamSignal = (2*signal_raw - (M-1)) / (M-1);

% 定义时间轴
t = (0: N-1)/Fs;

% 定义g(t)
signal_gt = cos(2*pi*fc*t);
Eg = sum(signal_gt(1:Tn).^2);

% 求出正交基信号
Energy = sum(signal_gt(1:Tn).^2);
ft = signal_gt(1:Tn)/sqrt(Energy/2);

% 初始化PAM
pamSignal_sin = zeros(1,N);

% 脉冲幅度调制
% raw_data = randi([1 M], 1, Rs);
disp('信号开始调制.../n');   
modulation_data = [-3 -1 1 3];
for c = 1:C
    for b = 1:B
        for a = 1:A
            for i = 1:4
                if signal_raw(1,2*i-1,a,b,c)==0 && signal_raw(1,2*i,a,b,c)==0
                    pamSignal_sin(Tn*((c-1)*(A*B*4)+(b-1)*A*4+(a-1)*4+(i-1))+1:Tn*((c-1)*(A*B*4)+(b-1)*A*4+(a-1)*4+i)) = signal_gt(1:Tn).*modulation_data(1);
                elseif signal_raw(1,2*i-1,a,b,c)==0 && signal_raw(1,2*i,a,b,c)==1
                    pamSignal_sin(Tn*((c-1)*(A*B*4)+(b-1)*A*4+(a-1)*4+(i-1))+1:Tn*((c-1)*(A*B*4)+(b-1)*A*4+(a-1)*4+i)) = signal_gt(1:Tn).*modulation_data(2);
                elseif signal_raw(1,2*i-1,a,b,c)==1 && signal_raw(1,2*i,a,b,c)==0
                    pamSignal_sin(Tn*((c-1)*(A*B*4)+(b-1)*A*4+(a-1)*4+(i-1))+1:Tn*((c-1)*(A*B*4)+(b-1)*A*4+(a-1)*4+i)) = signal_gt(1:Tn).*modulation_data(3);
                else
                    pamSignal_sin(Tn*((c-1)*(A*B*4)+(b-1)*A*4+(a-1)*4+(i-1))+1:Tn*((c-1)*(A*B*4)+(b-1)*A*4+(a-1)*4+i)) = signal_gt(1:Tn).*modulation_data(4);
                end
            end
        end
    end
end
disp('信号调制完成.../n');

% 定义能量Eb
Eb = sum(pamSignal_sin.^2)/N_symbol;

% % 求出正交基信号
% Energy = sum(pamSignal_sin.^2);
% ft = pamSignal_sin/sqrt(Energy);

% 通过AWGN信道传输
pamSignalNoisy = awgn(pamSignal_sin, SNR, 'measured');

% %% 时域波形绘制 - 噪声前后对比
% figure;
% plot(t, pamSignal_sin);
% title('M进制PAM信号的时域波形（无噪声）');
% xlabel('时间 (秒)');
% ylabel('幅度');
% axis tight;
% 
% figure;
% % plot(t, ft);
% plot(1:Tn, ft);
% title('M进制PAM信号的正交基信号波形');
% xlabel('时间 (秒)');
% ylabel('幅度');
% axis tight;
% 
% figure;
% subplot(2,1,1);
% plot(t, pamSignal_sin);
% title('M进制PAM信号的时域波形（无噪声）');
% xlabel('时间 (秒)');
% ylabel('幅度');
% 
% subplot(2,1,2);
% plot(t, pamSignalNoisy);
% title('M进制PAM信号的时域波形（有噪声）');
% xlabel('时间 (秒)');
% ylabel('幅度');

%% 频谱计算 - 噪声前后对比
f = Fs*(0:(N/2))/N;
PAM_FFT = fft(pamSignal_sin);
PAM_FFT_Noisy = fft(pamSignalNoisy);
P2 = abs(PAM_FFT/N);
P2_Noisy = abs(PAM_FFT_Noisy/N);
P1 = P2(1:N/2+1);
P1_Noisy = P2_Noisy(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);
P1_Noisy(2:end-1) = 2*P1_Noisy(2:end-1);

% % 频谱绘制 - 噪声前后对比
% figure;
% subplot(2,1,1);
% plot(f, P1);
% title('M进制PAM信号的频谱（无噪声）');
% xlabel('频率 (Hz)');
% ylabel('|P1(f)|');
% 
% subplot(2,1,2);
% plot(f, P1_Noisy);
% title('M进制PAM信号的频谱（有噪声）');
% xlabel('频率 (Hz)');
% ylabel('|P1(f)|');

% % figure;
% % subplot(2,1,1);
% % plot(f, P2);
% % title('M进制PAM信号的频谱（无噪声）');
% % xlabel('频率 (Hz)');
% % ylabel('|P2(f)|');
% % 
% % subplot(2,1,2);
% % plot(f, P2_Noisy);
% % title('M进制PAM信号的频谱（有噪声）');
% % xlabel('频率 (Hz)');
% % ylabel('|P2(f)|');


%% 最佳接收机设计
disp('信号开始解调.../n');
decode_signal = zeros(1,8,A,B,C);
receive_image = zeros(A,B,C);
for c = 1:C
    for b = 1:B
        for a = 1:A
            for i = 1:4
                receive_signal = pamSignalNoisy(Tn*((c-1)*(A*B*4)+(b-1)*A*4+(a-1)*4+(i-1))+1:Tn*((c-1)*(A*B*4)+(b-1)*A*4+(a-1)*4+i));
                channel = conv(conj(ft), receive_signal);
                result = channel(Tn);

                if result <= -2*sqrt(Eg/2)
                    decode_signal(1,2*i-1:2*i,a,b,c) = [0,0];
                elseif result > -2*sqrt(Eg/2) && result <= 0
                    decode_signal(1,2*i-1:2*i,a,b,c) = [0,1];
                elseif result > 0 && result <= 2*sqrt(Eg/2)
                    decode_signal(1,2*i-1:2*i,a,b,c) = [1,0];
                else
                    decode_signal(1,2*i-1:2*i,a,b,c) = [1,1];
                end
            end
            temp = bi2de(decode_signal(1,:,a,b,c), 'left-msb');
            receive_image(a,b,c) = uint8(temp); 
        end
    end
end
disp('信号解调完成.../n');
% figure;
% % plot(1:(N+Tn-1),channel);
% plot(1:(2*Tn-1),channel);
figure;
imshow(uint8(receive_image));
title('接收到的图片');


% 计算误码率（BER）
disp('计算误码率.../n');
errors = 0;
for c = 1:C
    for b = 1:B
        for a = 1:A
            errors = errors + sum(signal_raw(1,:,a,b,c) ~= decode_signal(1,:,a,b,c));
        end
    end
end
BER = errors/(N_symbol*2);
% 显示误码率
disp(['BER: ', num2str(BER)]);

% figure;
% plot(t, filteredSignal);
% title('相关接收后PAM信号的时域波形');
% xlabel('时间 (秒)');
% ylabel('幅度');

