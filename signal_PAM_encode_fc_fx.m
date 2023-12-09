clc,clear,close all;

% 参数设置
M = 4; % M进制PAM
total_t = 0.0001;
d = 2; % 
Rs = 100000; % 符号速率
Rb = Rs*log2(M); % 比特速率
fc = 1e6; % 载频
Fs = 10*fc; % 采样频率
T = 1/Rs; % 每符号脉冲持续时间
Tn = round(1/Rs*Fs); % 每符号数据点数
N = round(Fs*total_t); % 数据点数
N_symbol = total_t*Rs;
SNR = 10; % 信噪比，以分贝（dB）为单位

% 生成随机数据
signal_raw = randi([0 1], 1, N_symbol*2);

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
modulation_data = [-3 -1 1 3];
for i = 1:N_symbol
    if signal_raw(2*i-1)==0 && signal_raw(2*i)==0
        pamSignal_sin(Tn*(i-1)+1:Tn*i) = signal_gt(Tn*(i-1)+1:Tn*i).*modulation_data(1);
    elseif signal_raw(2*i-1)==0 && signal_raw(2*i)==1
        pamSignal_sin(Tn*(i-1)+1:Tn*i) = signal_gt(Tn*(i-1)+1:Tn*i).*modulation_data(2);
    elseif signal_raw(2*i-1)==1 && signal_raw(2*i)==0
        pamSignal_sin(Tn*(i-1)+1:Tn*i) = signal_gt(Tn*(i-1)+1:Tn*i).*modulation_data(3);
    else
        pamSignal_sin(Tn*(i-1)+1:Tn*i) = signal_gt(Tn*(i-1)+1:Tn*i).*modulation_data(4);
    end
end

% 定义能量Eb
Eb = sum(pamSignal_sin.^2)/N_symbol;

% % 求出正交基信号
% Energy = sum(pamSignal_sin.^2);
% ft = pamSignal_sin/sqrt(Energy);

% 通过AWGN信道传输
pamSignalNoisy = awgn(pamSignal_sin, SNR, 'measured');

%% 时域波形绘制 - 噪声前后对比
figure;
plot(t, pamSignal_sin);
title('M进制PAM信号的时域波形（无噪声）');
xlabel('时间 (秒)');
ylabel('幅度');
axis tight;

figure;
% plot(t, ft);
plot(1:Tn, ft);
title('M进制PAM信号的正交基信号波形');
xlabel('时间 (秒)');
ylabel('幅度');
axis tight;

figure;
subplot(2,1,1);
plot(t, pamSignal_sin);
title('M进制PAM信号的时域波形（无噪声）');
xlabel('时间 (秒)');
ylabel('幅度');

subplot(2,1,2);
plot(t, pamSignalNoisy);
title('M进制PAM信号的时域波形（有噪声）');
xlabel('时间 (秒)');
ylabel('幅度');

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

% 频谱绘制 - 噪声前后对比
figure;
subplot(2,1,1);
plot(f, P1);
title('M进制PAM信号的频谱（无噪声）');
xlabel('频率 (Hz)');
ylabel('|P1(f)|');

subplot(2,1,2);
plot(f, P1_Noisy);
title('M进制PAM信号的频谱（有噪声）');
xlabel('频率 (Hz)');
ylabel('|P1(f)|');

% figure;
% subplot(2,1,1);
% plot(f, P2);
% title('M进制PAM信号的频谱（无噪声）');
% xlabel('频率 (Hz)');
% ylabel('|P2(f)|');
% 
% subplot(2,1,2);
% plot(f, P2_Noisy);
% title('M进制PAM信号的频谱（有噪声）');
% xlabel('频率 (Hz)');
% ylabel('|P2(f)|');


%% 最佳接收机设计

decode_signal = zeros(1,2*N_symbol);
for i = 1:N_symbol
    receive_signal = pamSignalNoisy(Tn*(i-1)+1:Tn*i);
    channel = conv(conj(ft), receive_signal);
    result = channel(Tn);

%     match_result = max(max(max(result_1, result_2), result_3), result_4);
    if result <= -2*sqrt(Eg/2)
        decode_signal(2*i-1:2*i) = [0,0];
    elseif result > -2*sqrt(Eg/2) && result <= 0
        decode_signal(2*i-1:2*i) = [0,1];
    elseif result > 0 && result <= 2*sqrt(Eg/2)
        decode_signal(2*i-1:2*i) = [1,0];
    else
        decode_signal(2*i-1:2*i) = [1,1];
    end
end
figure;
% plot(1:(N+Tn-1),channel);
plot(1:(2*Tn-1),channel);
% subplot(2,2,2);
% plot(1:2*Tn-1,channel_2);
% subplot(2,2,3);
% plot(1:2*Tn-1,channel_3);
% subplot(2,2,4);
% plot(1:2*Tn-1,channel_4);


% 计算误码率（BER）
errors = sum(signal_raw ~= decode_signal);
BER = errors/(N_symbol*2);
% 显示误码率
disp(['BER: ', num2str(BER)]);

% figure;
% plot(t, filteredSignal);
% title('相关接收后PAM信号的时域波形');
% xlabel('时间 (秒)');
% ylabel('幅度');

