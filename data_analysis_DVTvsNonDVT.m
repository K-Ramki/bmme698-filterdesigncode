clear
close all
clc

simulink
filterDesigner

data_hieu_noDVT_run1 = xlsread("max30102_data_Hieu_noDVT_run1_02-18.xlsx");
data_hieu_noDVT_run2 = xlsread("max30102_data_Hieu_noDVT_run2_02-18.xlsx");
data_hieu_noDVT_run3 = xlsread("max30102_data_Hieu_noDVT_run3_02-18.xlsx");

data_hieu_DVT_run1 = xlsread("max30102_data_Hieu_DVT_run1_02-18.xlsx");
data_hieu_DVT_run2 = xlsread("max30102_data_Hieu_DVT_run2_02-18.xlsx");
data_hieu_DVT_run3 = xlsread("max30102_data_Hieu_DVT_run3_02-18.xlsx");

data_bryce_noDVT_run1 = xlsread("max30102_data_Bryce_noDVT_run1_02-18.xlsx");
data_bryce_noDVT_run2 = xlsread("max30102_data_Bryce_noDVT_run2_02-18.xlsx");
data_bryce_noDVT_run3 = xlsread("max30102_data_Bryce_noDVT_run3_02-18.xlsx");

data_bryce_DVT_run1 = xlsread("max30102_data_Bryce_DVT_run1_02-18.xlsx");
data_bryce_DVT_run2 = xlsread("max30102_data_Bryce_DVT_run2_02-18.xlsx");
data_bryce_DVT_run3 = xlsread("max30102_data_Bryce_DVT_run3_02-18.xlsx");

time_axis_Feb18 = zeros(size(data_hieu_noDVT_run1, 1), 1);

data_hieu_noDVT_run1_mod = data_hieu_noDVT_run1(:,2);
data_hieu_noDVT_run1_mod = data_hieu_noDVT_run1_mod(2:end);

data_hieu_noDVT_run2_mod = data_hieu_noDVT_run2(:,2);
data_hieu_noDVT_run2_mod = data_hieu_noDVT_run2_mod(2:end);

data_hieu_noDVT_run3_mod = data_hieu_noDVT_run3(:,2);
data_hieu_noDVT_run3_mod = data_hieu_noDVT_run3_mod(2:end);

data_hieu_DVT_run1_mod = data_hieu_DVT_run1(:,2);
data_hieu_DVT_run1_mod = data_hieu_DVT_run1_mod(2:end);

data_hieu_DVT_run2_mod = data_hieu_DVT_run2(:,2);
data_hieu_DVT_run2_mod = data_hieu_DVT_run2_mod(2:end);

data_hieu_DVT_run3_mod = data_hieu_DVT_run3(:,2);
data_hieu_DVT_run3_mod = data_hieu_DVT_run3_mod(2:end);

data_bryce_noDVT_run1_mod = data_bryce_noDVT_run1(:,2);
data_bryce_noDVT_run1_mod = data_bryce_noDVT_run1_mod(2:end);

data_bryce_noDVT_run2_mod = data_bryce_noDVT_run2(:,2);
data_bryce_noDVT_run2_mod = data_bryce_noDVT_run2_mod(2:end);

data_bryce_noDVT_run3_mod = data_bryce_noDVT_run3(:,2);
data_bryce_noDVT_run3_mod = data_bryce_noDVT_run3_mod(2:end);

data_bryce_DVT_run1_mod = data_bryce_DVT_run1(:,2);
data_bryce_DVT_run1_mod = data_bryce_DVT_run1_mod(2:end);

data_bryce_DVT_run2_mod = data_bryce_DVT_run2(:,2);
data_bryce_DVT_run2_mod = data_bryce_DVT_run2_mod(2:end);

data_bryce_DVT_run3_mod = data_bryce_DVT_run3(:,2);
data_bryce_DVT_run3_mod = data_bryce_DVT_run3_mod(2:end);

for n = 1:size(time_axis_Feb18, 1)
    time_axis_Feb18(n) = data_hieu_noDVT_run1(n, 1) - data_hieu_noDVT_run1(1, 1);
end

time_axis_Feb18 = time_axis_Feb18 .* 10^-3;
time_axis_Feb18 = time_axis_Feb18(2:end);

sampling_rate = 1/(0.04);

for n = 1:size(time_axis_Feb18, 1)
    time_axis_Feb18(n) = time_axis_Feb18(n, 1) - (1/sampling_rate);
end

close all

plot(time_axis_Feb18,data_hieu_noDVT_run1_mod);
hold on
plot(time_axis_Feb18,data_hieu_noDVT_run2_mod);
plot(time_axis_Feb18,data_hieu_noDVT_run3_mod);

plot(time_axis_Feb18,data_hieu_DVT_run1_mod);
plot(time_axis_Feb18,data_hieu_DVT_run2_mod);
plot(time_axis_Feb18,data_hieu_DVT_run3_mod);

xlabel("Time (s)")
ylabel("Amplitude")
title("Hieu DVT vs No DVT Feb 18")
legend("No DVT Trial 1", "No DVT Trial 2", "No DVT Trial 3", "DVT Trial 1", "DVT Trial 2", "DVT Trial 3")
hold off

freq_ax_Feb18 = ([1:size(time_axis_Feb18, 1)] - 1)/(size(time_axis_Feb18, 1)) .* sampling_rate;
freq_ax_Feb18 = freq_ax_Feb18';

data_hieu_noDVT_run1_fft = abs(fft(data_hieu_noDVT_run1_mod));
data_hieu_noDVT_run1_power = data_hieu_noDVT_run1_fft .* data_hieu_noDVT_run1_fft;

data_hieu_noDVT_run2_fft = abs(fft(data_hieu_noDVT_run2_mod));
data_hieu_noDVT_run2_power = data_hieu_noDVT_run2_fft .* data_hieu_noDVT_run2_fft;

data_hieu_noDVT_run3_fft = abs(fft(data_hieu_noDVT_run3_mod));
data_hieu_noDVT_run3_power = data_hieu_noDVT_run3_fft .* data_hieu_noDVT_run3_fft;

data_hieu_DVT_run1_fft = abs(fft(data_hieu_DVT_run1_mod));
data_hieu_DVT_run1_power = data_hieu_DVT_run1_fft .* data_hieu_DVT_run1_fft;

data_hieu_DVT_run2_fft = abs(fft(data_hieu_DVT_run2_mod));
data_hieu_DVT_run2_power = data_hieu_DVT_run2_fft .* data_hieu_DVT_run2_fft;

data_hieu_DVT_run3_fft = abs(fft(data_hieu_DVT_run3_mod));
data_hieu_DVT_run3_power = data_hieu_DVT_run3_fft .* data_hieu_DVT_run3_fft;

data_bryce_noDVT_run1_fft = abs(fft(data_bryce_noDVT_run1_mod));
data_bryce_noDVT_run1_power = data_bryce_noDVT_run1_fft .* data_bryce_noDVT_run1_fft;

data_bryce_noDVT_run2_fft = abs(fft(data_bryce_noDVT_run2_mod));
data_bryce_noDVT_run2_power = data_bryce_noDVT_run2_fft .* data_bryce_noDVT_run2_fft;

data_bryce_noDVT_run3_fft = abs(fft(data_bryce_noDVT_run3_mod));
data_bryce_noDVT_run3_power = data_bryce_noDVT_run3_fft .* data_bryce_noDVT_run3_fft;

data_bryce_DVT_run1_fft = abs(fft(data_bryce_DVT_run1_mod));
data_bryce_DVT_run1_power = data_bryce_DVT_run1_fft .* data_bryce_DVT_run1_fft;

data_bryce_DVT_run2_fft = abs(fft(data_bryce_DVT_run2_mod));
data_bryce_DVT_run2_power = data_bryce_DVT_run2_fft .* data_bryce_DVT_run2_fft;

data_bryce_DVT_run3_fft = abs(fft(data_bryce_DVT_run3_mod));
data_bryce_DVT_run3_power = data_bryce_DVT_run3_fft .* data_bryce_DVT_run3_fft;

figure
plot(freq_ax_Feb18(1:742), log(data_hieu_noDVT_run1_power(1:742)));
hold on
plot(freq_ax_Feb18(1:742), log(data_hieu_noDVT_run2_power(1:742)));
plot(freq_ax_Feb18(1:742), log(data_hieu_noDVT_run3_power(1:742)));

plot(freq_ax_Feb18(1:742), log(data_hieu_DVT_run1_power(1:742)));
plot(freq_ax_Feb18(1:742), log(data_hieu_DVT_run2_power(1:742)));
plot(freq_ax_Feb18(1:742), log(data_hieu_DVT_run3_power(1:742)));

xlabel("Freq (Hz)")
ylabel("Ln (amplitude)")
title("Hieu - Logarithmic Power Spectrum of DVT vs Non-DVT")
legend("No DVT 1", "No DVT 2", "No DVT 3", "DVT 1", "DVT 2", "DVT 3");
hold off

% figure
% plot(freq_ax_Feb18(1:742), log(data_bryce_noDVT_run1_power(1:742)));
% hold on
% plot(freq_ax_Feb18(1:742), log(data_bryce_noDVT_run2_power(1:742)));
% plot(freq_ax_Feb18(1:742), log(data_bryce_noDVT_run3_power(1:742)));
% 
% plot(freq_ax_Feb18(1:742), log(data_bryce_DVT_run1_power(1:742)));
% plot(freq_ax_Feb18(1:742), log(data_bryce_DVT_run2_power(1:742)));
% plot(freq_ax_Feb18(1:742), log(data_bryce_DVT_run3_power(1:742)));
% 
% xlabel("Freq (Hz)")
% ylabel("Log (amplitude)")
% title("Bryce - Logarithmic Power Sprectrum of DVT vs Non-DVT")
% legend("No DVT 1", "No DVT 2", "No DVT 3", "DVT 1", "DVT 2", "DVT 3");
% hold off