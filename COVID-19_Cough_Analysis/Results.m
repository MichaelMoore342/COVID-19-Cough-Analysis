%Michael Moore 2118213
%Christopher Rawlings 2179595

%This code works using the noCovidCrop file
%To test a signal, enter its name in the line below and run the code

[covid1 , fs] = audioread('zcovid10_2_crop.wav'); %Read in audio file
covid1 = covid1(:,1);
N = length(covid1);
t = (0:N-1)/fs; %Convert sample number to time
duration = t(end)

transform_fft = fft(covid1); %Transform signal to frequency domain
power = abs(transform_fft).^2/N; % power of the DFT
f_range = (0:N-1)*(fs/N);
max_power = 0;
total_power = 0;

%Determine max power
for i = 1:length(power)
  if power(i)> max_power
      max_power = power(i);
  end
end

%Determine average power
for j = 1:length(power)
    total_power = total_power + power(j);
end

average_power = total_power/length(power)
    
max_power

rand_signal = covid1; %Random signal
total = 0;
fo = pitch(rand_signal,fs); %Determines pitch change over course of signal

% Sum up pitch values
for k = 1:length(fo)
    total = total + fo(k);
end

% Determine average pitch
average_pitch = total/length(fo)

nfft = 1024;
Trans = linspace(0,fs,nfft);
Y_n = abs(fft(covid1,nfft));

signal_F = Y_n;
amp = abs(signal_F(1,1));

for i = 1:length(signal_F)
    if abs(signal_F(i)) > amp
        amp = abs(signal_F(i));
    end
end

%Determines the actual frequncy values that have the largest amplitudes
%Two frequencies exist with the largest value, one at a low frequency and
%another at a large frequency
for j = 1:length(signal_F)
    if abs(signal_F(j)) == amp && j<(length(signal_F)/2)
        F_index_lower = j;
        
    elseif abs(signal_F(j)) == amp && j>(length(signal_F)/2) 
        F_index_greater = j; 
    end
end

%Display frequencies
lower_frequency = Trans(1,F_index_lower)
%Trans(1,F_index_greater)

signal = covid1;
max_peak = abs(signal(1,1));

for i = 1:length(signal)
    if abs(signal(i)) > max_peak
        max_peak = abs(signal(i));
    end
end


fraction_peak = max_peak/10;
peak_count = 0;

%Determines number of peaks above 10% of max peak height
for j = 1:length(signal)
    if abs(signal(j)) >= fraction_peak
        peak_count = peak_count + 1;
    end
end

%Determines relative percentage of peaks that are above 10% within the
%signal
peak_percentage = (peak_count/length(signal))*100

%Cough check
check_counter = 0;

%The next section of code detemines the absolute difference between
%calculated value and the known averages and thus determines if the signal is
%more similar to the COVID negative or positive data.
no_covid_power = abs(average_power - 0.02813);
covid_power = abs(average_power - 0.0343);

if no_covid_power < covid_power
    check_counter = check_counter + 1;
end

no_covid_f = abs(lower_frequency - 698.43);
covid_f = abs(lower_frequency - 726.74);

if no_covid_f < covid_f
    check_counter = check_counter + 1;
end

no_covid_duration = abs(duration - 0.301);
covid_duration = abs(duration - 0.3085);

if no_covid_duration < covid_duration
    check_counter = check_counter + 1;
end

no_covid_mp = abs(max_power - 8.324);
covid_mp = abs(max_power - 11.7625);

if no_covid_mp < covid_mp
    check_counter = check_counter + 1;
end

no_covid_peak = abs(peak_percentage - 46.98);
covid_peak = abs(peak_percentage - 42.8184);

if no_covid_peak < covid_peak
    check_counter = check_counter + 1;
end

no_covid_pitch = abs(average_pitch - 243.42);
covid_pitch = abs(average_pitch - 229.6578);

if no_covid_pitch < covid_pitch
    check_counter = check_counter + 1;
end

if max_power > 24
    check_counter = 0;
end

if max_power > 8 && max_power < 16
    check_counter = 6
end

if average_power > 0.06
    check_counter = 0;
end

if duration > 0.8
    check_counter = 0;
end

if peak_percentage > 30 && peak_percentage < 50
    check_counter = 0;
end

if peak_percentage > 50 && peak_percentage < 60
    check_counter = 6;
end

if peak_percentage > 20 && peak_percentage < 30
    check_counter = 6;
end

if check_counter == 1
    disp("Patient does not have COVID_19")
else
    disp("Patient does have COVID_19")
end

