%Code was developed by  Michael Moore (2118213) and Christopher Rawlings (2179595) 
%The noCovidCrop file is used for this code and has been as indexed specifically for this code.

clc
n_fs = 48000;
nfft = 1024;

sAudioFiles = dir(fullfile('C:\Program Files (x86)\University Year 3\Semester 2\Biomedical Signals, Systems and Control\Assignment\Project\noCovidCrop'));  

arrAudio = cell(1, 70);
arrAudioNCF = cell(1, 11);
arrAudioCF = cell(1, 9);
arrDuration = cell(1, 51);
arrFreqMax = cell(1, 51);
arrFreqMin = cell(1, 51);
arrMaxPeak = cell(1, 51);
arrNumPeaks = cell(1, 51);
arrMaxPeakPercentage = cell(1, 51);
arrAvgPitch = cell(1, 51);
arrMaxValue = cell(1, 51);
arrMinValue = cell(1, 51);
arrMeanValue = cell(1, 51);
arrSTDValue = cell(1, 51);
arrPowerAvg = cell(1, 51);
arrPowerMax = cell(1, 51);


for iCount = 3:53
   
   %Calculate the duration of all of the non-covid files.
    arrAudio{iCount-2} = audioread(sAudioFiles(iCount).name);  
    info = audioinfo(sAudioFiles(iCount).name);
    arrDuration{iCount-2} = [info.Duration];
    
 %Calculate the frequency of the non-covid files.
    [no_covid,n_fs] = audioread(sAudioFiles(iCount).name);
    Trans = linspace(0,n_fs,nfft);
    Y_n = abs(fft(no_covid,nfft));
    signal_F = Y_n;
    amp = abs(signal_F(1,1));
    dt = 1/n_fs;
    t = 0:dt:(length(no_covid)*dt)-dt;
    %Increasing or decreasing the iCount variable allows for graph
    %selection
    if iCount == 3
        figure(1)
         if iCount < 33
             
             iCountSO = floor(iCount/3);
             [no_covidOS, n_fs] = audioread(sAudioFiles(53+iCountSO).name);
             tNCOS = 0:dt:(length(no_covidOS)*dt)-dt;
             Trans = linspace(0,n_fs,nfft);
             Y_nOS = abs(fft(no_covidOS,nfft));
             
             subplot(2,1,1)
             plot(Trans(1:nfft),Y_nOS(1:nfft));
             xlabel('Frequency');
             ylabel('Absolute Value'); 
             title('Original signal');
             
             subplot(2,1,2)
             plot(Trans(1:nfft),Y_n(1:nfft));
             xlabel('Frequency');
             ylabel('Absolute Value');
             title('Extracted cough')
             
             subplot(2,1,1)
             plot(tNCOS, no_covidOS)
             xlabel('Seconds')
             ylabel('Amplitude')
             title('Original Signal')
             
             subplot(2,1,2)
             plot(t, no_covid)
             xlabel('Seconds');
             ylabel('Amplitude')
             title('Extracted cough')
             
             
            
         elseif iCount >= 33
             
             if (iCount >= 33) && (iCount < 36)
                iCountSO = 1;
             elseif (iCount >= 36) && (iCount < 40)
                iCountSO = 2;
             elseif (iCount == 40)
                iCountSO = 3;
             elseif (iCount >= 41) && (iCount < 44)
                iCountSO = 4;
             elseif (iCount >= 44) && (iCount < 47 )
                iCountSO = 5;
             elseif (iCount >= 47) && (iCount < 51)
                iCountSO = 6;
             elseif (iCount == 51) 
                iCountSO = 7; 
             elseif (iCount >= 52) && (iCount < 54)
                iCountSO = 8;
             end   
             
             [covidOS,n_fs] = audioread(sAudioFiles(64 +iCountSO).name);
              tCOS = 0:dt:(length(covidOS)*dt)-dt;
              Trans = linspace(0,n_fs,nfft);
              Y_nOS = abs(fft(covidOS,nfft));
              
             subplot(2,1,1)
             plot(tCOS, covidOS)
             xlabel('Seconds')
             ylabel('Amplitude')
             title('Original Signal')
             
             subplot(2,1,2)
             plot(t, no_covid)
             xlabel('Seconds');
             ylabel('Amplitude')
             title('Extracted cough')
             
             subplot(2,1,1)
             plot(Trans(1:nfft),Y_nOS(1:nfft));
             xlabel('Frequency');
             ylabel('Absolute Value'); 
             title('Original signal');
             
             subplot(2,1,2)
             plot(Trans(1:nfft),Y_n(1:nfft));
             xlabel('Frequency');
             ylabel('Absolute Value');
             title('Extracted cough')
             
         end
        
    end
    
    
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

    %Record frequencies
    arrFreqMin{iCount-2} =  F_index_lower;

    arrFreqMax{iCount-2} = F_index_greater;
    
    
    
 %Calculate the number of peaks
    signal = no_covid;
    max_peak = max(abs(signal));
    arrMaxPeak{1, iCount -2} = max_peak; 


    fraction_peak = max_peak/10;
    peak_count = 0;

    %Determines number of peaks above 10% of max peak height
        for j = 1:length(signal)
            if abs(signal(j)) >= fraction_peak
             peak_count = peak_count + 1;
            end
        end

     arrNumPeaks{1, iCount -2} = peak_count;   
        
    %Determines relative percentage of peaks that are above 10% within the
    %signal
    peak_percentage = (peak_count/length(signal))*100;
    arrMaxPeakPercentage{1, iCount -2} = peak_percentage;
    
 %Calculate the pitch
         rand_signal = no_covid;
         total = 0;
         fo = pitch(rand_signal,n_fs); %Determines pitch change over course of signal

        % Sum up pitch values
        for k = 1:length(fo)
            total = total + fo(k);
        end

        % Determine average pitch
        average_pitch = total/length(fo);
        arrAvgPitch{iCount-2} = average_pitch;
        
        
 %Calculate some basic statistics
 
    %calculate the maximum value
    max_Value = max(abs(signal));
    arrMaxValue{1, iCount -2} = max_Value;
    
    %Calculate the min Value
    min_Value = min((signal));
    arrMinValue{1, iCount -2} = min_Value;
    
    %Calculate the mean
    mean_Value = mean(abs(signal));
    arrMeanValue{1, iCount -2} = mean_Value;
    
    %Calculate SD value 
    std_Value = std(abs(signal));
    arrSTDValue{1, iCount -2} = std_Value; 
    
    
 %Calculate the power
    transform_fft = fft(signal); %Transform signal to frequency domain
    N_crop = length(signal);
    power = abs(transform_fft).^2/N_crop; % power of the DFT
    f_range = (0:N_crop-1)*(n_fs/N_crop);
    %plot(f_range,power);
    max_power = 0;
    total_power = 0;

    %Determine max power
        for i = 1:length(power)
            if power(i)> max_power
                 max_power = power(i);
            end
        end
        
        arrPowerMax{1, iCount-2} = max_power;
        
    %Determine average power
        for j = 1:length(power)
         total_power = total_power + power(j);
        end
        
        arrPowerAvg{1, iCount-2} = total_power/length(power);   
end   

%Plot the duration histogram
figure(2);
histogram(cell2mat(arrDuration(1:30)),(0:0.2:1),'Normalization','Probability', 'Facecolor', [0 0.4470 0.7410], 'FaceAlpha', 0.75)
hold on
histogram(cell2mat(arrDuration(31:50)),(0:0.2:1),'Normalization','Probability', 'Facecolor', 'r', 'FaceAlpha', 0.75)
title('Duration of non-covid patients')
xlabel('Duration (s)')
ylabel('Fraction')
legend({'Non-Covid','Covid'})
hold off

%Plot the minimum frequency histogram
figure(3);
subplot(1,2,1)
histogram(cell2mat(arrFreqMin(1:30)),(0:5:50),'Normalization','Probability',  'Facecolor', [0 0.4470 0.7410], 'FaceAlpha', 0.75)
hold on
histogram(cell2mat(arrFreqMin(31:51)),(0:5:50),'Normalization','Probability', 'Facecolor', 'r', 'FaceAlpha', 0.75)
title('The minimum frequencies')
xlabel('Frequency')
ylabel('Fraction')
legend({'Non-Covid','Covid'})
hold off

%Plot the maximum frequency histogram
subplot(1,2,2)
histogram(cell2mat(arrFreqMax(1:30)),(975:10:1035) , 'Normalization','Probability', 'Facecolor', [0 0.4470 0.7410], 'FaceAlpha', 0.75)
hold on
histogram(cell2mat(arrFreqMax(31:51)),(975:10:1035) , 'Normalization','Probability', 'Facecolor', 'r', 'FaceAlpha', 0.75)
title('The maximum frequencies')
xlabel('Frequency')
ylabel('Fraction')
legend({'Non-Covid','Covid'})
hold off

%Plot the highest peak histogram
figure(4);
histogram(cell2mat(arrMaxPeak(1:30)),(0:0.2:1) , 'Normalization','Probability', 'Facecolor', [0 0.4470 0.7410], 'FaceAlpha', 0.75)
hold on
histogram(cell2mat(arrMaxPeak(31:51)),(0:0.2:1) , 'Normalization','Probability', 'Facecolor', 'r', 'FaceAlpha', 0.75)
title('Value of the highest peak')
xlabel('Ranges of the values')
ylabel('Fraction')
legend({'Non-Covid','Covid'})
hold off

%Plot the number of peaks histogram
figure(5)
histogram(cell2mat(arrNumPeaks(1:30)), (0:2000:20000) , 'Normalization','Probability', 'Facecolor', [0 0.4470 0.7410], 'FaceAlpha', 0.75)
hold on
histogram(cell2mat(arrNumPeaks(31:51)), (0:2000:20000) , 'Normalization','Probability', 'Facecolor', 'r', 'FaceAlpha', 0.75)
title('The number of peaks above 10% of max peak height')
xlabel('Number of peaks')
ylabel('Fraction')
legend({'Non-Covid','Covid'})
hold off

%Plot the percentage of peaks histogram
figure(6)
histogram(cell2mat(arrMaxPeakPercentage(1:30)), (0:10:100) , 'Normalization','Probability', 'Facecolor', [0 0.4470 0.7410], 'FaceAlpha', 0.75)
hold on
histogram(cell2mat(arrMaxPeakPercentage(31:51)), (0:10:100) , 'Normalization','Probability', 'Facecolor', 'r', 'FaceAlpha', 0.75)
title('The percentage of peaks above 10% of max peak height')
xlabel('Percentage of peaks (%)')
ylabel('Fraction')
legend({'Non-Covid','Covid'})
hold off

%Plot the average pitch histogram
figure(7)
histogram(cell2mat(arrAvgPitch(1:30)), (50:50:400), 'Normalization','Probability', 'Facecolor', [0 0.4470 0.7410], 'FaceAlpha', 0.75)
hold on
histogram(cell2mat(arrAvgPitch(31:51)) , (50:50:400) , 'Normalization','Probability',  'Facecolor', 'r', 'FaceAlpha', 0.75)
title('The average pitch')
xlabel('Pitch')
ylabel('Fractione')
legend({'Non-Covid','Covid'})
hold off

%Plot the maximum value histogram
figure(8)
histogram(cell2mat(arrMaxValue(1:30)),(0:0.2:1.2), 'Normalization','Probability', 'Facecolor', [0 0.4470 0.7410], 'FaceAlpha', 0.75)
hold on
histogram(cell2mat(arrMaxValue(31:51)) ,(0:0.2:1.2), 'Normalization','Probability',  'Facecolor', 'r', 'FaceAlpha', 0.75)
title('The Maximum Value')
xlabel('Value')
ylabel('Percentage')
legend({'Non-Covid','Covid'})
hold off

%Plot the minimum value histogram
figure(9)
histogram(cell2mat(arrMinValue(1:30)),(-1.4:0.2:0), 'Normalization','Probability', 'Facecolor', [0 0.4470 0.7410], 'FaceAlpha', 0.75)
hold on
histogram(cell2mat(arrMinValue(31:51)),(-1.4:0.2:0) , 'Normalization','Probability',  'Facecolor', 'r', 'FaceAlpha', 0.75)
title('The Minimum Value')
xlabel('Value')
ylabel('Percentage')
legend({'Non-Covid','Covid'})
hold off

%Plot the mean value histogram
figure(10)
histogram(cell2mat(arrMeanValue(1:30)),(0:0.02:0.2), 'Normalization','Probability', 'Facecolor', [0 0.4470 0.7410], 'FaceAlpha', 0.75)
hold on
histogram(cell2mat(arrMeanValue(31:51)) ,(0:0.02:0.2), 'Normalization','Probability',  'Facecolor', 'r', 'FaceAlpha', 0.75)
title('The Mean Value')
xlabel('Value')
ylabel('Percentage')
legend({'Non-Covid','Covid'})
ylim([0 0.5]);
hold off

%Plot the standard deviation histogram
figure(11)
histogram(cell2mat(arrSTDValue(1:30)),(0:0.05:0.3), 'Normalization','Probability', 'Facecolor', [0 0.4470 0.7410], 'FaceAlpha', 0.75)
hold on
histogram(cell2mat(arrSTDValue(31:51)),(0:0.05:0.3), 'Normalization','Probability',  'Facecolor', 'r', 'FaceAlpha', 0.75)
title('The STD Value')
xlabel('Value')
ylabel('Percentage')
legend({'Non-Covid','Covid'})
ylim([0 0.5]);
hold off

%Plot the average power histogram
figure(12)
subplot(1,2,1)
histogram(cell2mat(arrPowerAvg(1:30)),(0:0.03:0.15), 'Normalization','Probability', 'Facecolor', [0 0.4470 0.7410], 'FaceAlpha', 0.75)
hold on
histogram(cell2mat(arrPowerAvg(31:51)),(0:0.03:0.15), 'Normalization','Probability',  'Facecolor', 'r', 'FaceAlpha', 0.75)
title('The Average Power')
xlabel('Value')
ylabel('Percentage')
legend({'Non-Covid','Covid'})
hold off

%Plot the maximum power histogram
subplot(1,2,2)
histogram(cell2mat(arrPowerMax(1:30)),(0:8:40), 'Normalization','Probability', 'Facecolor', [0 0.4470 0.7410], 'FaceAlpha', 0.75)
hold on
histogram(cell2mat(arrPowerMax(31:51)),(0:8:40), 'Normalization','Probability',  'Facecolor', 'r', 'FaceAlpha', 0.75)
title('The Maximum Power')
xlabel('Value')
ylabel('Percentage')
legend({'Non-Covid','Covid'})
hold off

Avgpower = cell2table(arrPowerAvg)
DFreq = cell2table(arrFreqMin)
Duration = cell2table(arrDuration)
MPower = cell2table(arrPowerMax)
PPeaks = cell2table(arrMaxPeakPercentage)
Pitch = cell2table(arrAvgPitch)

MaxMag = cell2table(arrMaxValue)
NPeaks = cell2table(arrNumPeaks)
Mean = cell2table(arrMeanValue)
SD = cell2table(arrSTDValue)
MinFreq = cell2table(arrFreqMin)
MaxFreq = cell2table(arrFreqMax)