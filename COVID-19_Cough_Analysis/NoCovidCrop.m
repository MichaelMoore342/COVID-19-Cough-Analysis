%Michael Moore 2118213
%Christopher Rawlings 2179595


%This is the code used to separate the non-Covid-19 coughs from the original
%signal. Works with the noCovid folder

clc
sAudioFiles = dir(fullfile('C:\Program Files (x86)\University Year 3\Semester 2\Biomedical Signals, Systems and Control\Assignment\Project\noCovid'));  

n_fs = 48000;

for iCount = 3:13
    
    [no_covid, n_fs] = audioread(sAudioFiles(iCount).name);
    N = length(no_covid);
    dt = 1/n_fs;
    t = (0:N-1)/n_fs;
    
    figure(iCount +10); 
    no_covid =no_covid(:,1);
    plot(t, no_covid);
    xlabel('Seconds');
    ylabel('Amplitude');
 
    sTemp = imdilate(abs(no_covid), true(1000,1));
    inTemp = sTemp < 0.05;
    inTemp2 = inTemp(:).';
    A = strfind(inTemp2, [1 0]);
    B = strfind(inTemp2, [0 1]);
    
    iLength = length(A);
    
    C = B-A;
    
    D = C > 8000; 
    
    figure(20+iCount);
    iCounter = 0;
    for iCount2 = 1:iLength
        if D(1,iCount2) == 1
                samples = [A(iCount2),B(iCount2)];
                [no_covidCrop,n_fs] = audioread(sAudioFiles(iCount).name,samples);
                t_crop = (A(iCount2):B(iCount2))/n_fs;
                Filename = sprintf('no_covid%i_%i_crop.wav',iCount,iCount2-iCounter);
                audiowrite(Filename,no_covidCrop,n_fs); 

                
                subplot(2,2,iCount2-iCounter);
                plot(t_crop,no_covidCrop)
                titlename = sprintf('Covid Negative Patient %i', iCount -2);
                title(titlename);
                xlabel('Time (s)');
                ylabel('Amplitude');
        else
            iCounter = iCounter +1;
            
        end
               
    end      
 
   
end