function ppu4fMRI
%%%%
%AIM:       This function aims at loading ScanPhysLog files containing ppu data,
%           computing the tachogram (HRV) and extracting and resampling the HF-HRV 
%           and LF-HRV to be used as regressors in fMRI (task or rs)analysis
%
%INPUT :    - a log file named 'SCANPHYSLOGxxx.log' 
%           acquired by the Philips 3T scanner
%           - acquisition parameters
%
%OUTPUT :   - HR and HRV parameters
%           - LF-HRV and HF-HRV regressors for fMRI
%
%HISTORY :
%			Written by C.Delon-Martin, the 16-01-2024.Modified by E. Gourieux 07-2025 
%
%ALGORITHM :
% 			1st, select the log file and enter acquisition parameters
%			2nd, load the file (duration relates to file length: can be 10 min!
%           3rd, plot the gradients' time courses and find begin and end of
%           sequence acquisition (visual control)
%           4th, compute the peaks of PPU for computing tachogram (visual control)
%           5th, compute global frequency analysis
%           6th, compute LF and HF components of tachogram (visual control)
%           7th, resample at TR the LF and HF components of tachogram
%           8th, save HR and HRV values and LF-regressors for fMRI (.tsv file) 
%
%EXTERNAL FUNCTIONS :
%			spm routines: spm_select
%			
%KEYWORDS : Physiological Data, frequential analysis

%% Default value
ndyn_def = '126'; % Number of dynamics
ndum_def = '4'; % Number of dummies before fMRI acquisition
TR_def = '2'; % Repetition time (s)

%% processing parameters (do NOT modify)
fsamp = 500;    % (Hz) sampling frequency at acquisition
fmax  = 10;     % (Hz) resampling frequency for tachogram
fdisp = 1;      % (Hz) max frequency for frequency spectrum display
fLFl  = 0.04;   % (Hz) LF lower limit
fLFh  = 0.15;   % (Hz) LF upper limit = HF lower limit
fHFh  = 0.4;    % (Hz) HF upper limit

%% display parameters
close all;
monitor = groot;
scz     = monitor.ScreenSize;% taille ecran
hfig    = round(scz(4)/5);
sz      = 50; %scatter size

%% LOADING DATA

% interactive file selection and paramaters
curr_dir= pwd;
R       = spm_select(1,'^SCANPHYSLOG.*\.log$','Select log file');
[logdir, logname] = fileparts(R);

defval  = {ndyn_def,ndum_def, TR_def};
answer  = inputdlg({'Number of dynamics','Number of dummies','TR (s)'},...
              'Enter acquisition parameters', [1 50; 1 50; 1 50], defval);
ndyn    = str2num(answer{1}); % nb of volumes(dynamics) during the functional sequence
ndum    = str2num(answer{2}); % nb of dummies before fMRI acquisition
TR      = str2num(answer{3}); % (s) repetition time
seq_dur = (ndyn+ndum)*TR; % (s) sequence duration

% load data from file
disp(['Analyze of file: ',logname,'.log']);
disp(' '); 
disp(' Loading file... be patient: it can take few minutes ');
cd(logdir);

flog    = fopen([logname,'.log'],'rt');
tic
for iline=1:6
   newline=fgetl(flog);
   disp(newline)
end
iline=0; %index pour les lignes du tableau
newline=fgetl(flog);
while strcmp(newline(1),'#');
      disp(newline)
end
while ischar(newline)
      iline=iline+1;
      data(iline,1:11)=str2num(newline);
      newline=fgetl(flog);
end
fclose(flog);
toc
disp('Data file loaded ');

l_scan  = length(data); % length of data
tot_dur = l_scan/fsamp;% (s) total duration of the ScanPhysLog

%% FROM GRADIENTS INFORMATION, detecting the begining and end of acquisition
grad    = data(:,7:9); %correspond to gradx, grady and gradz
gradtot = sum(abs(grad),2); %

% plot gradient
fppu    = figure('Name',['Analysis of ',logname],'NumberTitle','off','Position',[1 scz(4)/2 scz(3)-10 scz(4)/2]); 
agrad   = subplot(5,1,1); 
xt      = 1/fsamp:1/fsamp:length(gradtot)/fsamp; % xt time vector in s
plot(xt,gradtot,'k') % to visualize the gradients of the sequence
title('Gradients with sequence acquisition','FontSize',14)

% get the beginning and the end of the acquistion
indgrad = find((gradtot~=0) & circshift(gradtot(:),-1)==0 & circshift(gradtot(:),-2)==0 & circshift(gradtot(:),-3)==0 ...
    & circshift(gradtot(:),-4)==0 & circshift(gradtot(:),-5)==0 & circshift(gradtot(:),-6)==0 ...
    & circshift(gradtot(:),1)~=0 & circshift(gradtot(:),2)~=0 & circshift(gradtot(:),3)~=0);
lgrad   = size(indgrad,1); 
iend    = indgrad(lgrad) -1; %ok iend

indgrad = find((gradtot==0) & circshift(gradtot(:),1)==0 & circshift(gradtot(:),2)==0 & circshift(gradtot(:),3)==0 ...
   & circshift(gradtot(:),4)==0 & circshift(gradtot(:),5)==0 & circshift(gradtot(:),6)==0 ...
   & circshift(gradtot(:),-1)~=0 & circshift(gradtot(:),-2)~=0 & circshift(gradtot(:),-3)~=0);
lgrad   = size(indgrad,1);
if indgrad(lgrad) == length(gradtot)
    ibegin = indgrad(lgrad-1);%ok ibegin
else
    ibegin = indgrad(lgrad); %ok ibegin
end
disp(['Sequence acquisition from index ', num2str(ibegin), ' to index ', num2str(iend)]);
tbegin  = ibegin / fsamp; %start acquisition time (s) from begining of the file
tend    = iend / fsamp; %end of acquisition time (s) from begining of the file

% visualization to check the begin and end of sequence
hold on; scatter(tbegin,0,sz,'xr'); scatter(tend,0,sz,'xr');

%% PPU ANALYSIS
ppu     = data(:,5);
% extract peak
[ppus, peak_val,peak_ind] = extract_peaks(ppu);

% plot ppu and peaks
appu    = subplot(5,1,2);  
plot(xt,ppus,'k');
xl = xlim;
title('PPU with peak extractions','FontSize',14)
hold on; 
plot(appu, peak_ind/fsamp, peak_val,'xg');

% compute raw tachogram
R      = peak_ind/fsamp; % time of R-peaks (in seconds)
RR     = R - circshift(R,1); % RR time interval
RR(1)  = RR(2);
tR     = cumsum(RR); % in seconds

% QC 
[mean_RR, std_RR, cvar_RR] = report_QC(RR,'Tachogram before correction ');

% resample tachogram for having a regular grid (necessary for fft)
mRR = mean(RR);
ty     = 0:1/fmax:round(max(tR));     % resampled time interval at fmax
yrr    = interp1(tR,RR-mRR,ty,'pchip')'; % resampled tachogram (less interp error at beginning than spline)
yrr    = yrr + mRR;
ty     = ty';

% plot tachogram
atach   = subplot(5,1,3); 
plot(ty(round(tbegin*fmax):min(length(ty),round(tend*fmax))),yrr(round(tbegin*fmax):min(length(yrr),round(tend*fmax))),'k'); %regular tachogram
xlim(xl);
title('Tachogram and (LF+HF) components');

% correct tachogram (peak not detected)
[tR_corr,RR_corr] = tacho_correct(tR,RR);

% QC 
[mean_RR_corr, std_RR_corr, cvar_RR_corr] = report_QC(RR_corr,'Tachogram after correction ');

% resample corrected tachogram for having a regular grid (necessary for fft)
mRR_corr = mean(RR_corr);
ty       = 0:1/fmax:round(max(tR_corr));     % resampled time interval at fmax
yrr_corr = interp1(tR_corr,RR_corr-mRR_corr,ty,'pchip')'; % resampled tachogram (less interp error at beginning than spline)
yrr_corr = yrr_corr + mRR_corr;
ty       = ty';
yrr      = yrr_corr;
myrr     = mean(yrr);

% plot corrected tachogram
hold on; 
plot(atach,ty(round(tbegin*fmax):min(length(ty),round(tend*fmax))),yrr_corr(round(tbegin*fmax):min(length(yrr_corr),round(tend*fmax))),'Color','#008000');


%% Heart Rate computation
HR     = 60./RR_corr;           % heart rate raw (bpm)
meanHR = mean(HR);              % mean heart rate (bpm)
stdHR  = std(HR);               % std heart rate (bpm)
HRsamp = 60./yrr;               % heart rate resampled (bpm)


%% FREQUENTIAL ANALYSIS OF HRV - usinf FTT

kLFl    = ceil(length(ty)*fLFl/fmax); % index of LFlow
kLFh    = ceil(length(ty)*fLFh/fmax); % index of LFhigh
kHFh    = ceil(length(ty)*fHFh/fmax); % index of HFhigh
kdisp   = floor(length(ty)*fdisp/fmax);% index of max freq to display
fscan   = 0:fmax/length(ty):fmax-fmax/length(ty); % frequency vector
[ftppus,tot_pow, LF_pow, HF_pow, LFnu, HFnu,LF_HF_ratio,ftppunom] = hrv_freq_analysis(yrr,ty,fmax,kLFl,kLFh,kHFh);

% plot frequency spectrum
fspec    = figure('Name',['Frequency spectrum of ',logname],'NumberTitle','off','Position',[scz(4)*0.6 scz(3)*0.4 scz(3)*0.3 scz(4)*0.3]); 
plot(fscan(1:kdisp+1),ftppus(1:kdisp+1)/max(ftppus),'k');
title('Frequency spectrum of PPU','FontSize',14);
xlabel('frequency (Hz)','FontSize',14);
xLFl     = [fscan(kLFl) fscan(kLFl)]; y = [0 1]; line(xLFl,y,'Color','r');
xLFh     = [fscan(kLFh) fscan(kLFh)]; line(xLFh,y,'Color','r');
xHFh     = [fscan(kHFh) fscan(kHFh)]; line(xHFh,y,'Color','r');

%% TIME ANALYSIS of HRV

[SDNN, RMSSD] = hrv_time_analysis(RR_corr);

%% BANDPASS FILTERING (LF- and HF- components) of the tachogram

% LF-component and HF-component
% Filtering parameters : butterworth
Wl = [fLFl fLFh]; % frequency band (Hz) for LF component
Wh = [fLFh fHFh]; % frequency band (Hz) for HF component
[bl, al] = butter(4, Wl/(fmax/2), 'bandpass'); % Butterworth filter of order 4 for LF filter
[bh, ah] = butter(4, Wh/(fmax/2), 'bandpass'); % Butterworth filter of order 4 for HF filter

% Apply filters for LF and HF components
yLF = filtfilt(bl, al, yrr-myrr); 
yHF = filtfilt(bh, ah, yrr-myrr); 

% plot
figure(fppu)
acomp   = subplot(5,1,4);
plot(ty(round(tbegin*fmax):min(length(ty),round(tend*fmax))),yLF(round(tbegin*fmax):min(length(yLF),round(tend*fmax))),'Color','m'); %display LF component of tachogram 
title('LF components of tachogram');
xlim(xl);

acomp2   = subplot(5,1,5);
plot(ty(round(tbegin*fmax):min(length(ty),round(tend*fmax))),yHF(round(tbegin*fmax):min(length(yHF),round(tend*fmax))),'Color','m'); %display HF component 
title('HF components of tachogram');
xlim(xl);
xlabel('Time (s)');

% plot, comparison of LF + HF with tachogram (should be very close)
axes(atach);
hold on; 
plot(ty(round(tbegin*fmax):min(length(ty),round(tend*fmax))),yLF(round(tbegin*fmax):min(length(yLF),round(tend*fmax)))+yHF(round(tbegin*fmax):min(length(yHF),round(tend*fmax)))+mean_RR,'m'); %regular tachogram
legend('Tachogram','Tachogram corrected','LF + HF components','NumColumns',1,'Location','northwest');

% create LF- and HF regressors for fMRI (use tbegin and tend)
LF_reg  = yLF(round(tbegin*fmax):min(round(tend*fmax),length(yLF)));
HF_reg  = yHF(round(tbegin*fmax):min(round(tend*fmax),length(yHF)));
t_reg   = round(tbegin*fmax):round(tend*fmax);

% resample at TR and display
[yLF_reg,treg] = resample(LF_reg, round(tbegin*fmax):min(round(tend*fmax),length(yLF)), 1/(fmax*TR)); %resampling at TR
[yHF_reg,treg] = resample(HF_reg, round(tbegin*fmax):min(round(tend*fmax),length(yHF)), 1/(fmax*TR)); 

% plot regressor
axes(acomp);
hold on; plot(treg/fmax,yLF_reg,'Color','b');
legend('LF_component','LF_regressor','NumColumns',2,'Location','northwest')

axes(acomp2);
hold on; plot(treg/fmax,yHF_reg,'Color','b');
legend('HF_component','HF_regressor','NumColumns',2,'Location','northwest')

%% SAVE FIGURE FOR VISUAL INSPECTION
saveas(fppu,[logdir,filesep,logname,'_QC.png'])% for Quality Control
saveas(fspec,[logdir,filesep,logname,'_frequency_spectrum.png'])% for spectrum
    
%% CREATE OUTPUT FILES (.tsv) 
fout = fopen([logdir,filesep,logname,'_results.tsv'],'w');
fprintf(fout,['TACHOGRAM ANALYSIS RESULTS of ',logname,'\n\n']);
fprintf(fout,'Heart Rate\n');
fprintf(fout,'Heart rate (bpm) (mean HR) \t %8.4f\n',meanHR);
fprintf(fout,'sdtHR \t %8.4f\n\n',stdHR);
fprintf(fout,'Tachogram Quality Control (before correction) \n');
fprintf(fout,'mean_RR \t %8.4f\n',mean_RR);
fprintf(fout,'std_RR \t %8.4f\n',std_RR);
fprintf(fout,'cv_RR \t %8.4f\n\n',cvar_RR);
fprintf(fout,'Tachogram Quality Control (after correction) \n');
fprintf(fout,'mean_RR_corr \t %8.4f\n',mean_RR_corr);
fprintf(fout,'std_RR_corr \t %8.4f\n',std_RR_corr);
fprintf(fout,'cv_RR_corr \t %8.4f\n\n',cvar_RR_corr);
fprintf(fout,'Time Analysis of HRV \n');
fprintf(fout,'RMSSD \t %8.4f\n',RMSSD);
fprintf(fout,'SDNN \t %8.4f\n\n',SDNN);
fprintf(fout,['Frequency Analysis of HRV \n']);
fprintf(fout,'Max ppu frequency \t %8.4f\n', ftppunom);
fprintf(fout,'Total power \t %8.4f\n',tot_pow);
fprintf(fout,'LF power \t %8.4f\n',LF_pow);
fprintf(fout,'HF power \t %8.4f\n',HF_pow);
fprintf(fout,'LF nu \t %8.4f\n',LFnu);
fprintf(fout,'HF nu \t %8.4f\n',HFnu);
fprintf(fout,'LF/HF ratio \t %8.4f\n\n',LF_HF_ratio);
fclose(fout);

fout = fopen([logdir,filesep,logname,'_reg_LF_HRV.tsv'],'w');
fprintf(fout,'LF-HRV\n');
fprintf(fout,'%8.5f\n',yLF_reg);  % print all regressors' values
fclose(fout);


end

%%%%%%%%%%%%%%%
% SUB-FUNCTIONS 
%%%%%%%%%%%%%%%
function [ppus,peakval,peakind] = extract_peaks(ppu)

b       = [0 0.2 0.4 0.6 0.8 1 0.8 0.6 0.4 0.2 0]; %small smoothing kernel
pct     = 10; % percentage of max amplitude seems OK for peak finding (can be adjusted)
minpeakdist = 250; % min peak distance between maxima seems ok
minpeakwidth = 50; % mean peak width 

ppus    = conv(ppu,b,'same')/sum(b); % small smoothing to avoid multiple peaks. No shift in peaks

minppus = min(ppus); % find max amplitude of smoothed ppu
maxppus = max(ppus);
minPeakHeight = pct/100*(maxppus - minppus); % pct% of max amplitude 

[peakval,peakind] = findpeaks(ppus-minppus,'MinPeakHeight',minPeakHeight - minppus,'MinPeakWidth',minpeakwidth, 'MinPeakDistance', minpeakdist); % search for maxima
peakval = peakval + minppus;

end
%%%%%%%%
function [tR_corr,RR_corr] = tacho_correct(tR,RR)

% correct for missed peaks if coefficient of variation > 7.5% (empirical value at rest)
q      = quantile(RR,[0.25 0.5 0.75]);% compute quartiles
medRR  = q(2);                        % median of RR
iqr    = q(3) - q(1);                 % interquartile interval
valoutsup = find(RR > q(3)+3*iqr);    % find missing peaks
RR_corr = RR;                         % initialize RR_corr

if length(valoutsup)>=1 % exist missed peaks: ok
    disp('Correction for missed peaks');
    RR_corr = RR;
    for i_missed = 1:length(valoutsup)
        missed_npeak = round(RR_corr(valoutsup(i_missed))/medRR) - 1; % nb of missed peaks
        RR_int   = zeros(length(RR_corr)+missed_npeak,1);  % initialize corrected RR vector
        RR_int(1:valoutsup(i_missed)-1) = RR_corr(1:valoutsup(i_missed)-1); % copy correct RR part of vector
        RR_int(valoutsup(i_missed):valoutsup(i_missed)+missed_npeak) = RR_corr(valoutsup(i_missed))/(missed_npeak + 1);% add missed peak
        RR_int(valoutsup(i_missed)+missed_npeak+1:length(RR_int)) = RR_corr(valoutsup(i_missed)+1:length(RR_corr)); % copy remaining part of vector
        clear RR_corr; % since RR_corr and RR_int do not have the same length
        RR_corr = RR_int;
        valoutsup = valoutsup+missed_npeak; % because the missed peak is moved
    end
    % tR_corr = cumsum(RR_corr)+ t_init; 
    tR_corr = cumsum(RR_corr);
end

% correct for misplaced peaks : that is one peak overestimated followed by
% one peak under-estimated => replace by 2 equal intervals
q      = quantile(RR_corr,[0.25 0.5 0.75]);% compute quartiles after missed peaks correction
medRR  = q(2);                        % new median of RR
iqr    = q(3) - q(1);                 % new interquartile interval
valoutinf = find(RR_corr < q(1)-3*iqr); % find misplaced peaks

if length(valoutinf)>=1 % exist misplaced peaks
    disp('Correction for misplaced peaks');
    for i_missed = 1:length(valoutinf)
        RR_int2 = (RR_corr(valoutinf(i_missed)) + RR_corr(valoutinf(i_missed)-1))/2;
        RR_corr(valoutinf(i_missed)-1) = RR_int2;
        RR_corr(valoutinf(i_missed)) = RR_int2;
        tR_corr = cumsum(RR_corr);
    end
    tR_corr = cumsum(RR_corr); 
end
tR_corr = cumsum(RR_corr); % tR_corr repeated (in case no correction to apply)
end    
%%%%%%%%
function [meandata, stddata, cvardata] = report_QC(data,dataname)
    %report mean, sttdev and coeff of variation of data 
    meandata = mean(data);         % in seconds
    stddata  = std(data);          % in seconds
    cvardata     = stddata/meandata;       % coeff variation
    disp(' '); disp(['          ',dataname,' Quality Control ']);
    disp(['mean ',dataname,' = ',num2str(meandata),' +/- ', num2str(stddata),' s'])
    disp(['coeff of variation = ',num2str(cvardata*100)])
end
%%%%%%%%
%%%%%%%%
function [SDNN, RMSSD] = hrv_time_analysis(tacho);
    % SDNN = std of RR intervals (after correction)
    SDNN       = std(tacho);
    % RMSSD calculation
    diff_tacho = abs(diff(tacho)); %successive RR diffs 
    RMSSD      = sqrt(sum(diff_tacho.^2)/length(diff_tacho)); % relates to parasympathetic
end
%%%%%%%%
function [ftppus, tot_pow, LF_pow, HF_pow, LFnu, HFnu,LF_HF_ratio,ftppunom] = hrv_freq_analysis(tacho,time,fmax,kLFl,kLFh,kHFh);
    
    ntot    = length(time); %nb of points after resampling
    dty     = 1/fmax;    % time interval after resampling
    tot_dur = ntot*dty;  % total duration after resampling
    dfy     = fmax/ntot; % elementary frequency interval "  "
    fscan   = 0:dfy:fmax-dfy; % frequency vector
    b       = [0 0.2 0.4 0.6 0.8 1 0.8 0.6 0.4 0.2 0]; %small smoothing kernel
        
    % frequential analysis of tachogram`
    ftppu   = abs(fft(tacho-mean(tacho)));      % amplitude of the spectrum
    ftppus  = conv(ftppu,b,'same')/sum(b);  % small smoothing
    [ftppumax,iftppumax] = max(ftppus(kLFl:kHFh)); % val and index of peak of the spectrum
    ftppunom = (iftppumax+kLFl)*dfy;               % index of the peak of frequency
    disp(['Max ppu frequency at ',num2str(ftppunom), ' Hz']);
    
    % HRV parameters
    tot_pow = sum(ftppus(1:kHFh+1).^2); % total power
    LF_pow  = sum(ftppus(kLFl:kLFh).^2);% LF power
    HF_pow  = sum(ftppus(kLFh:kHFh).^2);% HF power
    LFnu    = LF_pow/tot_pow;
    HFnu    = HF_pow/tot_pow;
    LF_HF_ratio = LFnu/HFnu;
end
