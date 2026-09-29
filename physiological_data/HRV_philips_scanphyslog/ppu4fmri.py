import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, butter, filtfilt, resample
from scipy.interpolate import PchipInterpolator
from scipy.fft import fft
import tkinter as tk
from tkinter import filedialog, simpledialog

def extract_peaks(ppu, fsamp):
    # small smoothing kernel
    b = np.array([0, 0.2, 0.4, 0.6, 0.8, 1, 0.8, 0.6, 0.4, 0.2, 0])
    pct = 10
    meanpeakdist_time = 0.5
    minpeakwidth_time = 0.1
    minpeakdist = meanpeakdist_time * fsamp
    minpeakwidth = minpeakwidth_time * fsamp
    
    ppus = np.convolve(ppu, b, 'same') / np.sum(b)
    print("ppus", ppus)
    
    minppus = np.min(ppus)
    maxppus = np.max(ppus)
    minPeakHeight = (pct / 100) * (maxppus - minppus)
    
    # search for maxima
    peaks, properties = find_peaks(ppus - minppus, height=minPeakHeight, width=minpeakwidth, distance=minpeakdist)
    peakval = ppus[peaks]
    peakind = peaks
    
    return ppus, peakval, peakind

def tacho_correct(tR, RR):
    # correct for missed peaks if coefficient of variation > 7.5%
    q = np.percentile(RR, [25, 50, 75])
    medRR = q[1]
    iqr = q[2] - q[0]
    valoutsup = np.where(RR > q[2] + 3 * iqr)[0]
    
    RR_corr = RR.copy()
    
    if len(valoutsup) >= 1:
        print('Correction for missed peaks')
        for i_missed in range(len(valoutsup)):
            missed_npeak = int(np.round(RR_corr[valoutsup[i_missed]] / medRR)) - 1
            RR_int = np.zeros(len(RR_corr) + missed_npeak)
            RR_int[:valoutsup[i_missed]] = RR_corr[:valoutsup[i_missed]]
            RR_int[valoutsup[i_missed]:valoutsup[i_missed]+missed_npeak+1] = RR_corr[valoutsup[i_missed]] / (missed_npeak + 1)
            RR_int[valoutsup[i_missed]+missed_npeak+1:] = RR_corr[valoutsup[i_missed]+1:]
            RR_corr = RR_int.copy()
            valoutsup += missed_npeak
            
        tR_corr = np.cumsum(RR_corr)
        
    q = np.percentile(RR_corr, [25, 50, 75])
    iqr = q[2] - q[0]
    valoutinf = np.where(RR_corr < q[0] - 3 * iqr)[0]
    
    if len(valoutinf) >= 1:
        print('Correction for misplaced peaks')
        for i_missed in range(len(valoutinf)):
            RR_int2 = (RR_corr[valoutinf[i_missed]] + RR_corr[valoutinf[i_missed]-1]) / 2
            RR_corr[valoutinf[i_missed]-1] = RR_int2
            RR_corr[valoutinf[i_missed]] = RR_int2
            
    tR_corr = np.cumsum(RR_corr)
    return tR_corr, RR_corr

def report_QC(data, dataname):
    meandata = np.mean(data)
    stddata = np.std(data)
    cvardata = stddata / meandata
    print(f'\n          {dataname} Quality Control ')
    print(f'mean {dataname} = {meandata:.4f} +/- {stddata:.4f} s')
    print(f'coeff of variation = {cvardata * 100:.4f}')
    return meandata, stddata, cvardata

def hrv_time_analysis(tacho):
    SDNN = np.std(tacho, ddof=1)
    diff_tacho = np.abs(np.diff(tacho))
    RMSSD = np.sqrt(np.sum(diff_tacho**2) / len(diff_tacho))
    return SDNN, RMSSD

def hrv_freq_analysis(tacho, time_vec, fmax, kLFl, kLFh, kHFh):
    ntot = len(time_vec)
    dfy = fmax / ntot
    fscan = np.arange(0, fmax, dfy)
    b = np.array([0, 0.2, 0.4, 0.6, 0.8, 1, 0.8, 0.6, 0.4, 0.2, 0])
    
    ftppu = np.abs(fft(tacho - np.mean(tacho)))
    ftppus = np.convolve(ftppu, b, 'same') / np.sum(b)
    
    kLFl, kLFh, kHFh = int(kLFl), int(kLFh), int(kHFh)
    
    iftppumax = np.argmax(ftppus[kLFl:kHFh+1])
    ftppunom = (iftppumax + kLFl) * dfy
    print(f'Max ppu frequency at {ftppunom:.4f} Hz')
    
    tot_pow = np.sum(ftppus[:kHFh+1]**2)
    LF_pow = np.sum(ftppus[kLFl:kLFh]**2)
    HF_pow = np.sum(ftppus[kLFh:kHFh]**2)
    
    LFnu = LF_pow / tot_pow if tot_pow else 0
    HFnu = HF_pow / tot_pow if tot_pow else 0
    LF_HF_ratio = LFnu / HFnu if HFnu else 0
    
    return ftppus, tot_pow, LF_pow, HF_pow, LFnu, HFnu, LF_HF_ratio, ftppunom

def ppu4fMRI():
    # Default values
    ndyn_def = '126'
    ndum_def = '4'
    TR_def = '2'
    
    fsamp = 500
    fmax = 10
    fdisp = 1
    fLFl = 0.04
    fLFh = 0.15
    fHFh = 0.4
    
    # Hide Tkinter root
    root = tk.Tk()
    root.withdraw()
    
    logpaths = filedialog.askopenfilenames(
        title='Select log files',
        filetypes=[("Log files", "*.log")]
    )

    if not logpaths:
        print("No files selected.")
    
    logpaths = list(logpaths)

    print(f"Selected files: {logpaths}")

    ndyn = int(simpledialog.askstring("Input", "Number of dynamics", initialvalue=ndyn_def))
    ndum = int(simpledialog.askstring("Input", "Number of dummies", initialvalue=ndum_def))
    TR = float(simpledialog.askstring("Input", "TR (s)", initialvalue=TR_def))
        

    for logpath in logpaths:
        print(f"\n {logpath}")
        logdir = os.path.dirname(logpath)
        logname_ext = os.path.basename(logpath)
        logname, _ = os.path.splitext(logname_ext)
        
        print(f'Analyze of file: {logname}.log\n')
        print(' Loading file... be patient: it can take few minutes ')
        
        # Parse file skipping header
        data = []
        with open(logpath, 'r') as flog:
            for line in flog:
                if not line.startswith('#'):
                    try:
                        data.append([float(x) for x in line.split()[:11]])
                    except ValueError:
                        continue
        data = np.array(data)
        print('Data file loaded ')
        print(data)
        
        l_scan = len(data)
        tot_dur = l_scan / fsamp
        
        grad = data[:, 6:9]
        gradtot = np.sum(np.abs(grad), axis=1)
        
        fig_ppu, axes_ppu = plt.subplots(5, 1, figsize=(12, 12), num=f'Analysis of {logname}')
        xt = np.arange(1, len(gradtot) + 1) / fsamp
        
        axes_ppu[0].plot(xt, gradtot, 'k')
        axes_ppu[0].set_title('Gradients with sequence acquisition', fontsize=14)
            
        # Get the beginning and the end of the acquisition
        indgrad = np.where(
            (gradtot != 0) &
            (np.roll(gradtot, -1) == 0) &
            (np.roll(gradtot, -2) == 0) &
            (np.roll(gradtot, -3) == 0) &
            (np.roll(gradtot, -4) == 0) &
            (np.roll(gradtot, -5) == 0) &
            (np.roll(gradtot, -6) == 0) &
            (np.roll(gradtot, 1) != 0) &
            (np.roll(gradtot, 2) != 0) &
            (np.roll(gradtot, 3) != 0)
        )[0]

        lgrad = len(indgrad)
        iend = indgrad[lgrad - 1] - 1  # ok iend
        indgrad = np.where(
            (gradtot == 0) &
            (np.roll(gradtot, 1) == 0) &
            (np.roll(gradtot, 2) == 0) &
            (np.roll(gradtot, 3) == 0) &
            (np.roll(gradtot, 4) == 0) &
            (np.roll(gradtot, 5) == 0) &
            (np.roll(gradtot, 6) == 0) &
            # (np.roll(gradtot, 7) == 0) &
            # (np.roll(gradtot, 8) == 0) &
            # (np.roll(gradtot, 9) == 0) &
            # (np.roll(gradtot, 10) == 0) &
            # (np.roll(gradtot, 11) == 0) &
            # (np.roll(gradtot, 12) == 0) &
            # (np.roll(gradtot, 13) == 0) &
            # (np.roll(gradtot, 14) == 0) &
            # (np.roll(gradtot, 15) == 0) &
            # (np.roll(gradtot, 16) == 0) &
            # (np.roll(gradtot, 17) == 0) &
            # (np.roll(gradtot, 18) == 0) &
            # (np.roll(gradtot, 19) == 0) &
            # (np.roll(gradtot, 20) == 0) &
            (np.roll(gradtot, -1) != 0) &
            (np.roll(gradtot, -2) != 0) &
            (np.roll(gradtot, -3) != 0) #&
            # (np.roll(gradtot, -4) != 0) &
            # (np.roll(gradtot, -5) != 0) &
            # (np.roll(gradtot, -6) != 0)
        )[0]

        lgrad = len(indgrad)
        if indgrad[lgrad - 1] == len(gradtot) - 1:
            ibegin = indgrad[lgrad - 2]  # ok ibegin
        else:
            ibegin = indgrad[lgrad - 1]  # ok ibegin
            
        print(f'Sequence acquisition from index {ibegin} to index {iend}')
        tbegin = ibegin / fsamp
        tend = iend / fsamp
        
        axes_ppu[0].plot(tbegin, 0, 'xr', markersize=10)
        axes_ppu[0].plot(tend, 0, 'xr', markersize=10)
        
        ppu = data[:, 4]
        print("ppu", ppu)
        ppus, peak_val, peak_ind = extract_peaks(ppu, fsamp)
        
        axes_ppu[1].plot(xt, ppus, 'k')
        axes_ppu[1].set_title('PPU with peak extractions', fontsize=14)
        axes_ppu[1].plot(peak_ind / fsamp, peak_val, 'xg')
        plt.show()
        print(peak_ind)
        R = peak_ind / fsamp
        print("R", R)
        RR = R - np.roll(R, 1)
        RR[0] = RR[1]
        tR = np.cumsum(RR)
        
        mean_RR, std_RR, cvar_RR = report_QC(RR, 'Tachogram before correction ')
        
        mRR = np.mean(RR)
        ty = np.arange(0, np.round(np.max(tR)) + 1/fmax, 1/fmax)
        
        pchip = PchipInterpolator(tR, RR - mRR)
        yrr = pchip(ty) + mRR
        t_start_idx = int(np.round(tbegin * fmax))
        t_end_idx = int(np.min([len(ty), np.round(tend * fmax)]))
        axes_ppu[2].plot(ty[t_start_idx:t_end_idx], yrr[t_start_idx:t_end_idx], color="#90C1EE")
        
        # correct tachogram
        tR_corr, RR_corr = tacho_correct(tR, RR)
        mean_RR_corr, std_RR_corr, cvar_RR_corr = report_QC(RR_corr, 'Tachogram after correction ')
        
        mRR_corr = np.mean(RR_corr)
        ty = np.arange(0, np.round(np.max(tR_corr)) + 1/fmax, 1/fmax)
        pchip_corr = PchipInterpolator(tR_corr, RR_corr - mRR_corr)
        yrr_corr = pchip_corr(ty) + mRR_corr
        yrr = yrr_corr
        myrr = np.mean(yrr)
        
        
       
        axes_ppu[2].plot(ty[t_start_idx:t_end_idx], yrr_corr[t_start_idx:t_end_idx], color='#008000')
        axes_ppu[2].set_title('Tachogram and (LF+HF) components')
        
        HR = 60. / RR_corr
        meanHR = np.mean(HR)
        stdHR = np.std(HR, ddof=1)
        HRsamp = 60. / yrr
        
        kLFl = np.ceil(len(ty) * fLFl / fmax)
        kLFh = np.ceil(len(ty) * fLFh / fmax)
        kHFh = np.ceil(len(ty) * fHFh / fmax)
        kdisp = int(np.floor(len(ty) * fdisp / fmax))
        fscan = np.arange(0, fmax, fmax / len(ty))
        
        ftppus, tot_pow, LF_pow, HF_pow, LFnu, HFnu, LF_HF_ratio, ftppunom = hrv_freq_analysis(yrr, ty, fmax, kLFl, kLFh, kHFh)
        
        fig_spec, axes_spec = plt.subplots(1, 1, figsize=(8, 6), num=f'Frequency spectrum of {logname}')
        axes_spec.plot(fscan[:kdisp+1], ftppus[:kdisp+1] / np.max(ftppus), 'k')
        axes_spec.set_title('Frequency spectrum of PPU', fontsize=14)
        axes_spec.set_xlabel('frequency (Hz)', fontsize=14)
        axes_spec.axvline(x=fscan[int(kLFl)], color='r')
        axes_spec.axvline(x=fscan[int(kLFh)], color='r')
        axes_spec.axvline(x=fscan[int(kHFh)], color='r')
        
        SDNN, RMSSD = hrv_time_analysis(RR_corr)
        
        bl, al = butter(4, [fLFl/(fmax/2), fLFh/(fmax/2)], 'bandpass')
        bh, ah = butter(4, [fLFh/(fmax/2), fHFh/(fmax/2)], 'bandpass')
        
        yLF = filtfilt(bl, al, yrr - myrr)
        yHF = filtfilt(bh, ah, yrr - myrr)
        
        axes_ppu[3].plot(ty[t_start_idx:t_end_idx], yLF[t_start_idx:t_end_idx], color='m')
        axes_ppu[3].set_title('LF components of tachogram')
        
        axes_ppu[4].plot(ty[t_start_idx:t_end_idx], yHF[t_start_idx:t_end_idx], color='m')
        axes_ppu[4].set_title('HF components of tachogram')
        axes_ppu[4].set_xlabel('Time (s)')
        
        axes_ppu[2].plot(ty[t_start_idx:t_end_idx], yLF[t_start_idx:t_end_idx] + yHF[t_start_idx:t_end_idx] + mean_RR, 'm')
        axes_ppu[2].legend(['Tachogram ', 'Tachogram corrected', 'LF + HF components'])
        
        LF_reg = yLF[t_start_idx:t_end_idx]
        HF_reg = yHF[t_start_idx:t_end_idx]
        
        num_resample_pts = int((t_end_idx - t_start_idx) / (fmax * TR))
        yLF_reg = resample(LF_reg, num_resample_pts)
        yHF_reg = resample(HF_reg, num_resample_pts)
        treg = np.linspace(tbegin, tend, num_resample_pts)
        
        axes_ppu[3].plot(treg, yLF_reg, color='b')
        axes_ppu[3].legend(['LF_component', 'LF_regressor'])
        
        axes_ppu[4].plot(treg, yHF_reg, color='b')
        axes_ppu[4].legend(['HF_component', 'HF_regressor'])
        
        fig_ppu.tight_layout()
        fig_ppu.savefig(os.path.join(logdir, f'{logname}_QC.png'))
        fig_spec.tight_layout()
        fig_spec.savefig(os.path.join(logdir, f'{logname}_frequency_spectrum.png'))
        
        with open(os.path.join(logdir, f'{logname}_results.tsv'), 'w') as fout:
            fout.write(f'TACHOGRAM ANALYSIS RESULTS of {logname}\n\n')
            fout.write('Heart Rate\n')
            fout.write(f'Heart rate (bpm) (mean HR)\t{meanHR:.4f}\n')
            fout.write(f'sdtHR\t{stdHR:.4f}\n\n')
            fout.write('Tachogram Quality Control (before correction)\n')
            fout.write(f'mean_RR\t{mean_RR:.4f}\n')
            fout.write(f'std_RR\t{std_RR:.4f}\n')
            fout.write(f'cv_RR\t{cvar_RR:.4f}\n\n')
            fout.write('Tachogram Quality Control (after correction)\n')
            fout.write(f'mean_RR_corr\t{mean_RR_corr:.4f}\n')
            fout.write(f'std_RR_corr\t{std_RR_corr:.4f}\n')
            fout.write(f'cv_RR_corr\t{cvar_RR_corr:.4f}\n\n')
            fout.write('Time Analysis of HRV\n')
            fout.write(f'RMSSD\t{RMSSD:.4f}\n')
            fout.write(f'SDNN\t{SDNN:.4f}\n\n')
            fout.write('Frequency Analysis of HRV\n')
            fout.write(f'Max ppu frequency\t{ftppunom:.4f}\n')
            fout.write(f'Total power\t{tot_pow:.4f}\n')
            fout.write(f'LF power\t{LF_pow:.4f}\n')
            fout.write(f'HF power\t{HF_pow:.4f}\n')
            fout.write(f'LF nu\t{LFnu:.4f}\n')
            fout.write(f'HF nu\t{HFnu:.4f}\n')
            fout.write(f'LF/HF ratio\t{LF_HF_ratio:.4f}\n\n')
            
        with open(os.path.join(logdir, f'{logname}_reg_LF_HRV_with-dummy.txt'), 'w') as fout:
            #fout.write('LF-HRV\n')
            for v in yLF_reg:
                fout.write(f'{v:.5f}\n')
        
        with open(os.path.join(logdir, f'{logname}_reg_LF_HRV.txt'), 'w') as fout:
            #fout.write('LF-HRV\n')
            yLF_reg_without_dum = yHF_reg[ndum :]
            for v in yLF_reg_without_dum:
                fout.write(f'{v:.5f}\n')
                
        plt.show()

if __name__ == '__main__':
    ppu4fMRI()