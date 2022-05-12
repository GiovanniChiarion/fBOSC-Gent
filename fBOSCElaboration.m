% Add paths
addpath("C:\Users\Giovanni\OneDrive - Politecnico di Torino\PhD\Lavoro\Addons\FieldTrip")
addpath(genpath("fBOSC-main"))

%%

% clear
freq = 40; % Hz
search_precision = 1; % Hz

folder = 'newdata';
load(fullfile(folder,'data_test_ASSR_M1_wideband.mat'))

fs = data_test.fsample;
times = data_test.time;
trials = data_test.trial;

plotFactor = std(data_test.trial{1})*10;


%% Visualization

figure
for i=1:10
plot(data_test.time{i},data_test.trial{i}+i*plotFactor,'k');
hold on 
end
ylim tight
xline(0,'r','LineWidth',2)
title('First 10 trials')
xlabel('seconds')

%% Data preprocessing

data = ft_preprocessing([],data_test);

% cfg = [];
% cfg.resamplefs = 1000;
% data = ft_resampledata(cfg,data);
% fs=1000;

% cfg = [];
% cfg.hpfilter = "yes";
% cfg.hpfreq = 1;
% cfg.hpfiltord = 3;
% data = ft_preprocessing(cfg,data);

% Selection of part of data
% cfg = [];
% cfg.toilim = [0,1];
% data = ft_redefinetrial(cfg,data);

% Obtaining one averaged trial
% new_data = ft_timelockanalysis([],data);
% figure
% plot(new_data.time,new_data.avg)

%% Set-up fBOSC parameters

BOSC_type = 'fBOSC'; % or eBOSC

start_fBOSC;

% general setup
cfg = [];
cfg.(BOSC_type).F                 = [39:.1:41];
cfg.(BOSC_type).wavenumber        = 6;           % wavelet family parameter (time-frequency tradeoff)
cfg.(BOSC_type).fsample           = fs;         % current sampling frequency of EEG data

% padding
cfg.(BOSC_type).pad.tfr_s         = 0.1;      % padding following wavelet transform to avoid edge artifacts in seconds (bi-lateral)
cfg.(BOSC_type).pad.detection_s   = 0.1;      % padding following rhythm detection in seconds (bi-lateral); 'shoulder' for BOSC eBOSC.detected matrix to account for duration threshold
cfg.(BOSC_type).pad.background_s  = 0.1;      % padding of segments for BG (only avoiding edge artifacts)

% fooof parameters - fit with fixed line or allow a knee
cfg.(BOSC_type).fooof.aperiodic_mode    = 'fixed';
cfg.(BOSC_type).fooof.version           = 'python';
% cfg.(BOSC_type).fooof.max_n_peaks       = 1;
cfg.(BOSC_type).fooof.verbose           = true;

% threshold settings
cfg.(BOSC_type).threshold.excludePeak = [];                           % lower and upper bound of frequencies to be excluded during background fit (Hz) (previously: LowFreqExcludeBG HighFreqExcludeBG)
cfg.(BOSC_type).threshold.duration	= repmat(3, 1, numel(cfg.(BOSC_type).F)); % vector of duration thresholds at each frequency (previously: ncyc)
cfg.(BOSC_type).threshold.percentile  = .95;                              % percentile of background fit for power threshold

% episode post-processing
cfg.(BOSC_type).postproc.use      = 'no';        % Post-processing turned off for now

% general processing settings
cfg.(BOSC_type).channel           = [1]; % select posterior channels (default: all)
cfg.(BOSC_type).trial             = []; % select trials (default: all)
cfg.(BOSC_type).trial_background  = []; % select trials for background (default: all)

% Run fBOSC
clear BOSC
[BOSC, cfg] = eval([BOSC_type,'_wrapper(cfg, data)']);

% Plot the Results of the 1/f fit
% cfg.log_freqs = 0;
% cfg.plot_old = 0;
% fBOSC_fooof_plot(cfg,BOSC)

%% MY PLOT
off_before_cue = find(times{1}>=0 ,1); % samples
figure(11);
episodes = BOSC.episodes;
episodes = episodes(episodes.FrequencyMean >= freq-search_precision & episodes.FrequencyMean <= freq+search_precision,:);

for indTrial = 1:length(trials)
    origData = data.trial{indTrial}(cfg.(BOSC_type).channel(1),...
        cfg.(BOSC_type).pad.total_sample+1:end-cfg.(BOSC_type).pad.total_sample);

    origData_time = data.time{indTrial}(...
        cfg.(BOSC_type).pad.total_sample+1:end-cfg.(BOSC_type).pad.total_sample);
    
    % episoded selection
    ep = episodes(episodes.Trial==indTrial,:);

    % Plot
    hold on;
    plot(origData_time,squeeze(origData)+indTrial*plotFactor, 'k');
    ylabel({'Power';'(a.u.)'});
    set(get(gca,'YLabel'),'Rotation',45)
    
    for k=1:size(ep,1)
        onset = fix(ep.Onset(k)*fs) + off_before_cue; % offset prima dello 0
        durationS = fix(ep.DurationS(k)*fs);
        try
            eSel = origData(onset:onset+durationS);
            tSel = origData_time(onset:onset+durationS);
        catch
            eSel = origData(onset:end);
            tSel = origData_time(onset:end);
        end
        plot(tSel,eSel+indTrial*plotFactor,'r','LineWidth',2);
    end

    if indTrial == length(data.trial)
        xlabel('Time (s)');
    end
    set(findall(gcf,'-property','FontSize'),'FontSize',12)

end

xline(0,'r','LineWidth',2)
title(['Identified Episodes around ',num2str(freq),' Hz'])
print([BOSC_type, ' episodes'],'-dpng','-r300');

%% Checks

for i=1:length(data.trial)
    [pxx1(:,i),fxx]=pwelch(data.trial{i}(1:off_before_cue),[],[],[],fs);
end
m1 = mean(pxx1,2);
% figure
% semilogy(fxx,m1)

for i=1:length(data.trial)
    [pxx2(:,i),fxx]=pwelch(data.trial{i}(off_before_cue:end),[],[],[],fs);
end
m2 = mean(pxx2,2);
% figure
% semilogy(fxx,m2)

figure
nexttile;
semilogy(fxx,m1); hold on; semilogy(fxx,m2);legend('Before CUE', 'After CUE'); 
title('Averaged PSD over trials')
for i=1:5
    nexttile;
    plot(data.time{i},data.trial{i}); hold on;
    xline(0,'r')
    title(['Original Unprocessed Trial ',num2str(i)])
    xlabel('seconds')
end
