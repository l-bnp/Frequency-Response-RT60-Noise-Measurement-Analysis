% Source file to evaluate noise level
source = dsp.AudioFileReader('Noise_Aud.wav');
fs = source.SampleRate;

player = audioDeviceWriter('SampleRate',fs);

% Osciloscope vision
scope  = timescope('SampleRate',fs, ...
    'TimeSpanOverrunAction','Scroll', ...
    'TimeSpanSource','Property','TimeSpan',20,'ShowGrid',true, ...
    'YLimits',[20 80],'AxesScaling','Auto', ...
    'ShowLegend',true,'BufferLength',4*3*fs, ...
    'ChannelNames', ...
    {'Lt_CF','Leq_C','Lpeak_C','Lmax_CF'}, ...
    'Name','Sound Pressure Level C-Weighted');

% SPL meter object to make measures
SPL = splMeter('TimeWeighting','Fast', ...
    'FrequencyWeighting','C-weighting', ...
    'SampleRate',fs, ...
    'TimeInterval',2);

% 1 kHz Reference tone to calibrate
[test_sig,test_fs] = audioread(['test_tone_1khz.wav']);
calibrate(SPL, test_sig, 91);

% Plot values
while ~isDone(source)
    x = source();
    player(x);
    [Lt,Leq,Lpeak,Lmax] = SPL(x);
    scope([Lt,Leq,Lpeak,Lmax])
end

release(source)
release(player)
release(SPL)
release(scope)

