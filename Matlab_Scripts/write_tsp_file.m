% Returned values
%  x  -> TSP generated by gen_tsp2
%  t  -> time vector associated with x

T = 4800;       % duration of the TSP. The default is 1200ms.
fs = 44.1;      % sampling frequency (in kHz). The default is 44.1kHz.
bw = 22.05;     % bandwidth (in kHz) of the TSP. The default is the Nyquist frequency,
                % i.e. half of the sampling frequency.
bs = 0;         % bandshift (in kHz) of the TSP. The default is 0kHz.
ta = 240;       % starting time of the TSP (in ms). The default is T/10 ms.
tb = 200;        % group delay growing rate (in ms/kHz). The default is 50ms/kHz.
                % Note that, if a TSP with monotonically growing (or decreasing)
                % frequency, is desired, the condition T >= ta + tb*(bs+bw) must 
                % be satisfied.

[x,t] = gen_tsp3(T,fs,bw,bs,ta,tb);

figure(2);
spectrogram(x,256,250,256,44100,'yaxis');
set(gca,'FontSize',24);
%audiowrite('mytsp.wav',x,44100);
%sound(x,44100);
