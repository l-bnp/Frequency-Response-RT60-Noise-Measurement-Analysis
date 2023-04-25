

%% 0.1 DEFINICAO DOS DIRETORIOS E ARQUIVOS %%
% Diretorio e arquivos com as respostas gravadas do local
% Usar 3 gravações para se ter uma média?
% Separar o arquivo de TSP original das gravacoes??
dir_name = '../Recorded_Signals/';
file_names = ['tsp44-rec1.wav';
              'tsp44-rec2.wav';
              'tsp44-rec3.wav'];
     
%% 0.2 PARAMETROS CONSTANTES %%
% Frequencia de resample desejada em hz
fs = 44100;
% Periodo total de avaliacao da resposta em s
T = 1.5;
% Comprimento do sinal a ser avaliado em numero de samples
L = fs * T;
% Numero de samples da resposta usados para regressão linear (0.5s)
n0 = round(fs/2);


%% 1.1 Definicao dos polos, zeros e ganho %%

% Curva A
A_zeros = [0 0 0 0];
A_poles = [-129.4 -129.4 -676.7 -4636 -76655 -76655];
A_gain = 7.39705*10^9;

% Curva C
C_zeros = [0 0];
C_poles = [-129.4 -129.4 -76655 -76655];
C_gain = 5.91797*10^9;

%% 1.2 Criacao das funcoes de transferencia, a partir dos polos e zeros %%

A_sys = zpk(A_zeros,A_poles,A_gain);
C_sys = zpk(C_zeros,C_poles,C_gain);

[Ha,wouta] = freqresp(A_sys);
fouta = wouta/(2*pi);
maga = 20*log10(squeeze(abs(Ha)));

[Hc,woutc] = freqresp(C_sys);
foutc = woutc/(2*pi);
magc = 20*log10(squeeze(abs(Hc)));

figure(1);
semilogx(fouta, maga,'LineWidth',2.0);
set(gca,'FontSize',15);
title('Função de Transferência - Curva A');
xlabel('Frequência (Hz)');
ylabel('Magnitude (dB)');
axis([10 20000 -80 20]);
grid on;

figure(2);
semilogx(foutc, magc,'LineWidth',2.0);
set(gca,'FontSize',15);
title('Função de Transferência - Curva C');
xlabel('Frequência (Hz)');
ylabel('Magnitude (dB)');
axis([10 20000 -80 20]);
grid on;

%% 1.3 Geracao do diagrama de bode %%

% % Curva A
% figure(1);
% A_h = bodeplot(A_sys);
% A_opt = getoptions(A_h);
% A_opt.FreqUnits = 'Hz';
% A_opt.Grid = 'on';
% A_opt.XLim = {[10 20000]};
% A_opt.YLim{1} = [-80 0];
% A_opt.YLim{2} = [-360 360];
% A_opt.XLabel.FontSize = 20;
% A_opt.YLabel.FontSize = 20;
% A_opt.Title.FontSize = 20;
% A_opt.TickLabel.FontSize = 15;
% A_opt.Title.String = 'Função de Transferência - Curva A';
% setoptions(A_h,A_opt);
% 
% % Curva C
% figure(2);
% C_h = bodeplot(C_sys);
% C_opt = getoptions(C_h);
% C_opt.FreqUnits = 'Hz';
% C_opt.Grid = 'on';
% C_opt.XLim = {[10 20000]};
% C_opt.YLim{1} = [-80 0];
% C_opt.YLim{2} = [-360 360];
% C_opt.XLabel.FontSize = 20;
% C_opt.YLabel.FontSize = 20;
% C_opt.Title.FontSize = 20;
% C_opt.TickLabel.FontSize = 15;
% C_opt.Title.String = 'Função de Transferência - Curva C';
% setoptions(C_h,C_opt);


%% 2.1 















