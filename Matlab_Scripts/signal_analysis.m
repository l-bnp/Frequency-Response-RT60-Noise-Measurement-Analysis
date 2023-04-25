%% 0.1 DEFINICAO DOS DIRETORIOS E ARQUIVOS %%
% Diretorio e arquivos com as respostas gravadas do local
% Usar 3 gravações para se ter uma média?
% Separar o arquivo de TSP original das gravacoes??
dir_name = '../Recorded_Signals/';
file_names = ['TSP_Resp_1_Aud.wav';
              'TSP_Resp_2_Aud.wav';
              'TSP_Resp_3_Aud.wav'];
     
%% 0.2 PARAMETROS CONSTANTES %%
% Frequencia de resample desejada em hz
fs = 44100;
% Periodo total de avaliacao da resposta em s
T = 3;
% Comprimento do sinal a ser avaliado em numero de samples
L = fs * T;
% Numero de samples da resposta usados para regressão linear (0.5s)
n0 = round(fs/2);
          
%% 0.3 FREQUENCIAS DE AVALIACAO DA RESPOSTA %%
% Faixas de frequência a terem sua resposta individual analisada
% Separação equivalente a terços de oitava
%fc = 100*nthroot(2,3).^(0:23);
%fc = round(fc);
%fc(24) = 19599;
fc = [100 125 160 200 250 315 400 500 630 800 1000 1250 1600 2000 2500 3150 4000 5000 6300 8000 10000 12500 16000 19599];

%% 0.4 GERAR FILTROS E GRAFICO DE FUNCAO DE TRANSFERENCIA DOS FILTROS %%
% figure(1);
for n = 1:length(fc)
    % Cria os coeficientes de um filtro de butterworth em torno da fc
    [b(n,:),a(n,:)] = butter(3,fc(n)*[8/9 9/8]/(fs/2));
    % Salva a resposta em frequência do filtro numa matriz
    [h(n,:), w(n,:)] = freqz(b(n,:),a(n,:), fc, fs);
    % Plota a resposta em frequência em dB
    figure(1);
    subplot(5,5,n);
    semilogx(fc, 20*log10(abs(h(n,:))),'LineWidth',2.0);
    %set(gca,'FontSize',12);
    axis([100 20000 -60 0]);
    xticks([100 1000 10000]);
    yticks([-60 -40 -20 0]);
    title(['fc = ' num2str(fc(n)) ' Hz']);
    grid on;
end
  
%% LOOP PARA TODAS AS GRAVACOES %%
for k = 1:3
    
    %% 1.1 LEITURA DE ARQUIVOS E SAMPLING %%
    % Leitura do arquivo
    [x1,fs1] = audioread([dir_name,file_names(k,:)]);
    % Resample na frequencia desejada
    x = resample(x1,fs,fs1);    
    %    [b,a] = butter(3,fc(1)/2/(fs/2),'high');
    %    x = filtfilt(b,a,x);
    
    %% 1.2 CONVOLUCAO E MAGNITUDE DA RESPOSTA TEMPORAL %%
    % Convolui g(t) com f(-t). Resultado é a resposta ao impulso 
    % no dominio do tempo
    y_lin = conv(x(:,1),flipud(x(:,2))); 
    % Magnitude em escala de dB do sinal convoluido
    y = 20*log10(abs(y_lin));
    
    %% 1.3 CÁLCULO DA RESPOSTA EM FREQUÊNCIA %%
    % Compute the Fourier transform of the signal
    Y_dft = fft(y_lin);
    % Compute the two-sided spectrum P2. Then compute the single-sided 
    % spectrum P1 based on P2 and the even-valued signal length L
    P2 = 20*log10(abs(Y_dft/L));
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    % Define the frequency domain f and plot the single-sided 
    % amplitude spectrum P1
    f = fs*(0:(L/2))/L;
    figure(2)
    subplot(3, 1, k);
    semilogx(f,P1) 
    set(gca,'FontSize',15);
    title(['Resposta em Frequência - Medição ' num2str(k)])
    xlabel('Frequência (Hz)')
    ylabel('Amplitude (dB)');
    xlim([100 20000]);
    ylim([-80 40]);
    yticks([-80 -60 -40 -20 0 20 40])
    grid on;
    
    %% 1.4 VALOR MAXIMO E SELECAO DE REGIAO DA AMOSTRA %%
    % Valor maximo de y e indice correspondente a posicao desse valor 
    [y_max,idx_max] = max(y);    
    % Coloca a amostra de valor max no índice 1 e seleciona as proximas
    % T*fs amostras, que correspondem a T segundos
    y = y(idx_max+(0:T*fs));
    
    %% 1.5 VETOR TEMPORAL E REGRESSAO LINEAR %%
    % Vetor de tempo
    t = (0:(length(y)-1))/fs;
    % Coeficientes do polinomio de regressao linear da  
    % resposta temporal y do ambiente no intervalo de 1 a n0
    coef = polyfit(t(1:n0)',y(1:n0),1);
    % Polinomio correspondente
    p = polyval(coef,t);
    
    %% 1.6 GRAFICO DE RESPOSTA GERAL DO AMBIENTE %%
    % Figura: Resposta temporal do ambiente considerando todo o espectro
    figure(3);
    subplot(3, 1, k);
    % Plot de y e p normalizados em relacao a amplitude de p(1)
    plot(t,y-p(1),'b',t, p-p(1),'r','LineWidth',2.0) ;   
    % Parametros do grafico
    set(gca,'FontSize',15);
    axis([0 T -80 20]);
    title(['Resposta ao Impulso - Medição ' num2str(k)]);
    xlabel('Tempo (s)');
    ylabel('Amplitude (dB)');
    legend({'Resposta Temporal Medida','Polinômio Aproximado'},'Location','northeast');
    grid on;
    
    %% 1.7 CALCULO DO RT60 GERAL %%
    % Tempo de reverberação geral
    % Calcula em qual amostra são atingidos os 20dB de decaimento
    n20 = length(t);
    for m = 1:length(t)
        if (p(1)-p(m)) >= 20
            n20=m;
            break;
        end
    end
    % T20
    T20(k) = t(n20)-t(1);
    % T60: Multiplica-se o T20 por 3
    T60(k) = T20(k)*3;
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %% LOOP PARA TODAS AS FREQUENCIAS DE AVALIACAO %%
    for n = 1:length(fc)
        
        %% 2.1 FILTRAGEM DA RESPOSTA COMPLETA %%
        % Filtra o sinal da resposta do ambiente em fc
        y_lin_fc = filtfilt(b(n,:),a(n,:),y_lin);
        % Passa pra escala de dB
        y_fc = 20*log10(abs(y_lin_fc));
        
        %% 2.2 VALOR MAXIMO, SELECAO DE REGIAO DA AMOSTRA %%
        % Valor maximo de y e índice correspondente à posição desse valor 
        [y_fc_max,idx_fc_max] = max(y_fc);    
        % Coloca a amostra de valor max no índice 1 e seleciona as proximas
        % T*fs amostras, que correspondem a T segundos
        y_fc = y_fc(idx_fc_max+(0:T*fs));        

        %% 2.3 VETOR TEMPORAL E POLINOMIO DE APROXIMACAO %%
        % Vetor de tempo
        t = (0:(length(y_fc)-1))/fs;        
        % Coeficientes do polinomio que descreve (Aproximadamente) a  
        % resposta temporal y do ambiente no intervalo de 1 a n0
        coef_fc = polyfit(t(1:n0)',y_fc(1:n0),1);
        % Polinomio correspondente
        p_fc = polyval(coef_fc,t);
        
        %% 2.4 PLOT DA RESPOSTA ESPECIFICA DA FREQUENCIA %%
        % Figura: Resposta temporal do ambiente em cada frequência fc
        figure(3+k);
        % Posiciona no subplot de 4x4
        subplot(8,3,n);
        % Plota  a resposta temporal para a frquencia fc em azul e a
        % aproximacao polinomial em vermelho
        plot(t,y_fc-p_fc(1),'b',t,p_fc-p_fc(1),'r','LineWidth',2.0);       
        % Detalhes do grafico
        axis([0 T -80 20]);
        title(['fc = ' num2str(fc(n)) ' Hz']);
        grid on;
        
        %% 2.5 CALCULO DO RT60 ESPECIFICO DA FREQUENCIA %%
        % Calcula em qual amostra são atingidos os 20dB de decaimento
        n20 = length(t);
        for m = 1:length(t)
            if (p_fc(1)-p_fc(m)) >= 20
                n20=m;
                break;
            end
        end
        % T20
        T20_fc(n, k) = t(n20)-t(1);
        % T60: Multiplica-se o T20 por 3
        T60_fc(n, k) = T20_fc(n, k)*3;
        
    end
end

%% 3.1 GRAFICO DE T60 POR FREQUENCIA %%
T60_fc_transp = transpose(T60_fc);
T60_fc_mean = mean(T60_fc_transp);
T60_fc_std = std(T60_fc_transp);
figure(7);
semilogx(fc, T60_fc_mean,'b', fc, T60_fc_std, 'r','LineWidth',2.0);
% Parametros do grafico
set(gca,'FontSize',15);
axis([100 20000 0 max(T60_fc_mean)]);
title('RT60 Médio em Função da Frequência');
xlabel('Frequência (Hz)');
ylabel('Segundos');
legend({'RT60 Médio','Desvio Padrão'},'Location','northeast');
grid on;
    
    