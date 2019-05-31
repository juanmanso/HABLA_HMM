%% inic_hmm.m

close all
clear all

addpath('./mis_funciones')
load data;
%colordef none;
%whos

% hmm tiene: hmm.vars, hmm.means y hmm.trans

% Ahora calculamos p(x,s) para cada modelo y vemos cuál es el más verosimil.

hmms = [hmm1; hmm2; hmm3; hmm4; hmm5; hmm6];
p_hmm = zeros(1,length(hmms));

for i = 1:length(X2)

	qt = ST2(i+1);
	qprev = ST2(i);

	for k = 1:length(hmms)
		hmm = hmms(k);
		p_hmm(k) += log(mvnpdf(X2(i,:), hmm.means{qt}, hmm.vars{qt})*hmm.trans(qprev,qt));
	end
end
% Le sumo (porque estoy log) la prob de transición al último estado.
for k = 1:length(hmms)
	p_hmm(k) += log(hmms(k).trans(end-1,end));
end

figure
hold on
stem(p_hmm)
[maximo, estado_max] = max(p_hmm);
stem(estado_max, maximo,'r')
xlabel('Estados')
ylabel('Log(Likelihood)')
legend('Todas las p(x,\lambda', 'p(x,\lambda) máxima')
legend('location','SouthWest')

figure
plotseq2(X2, ST2, hmm3);

% Otra forma de hacerlo es
% log(bj(x)) = -log(2pi) -1/2*det(covj) - 1/2(x-muj)covj^-1(x-muj)'

%% Forward
% Para calcular el P(x) que escribiremos como P(x_1^t, x_t+1^T)
% Con la recursión forward se escribe P(x_1^t, s_t = j) y la sumatoria de todos esos se llega a P(x)
% En habla, como está el estado 5, hay que decirles que termina en 9.

% P(x, S[qué es max]) = Viterbi.

% P(x, Sm) ~ P(X) si se supone que la máxima es mucho mayor que las demás.
% En la realidad no se suele cumplir.

% P(x2) [en logfwd] >~ P(X2,SM) [en Viterbi]

%% Vamos a calcular P(X) con forward y se compara con P(X,Sm) a ver si da igual o no

p_fwd = zeros(1,length(hmms));
alfa = zeros(length(X2), max(ST2)-1,length(hmms));

for k = 1:length(hmms)
	[p_fwd(k), alfa(:,:,k)] = logfwd(X2, hmms(k));
end

alfa_hmm1 = alfa(:,:,1);
alfa_hmm2 = alfa(:,:,2);
alfa_hmm3 = alfa(:,:,3);
alfa_hmm4 = alfa(:,:,4);
alfa_hmm5 = alfa(:,:,5);
alfa_hmm6 = alfa(:,:,6);

% Le falta sumarle el último
p_fwd2(1) = logsum(alfa(end,2:end,1));
p_fwd2(2) = logsum(alfa(end,2:end,2));
p_fwd2(3) = logsum(alfa(end,2:end,3));
p_fwd2(4) = logsum(alfa(end,2:end,4));
p_fwd2(5) = logsum(alfa(end,2:end,5));
p_fwd2(6) = logsum(alfa(end,2:end,6));


%% Viterbi
%% Viterbi se hace idem pero calculando P(X,S) con viterbi calculando con P(X,SM)

secuencia_Vit = zeros(length(ST2),length(hmms));
logpx = zeros(1, length(hmms));

for k = 1:length(hmms)
	[secuencia_Vit(:,k), logpx(k)] = logvit(X2, hmms(k));
end

%% Comparación de ambos
figure
hold on
stem(p_hmm)
stem(p_fwd)
stem(logpx)
legend('Probabilidad hecha a mano', 'Probabilidad con Forward', 'Probabilidad con Viterbi')
legend('location', 'SouthWest')

[p_hmm; p_fwd; logpx]

return

%% Entrenamiento
