%% Algoritmo de Baum-Wech (EM + HMM)

inic_hmm

% hmm tiene: hmm.vars, hmm.means y hmm.trans

%% Qué modelo corresponde con ST2
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
legend('Todas las p(x,\theta)', 'p(x,\theta) máxima')
legend('location','SouthWest')

figure
plotseq2(X2, ST2, hmm3);

% Otra forma de hacerlo es
% log(bj(x)) = -log(2pi) -1/2*det(covj) - 1/2(x-muj)covj^-1(x-muj)'


%---------------------------------------------------------%

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

for k = 1:length(hmms)
	% Se obtienen los alfa y P(X) por medio de logfwd()
	[p_fwd(k), logAlfa] = logfwd(X2, hmms(k));
	hmms(k).logAlfa = logAlfa(:, 2:end);

	% Cálculo del P(X) utilizando los alfa a mano.
		logTrans = hmms(k).trans(2:4,end);
		logTrans(logTrans<1e-100)=1e-100;
		logTrans = log(logTrans);
	p_fwd2(k) = logsum(hmms(k).logAlfa(end,:) + logTrans');
end


%% Viterbi
%% Viterbi se hace idem pero calculando P(X,S) con viterbi calculando con P(X,SM)

secuencia_Vit = zeros(length(ST2),length(hmms));
logpx = zeros(1, length(hmms));

for k = 1:length(hmms)
	[secuencia_Vit(:,k), logpx(k)] = logvit(X2, hmms(k));
end


%% Idem P(X) pero backward

p_bwd = zeros(1,length(hmms));

for k = 1:length(hmms)
	% Se obtienen los beta y P(X) por medio de logbwd()
	[p_bwd(k), logBeta] = logbwd(X2, hmms(k));
	hmms(k).logBeta = logBeta;
end

%% Comparación de ambos
figure
hold on
stem(p_hmm)
stem(p_fwd)
stem(p_fwd2)
stem(p_bwd)
stem(logpx)
legend('Probabilidad hecha a mano', 'Probabilidad con Forward', 'Probabilidad con Forward a mano', 'Probabilidad con Backward', 'Probabilidad con Viterbi')
legend('location', 'SouthWest')
xlabel('Estados')
ylabel('Log(Likelihood)')
title('Comparación entre métodos de cálculo de P(X)')


[p_hmm; p_fwd; p_fwd2; p_bwd; logpx]

% Verifico alfa y beta


for k = 1:6
	alfa = hmms(k).logAlfa;
	beta = hmms(k).logBeta;
	p_verif(k,:) = zeros(1,length(alfa));

	for i = 1:length(alfa)
		p_verif(k,i) = logsum(alfa(i,:) + beta(i,:));
	end
end

if(p_verif&&p_verif(:,1))
	puts("Ok Alfa y Beta \n");
else
	puts("Alfa y Beta MAL \n");
end

[gama, xi] = calc_gamma_xi(X2, hmms(1));


%% Entrenamiento
