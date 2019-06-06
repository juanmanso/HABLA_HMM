
function [gama, xi] = calc_gamma_xi(x, hmm)

	means = hmm.means(2:end-1);
	vars = hmm.vars(2:end-1);
	hmm.trans(hmm.trans<1e-100) = 1e-100;
	logTrans = log(hmm.trans(2:end-1 ,2:end-1));
	logAlfa = hmm.logAlfa;
	logBeta = hmm.logBeta;

	[T, K] = size(logAlfa);	% T = tiempo_max; K = cant_clases

	gama = zeros(T,K);
	xi = zeros(K,K,T);

	% Calculo gama
	for t = (1:(T))
		den = logsum(logAlfa(t,:) + logBeta(t,:));
		for j = 1:K
			gama(t, j) = (logAlfa(t,j)+logBeta(t,j))-den;
		end
	end


	% Calculo xi
	log2pi = log(2*pi);
	for t = 2:T
		den = logsum(logAlfa(t,:) + logBeta(t,:));
		for j = 1:K
			for k = 1:K
				constante = -1/2 * log(det(vars{k})) - log2pi;
				invSig{k} = inv(vars{k});
				X = x(t,:) - means{k}';
				log_bj = constante - 1/2* (X*invSig{k})*X';

				xi(j,k,t) = (logAlfa(t-1,j)+logBeta(t,k)+log_bj+logTrans(j,k)) - den;
			end
		end
	end

	% Verifico que con xi obtengo gama "integrando xi con respecto a k"
%	t = 3; j = 3;
%	gama_t_j = [logsum(xi(j,:,t)) gama(t-1,j)]

	for t = 2:T
		for j = 1:K
			verif(t,k) = logsum(xi(j,:,t)) - gama(t-1,j);
		end
	end

	if(verif < 1e-10)
		puts("Ok Gama y Psi \n");
	else
		puts("Gama y Psi MAL \n");
	end


	gama = exp(gama);
	xi = exp(xi);
end
