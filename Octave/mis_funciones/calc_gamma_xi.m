
function [gama, xi] = calc_gamma_xi(x, hmm)

	means = hmm.means(2:end-1);
	vars = hmm.vars(2:end-1);
	hmm.trans(hmm.trans<1e-100) = 1e-100;
	logTrans = log(hmm.trans);
	logAlfa = hmm.logAlfa;
	logBeta = hmm.logBeta;

	[T, K] = size(logAlfa);	% T = tiempo_max; K = cant_clases

	gama = zeros(T,K);
	xi = zeros(T,K,K);

	% Calculo gama
	for t = (1:(T))
		den = logsum(logAlfa(t,:) + logBeta(t,:));
		for j = 1:K
			gama(t, j) = (logAlfa(t,j)+logBeta(t,j))-den;
		end
	end


	% Calculo xi
	log2pi = log(2*pi);
	for t = (2:(T))
		den = logsum(logAlfa(t,:) + logBeta(t,:));
		for i = 1:K
			for j = 1:K
				constante = -1/2 * log(det(vars{j})) - log2pi;
				invSig{j} = inv(vars{j});
				X = x(t-1,:) - means{j}';
				log_bj = constante - 1/2* (X*invSig{j})*X';

				xi(t,i,j) = (logAlfa(t-1,i)+logBeta(t,j)+log_bj+logTrans(i,j)) - den;
			end
		end
	end

	% Verifico que con xi obtengo gama "integrando xi con respecto a j"
	t = 3; i = 3;
	gama_t_i = [logsum(xi(i,:,t)) gama(i,t)]


	gama = exp(gama);
	xi = exp(xi);
end
