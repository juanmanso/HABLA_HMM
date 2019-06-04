
function [logProb, logBeta] = logbwd(x, means, vars, trans, logAlfa)

% LOGBWD Log version of the backward procedure
%
%    LOGPROB = LOGBWD(X,MEANS,VARS,TRANSITIONS,ALFA) returns the likelihood of
%    the 2-dimensional sequence X (one observation per row) with respect to
%    a Markov model with N states having means MEANS and variances VARS
%    (stored in N elements lists with empty matrices as first and last
%    elements to symbolize the entry and exit states) and transition matrix
%    TRANSITIONS.
%      Alternately, LOGFWD(X,HMM) can be used, where HMM is an object of the
%    form:
%       HMM.means = MEANS;
%       HMM.vars = VARS;
%       HMM.trans = TRANSITIONS;
%				HMM.alfa = ALFA;
%

if nargin == 2,
  model = means;
  means = model.means;
  vars = model.vars;
	logAlfa = model.logAlfa;
  model.trans(model.trans<1e-100) = 1e-100;
  logTrans = log(model.trans);
end;

nEst = length(means);	% Cantidad de estados
nEstInic = 2;
nEstFinal = nEst-1;	% Indice donde está el último estado al que puede ir SIN TERMINAR
[cant_pts, dim] = size(x);

% Inicializo beta
log2pi = log(2*pi);

logBeta = logTrans(1:end-1,end)';

%for i = nEstInic:nEstFinal
%	constante = -1/2 * log(det(vars{i})) - log2pi;
%	invSig{i} = inv(vars{i});
%	X = x(1,:) - means{i}';
%	log_bj = constante - 1/2* (X*invSig{i})*X';
%	beta(i) = log_bj + logTrans(end,i);	% VA CON LOGBJ o no?!
%end
%
%logBeta = beta(:)';

beta = logBeta';
% Hago la recursión backward
for t = fliplr(1:cant_pts-1)	% Itero en las observaciones
	beta_next = beta;

	for i = nEstInic:nEstFinal
		constante = -1/2 * log(det(vars{i})) - log2pi;
		invSig{i} = inv(vars{i});
		X = x(t+1,:) - means{i}';
		log_bj(i) = constante - 1/2* (X*invSig{i})*X';
	end

	for i = nEstInic:nEstFinal % Itero en los estados
		beta(i) = logsum(beta_next(nEstInic:nEstFinal)' + log_bj(2:end) + logTrans(i,nEstInic:nEstFinal));
	end

	logBeta = [beta(:)'; logBeta];	% Almaceno
end


% Calculo logProb

%	logProb = 0;
%	for k = nEstInic:nEstFinal
%		logProb += logsum(logBeta(:,k) + logAlfa(:,k-1));
%	end

logProb = logsum(log_bj(2:end) + logBeta(1,2:end) + logTrans(1,2:4));

end
