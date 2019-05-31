
function [logProb, beta] = logbwd(x, means, vars, trans)

% LOGBWD Log version of the backward procedure
%
%    LOGPROB = LOGFWD(X,MEANS,VARS,TRANSITIONS) returns the likelihood of
%    the 2-dimensional sequence X (one observation per row) with respect to
%    a Markov model with N states having means MEANS and variances VARS
%    (stored in N elements lists with empty matrices as first and last
%    elements to symbolize the entry and exit states) and transition matrix
%    TRANSITIONS.
%      Alternately, LOGFWD(X,HMM) can be used, where HMM is an object of the
%    form:
%       HMM.means = MEANS;
%       HMM.vars = VARS;
%       HMM.trans = TRANSITIONS
%

if nargin == 2,
  model = means;
  means = model.means;
  vars = model.vars;
  model.trans(model.trans<1e-100) = 1e-100;
  logTrans = log(model.trans);
end;

nEst = length(means);	% Cantidad de estados
nEstInic = 2;
nEstFinal = nEst-1;	% Indice donde está el último estado al que puede ir SIN TERMINAR
[cant_pts, dim] = size(x);

% Inicializo beta
log2pi = log(2*pi);
for i = nEstInic:nEstFinal
	constante = -1/2 * log(det(vars{i})) - log2pi;
	invSig{i} = inv(vars{i});
	X = x(1,:) - means{i}';
	log_bj = constante - 1/2* (X*invSig{i})*X';
	beta(i) = log_bj + logTrans(1,i);
end


% Hago la recursión backward
for t = fliplr(1:cant_pts)	% Itero en las observaciones
	beta_prev = beta;
	for i = nEstInic:nEstFinal % Itero en los estados
		constante = -1/2 * log(det(vars{i})) - log2pi;
		invSig{i} = inv(vars{i});
		X = x(t+1,:) - means{i}';
		log_bj = constante - 1/2* (X*invSig{i})*X';
		beta(i) = logsum(beta_prev(nEstInic:nEstFinal) + log_bj + logTrans(2:nEstFinal,i));
	end
end


