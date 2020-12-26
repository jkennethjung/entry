clc;
clear;
close all;
echo off;

diary ../output/analysis.log
diary on;
rng(1);

%%% I. Globals %%%

global E mu sigma M alpha beta mu_eps beta_eps mean_eps Q gamma eta Nq_max;

E = 50;
mu = 1;
sigma = 1;
M = 10;
alpha = 1;
beta = 1;
gamma = 1;
eta = 1;
Nq_max = 4;

mu_eps = 0;
beta_eps = 1;
mean_eps = mu_eps + beta_eps*0.57721;
Q = 6;
Qs = zeros(Q, 1);
Qv = zeros(Q, 1);
for q = 0:(Q-1)
    Qs(q + 1) = 1/(2*Q) + q/Q;
    Qv(q + 1) = icdf('Lognormal', Qs(q + 1), mu, sigma);
end

%%% II. Draw Characteristics %%%

A = random('Lognormal', mu, sigma, E, 1);
Aq = zeros(E, 1);
Aq_idx = zeros(E, 1);
for j = 1:E
    q = 1;
    qj = 0;
    while qj == 0
        if q < Q
            if (A(j) - Qv(q) > 0) & (A(j) - Qv(q + 1) <= 0)
                d1 = A(j) - Qv(q);
                d2 = Qv(q) - A(j);
                if d1 == min(d1, d2)
                    qj = q;   
                else
                    qj = q + 1;
                end
            end
        else
            qj = q;
        end
        q = q + 1;
    end 
    Aq_idx(j) = qj;
    Aq(j) = Qv(qj);
end
A_hist = zeros(Q,1);
for q = 1:Q
    flags = (Aq_idx == q);
    A_hist(q) = sum(flags);
end

epsilon = -evrnd(mu_eps, beta_eps, [M,E]);

W = random('Lognormal', mu, sigma, M, 1);
R = random('Lognormal', mu, sigma, M, 1);
X = random('Lognormal', mu, sigma, M, 1);

%%% III. Construct State Space %%%

global S States A_hist PiV;

S = (Nq_max)^Q;
States = zeros(S, Q);
State = zeros(1, Q);
for s = 2:S
    maxed = (State == Nq_max);
    if maxed(Q) == 0
        State(Q) = State(Q) + 1;
    else
        q = Q;
        q0 = 0;
        while q0 == 0 & q > 0
            if maxed(q) == 1
                q = q - 1;
            else 
                q0 = q;
            end
        end
        State(q0) = State(q0) + 1;
        for q = (q0+1):Q
            State(q) = 0;
        end
    end
    States(s, :) = State;
end

%%% IV. Calculate Cournot Payoffs in Each State %%%

P = zeros(M, S);
Y = zeros(M, S, Q);
PiV = zeros(M, S, Q);

for m = 1:M
    w = W(m);
    r = R(m);        
    x = X(m);
    S_neg = zeros(S, 1);
    MC = zeros(Q, 1);
    for q = 1:Q
        a = Qv(q);
        MC(q) = marginal_cost(a, w, r, alpha);
    end
    for s = 2:S
        State = States(s, :);
        C = zeros(Q, 1);
        Qs = [0];    % Number of firms by quantile (excludes zeros)
        Cs = [0];    % Constant by quantile (excludes zeros)
        for q = 1:Q
            Qn = State(q);
            if Qn > 0
                a = Qv(q);
                c = (x*gamma - MC(q))/eta;
                if c > 0
                    C(q) = c;
                    Qs = [Qs; Qn];
                    Cs = [Cs; c];
                end
            end
        end

        Iq = (C > 0); % Dummies for quantile inclusion
        Nq = sum(Iq); % Number of quantiles represented in eqm
        if Nq > 0
            Qs = Qs(2:(Nq+1)); 
            Cs = Cs(2:(Nq+1)); 
            B = zeros(Nq);
            for r = 1:Nq
                row = Qs;
                row(r) = row(r) + 1;
                B(r, :) = row;
            end
            y = inv(B)*Cs;
            if sum(y < 0) == 0
                S_neg(s) = 1;
                p = x*gamma - y' * Qs;
                P(m, s) = p;
                j = 1;
                for q = 1:Q
                    if Iq(q) == 1
                        Y(m, s, q) = y(j);
                        PiV(m, s, q) = (p - MC(q)) * y(j);
                        j = j + 1;
                    end
                end
            end
        end 
    end
    disp('% of states with negative output');
    mean(S_neg)
end

cp_mat = ones(M, Q)/(M*Q); 
cp = cp_mat(:);
lb = zeros(M, Q);
ub = ones(M, Q);
sparse = zeros(M*Q, M*Q);
k = 1;
for q = 1:Q
    for m = 1:M
        for t = 1:Q
            for l = 1:M
                if (m == l) | (q == t)
                    sparse(k) = 1;
                end
                k = k + 1;
            end
        end
    end
end
options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective', ...
    'JacobPattern', sparse);
[cp, resnorm, residual, exitflag, output] = lsqnonlin(@fpe, cp, ...
    lb, ub, options);


%%% V. Construct Firm Beliefs %%%

function z = fpe(cp)
    global M Q States A_hist PiV;
    cp_mat = reshape(cp, [M, Q]);
    z = choice_prob(exp_profit(cp_mat, States, A_hist, PiV)) - cp_mat;
    disp(z);
end

function cp = choice_prob(EPi)
    global Q M;
    cp = zeros(M, Q);
    exp_EPi = exp(EPi);
    for q = 1:Q
        denom = sum(exp_EPi(:, q));
        cp(:, q) = exp_EPi(:, q)/denom;
    end
end

function EPi = exp_profit(p, States, A_hist, PiV) 
    global M Q S;
    EPi = zeros(M, Q);
    for q = 1:Q
        for m = 1:M
            p_s = zeros(S, 1);
            for s = 1:S
                p_s(s) = pr_state(p(m, :), s, States, A_hist);
            end
            EPi(m, q) = PiV(m, :, q)*p_s;
        end
    end
end

function p_s = pr_state(pr, s, States, A_hist) 
    global Q;
    p_s = 1;
    impossible = 0;
    for q = 1:Q
        n = A_hist(q);
        k = States(s, q);
        if (n < k) | impossible
            p_s = 0;
            impossible = 1;
        else 
            pr_q = pr(q);
            p_s = p_s*nchoosek(n, k)*pr_q^k*(1-pr_q)^(n-k); 
        end
    end
end

function mc = marginal_cost(a, w, r, alpha)
     mc = [w*(alpha*r/((1-alpha)*w))^(1-alpha) + ...
           r*((1-alpha)*w/(alpha*r))^alpha]/a;
end 
