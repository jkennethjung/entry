clc;
clear;
close all;
echo off;

diary ../output/analysis.log
diary on;
rng(1);

%%% I. Initialize  %%%

parpool(36);
load_as = '../temp/data.csv';
save_as = '../output/estimates.csv';
NS = 10;
NBS = 20;
mu = 1;
sigma = 1;
alpha = 0.6;
gamma = 1;
eta = 1;
Nq_max = 9;

mu_eps = 0;
beta_eps = 1;
mean_eps = mu_eps + beta_eps*0.57721;
Q = 4;
Qv = zeros(Q, 1);
for q = 1:Q
    Qv(q + 1) = 1.5*q;
end
ncol = 9;

data = readmatrix(load_as);
T = max(data(:, 1));
M = zeros(T, 1);
E = zeros(T, 1);
W = cell(T, 1);
R = cell(T, 1);
X = cell(T, 1);
for t = 1:T
    t_rows  = (data(:, 1) == t);
    data_t = data(t_rows, :);
    M(t) = max(data_t(:, 2));
    E(t) = sum(t_rows);
    disp('Number of entrants')
    disp(E(t))
    W{t} = zeros(M(t), 1);
    R{t} = zeros(M(t), 1);
    X{t} = zeros(M(t), 1);
    for m = 1:M(t)
        m_rows = (data_t(:, 2) == m);
        data_mt = data_t(m_rows, :);
        W{t}(m) = data_mt(4);
        R{t}(m) = data_mt(5);
        X{t}(m) = data_mt(6);
    end
end

%%% II. Draw Characteristics %%%
epsilon = cell(T, NS);
for t = 1:T
    for ns = 1:NS
        epsilon{t, ns} = -evrnd(mu_eps, beta_eps, [M(t), E(t)]);
    end
end

%%% III. Construct State Space %%%

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

disp('True auxiliary model parameters');
Z = data(:, 3:6);
y = data(:, 7);
Beta_1 = [mean(y); var(y)];
Beta_2 = (Z.'*Z)^(-1)*(Z.'*y);
Beta_0 = [Beta_1; Beta_2];
disp(Beta_0);

Theta = zeros(4, NBS);
%nbs = 1;
%for nbs = 1:NBS
%    msg = 'Initiating estimation for bootstrap = %d';
%    str = sprintf(msg, nbs)
    theta = [gamma; eta; mu; sigma];
    c = clock;
    fix(c)
    aux = @(theta) auxiliary(theta, Beta_0, NBS, NS, alpha, M, E, W, R, ...
            X, Q, Qv, T, S, States, epsilon, ncol);
    [theta, dBeta, exitflag, output] = fmincon(aux, theta, [], [], [], [], ...
        theta*1e-2, theta*1e2);
    disp('Outer loop optimization finished');
    disp(theta);
    c = clock;
    fix(c)
%    Theta(:, nbs) = theta;
%end
writematrix(Theta, save_as);

%%% SIMULATE() MUST RETURN A DATASET INSTEAD OF A BETA
%%% THIS DATASET MUST BE ASSEMBLED OVER ALL BOOTSTRAPS
%%% AND THEN AUXILIARY() CAN RETURN BETA

%%% NOTE THAT NBS INSTEAD OF nbs IS NOW THE ARGUMENT OF AUXILIARY()

function dBeta = auxiliary(theta, Beta_0, NBS, NS, alpha, M, E, W, R, X, Q, ...
        Qv, T, S, States, epsilon, ncol)
    disp('theta:')
    disp(theta)
    gamma = theta(1);
    eta = theta(2);
    mu = theta(3);
    sigma = theta(4);
    A_draw = draw_firms(mu, sigma, NBS, E, Q, Qv, T);
    A_hist = A_draw{1};
    Aq_idx = A_draw{2};
    data_sim = zeros(1, ncol);
    for nbs = 1:NBS
        data_nbs = simulate(theta, nbs, A_hist, Aq_idx, NS, alpha, M, E, W, R, X, Q, Qv, T, S, ...
            States, epsilon, ncol);
        data_sim = [data_sim; data_nbs];
    end
    data_sim(1, :) = [];
    Z = data_sim(:, 3:6);
    y = data_sim(:, 7);

    Beta_1 = [mean(y); var(y)];
    Beta_2 = (Z.'*Z)^(-1)*(Z.'*y);
    Beta = [Beta_1; Beta_2];
    dBeta = (Beta_0 - Beta).'*(Beta_0 - Beta);
    disp('Objective function');
    disp(dBeta);
end

function data_sim = simulate(theta, nbs, A_hist, Aq_idx, NS, alpha, M, E, W, R, X, Q, Qv, T, S, ...
        States, epsilon, ncol)
    gamma = theta(1);
    eta = theta(2);
    mu = theta(3);
    sigma = theta(4);

    data_sim = zeros(1, ncol);
    for t = 1:T
        %%% IV. Calculate Cournot Payoffs in Each State %%%
        % This may require large amount of memory, so placing in inner loop of t
        % rather than over all t beforehand
        
        M_t = M(t);
        E_t = E(t);
        P = zeros(M_t, S);
        Y = zeros(M_t, S, Q);
        L = zeros(M_t, S, Q);
        K = zeros(M_t, S, Q);
        PiV = zeros(M_t, S, Q);
        W_t = W{t};
        R_t = R{t};
        X_t = X{t};
        
        for m = 1:M_t
            w = W_t(m);
            r = R_t(m);        
            x = X_t(m);
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
            %disp('% of states with negative output');
            %mean(S_neg)
        end
        
        for m = 1:M_t
            L(m, :, :) = (1-alpha)*Y(m, :, :)/W_t(m);
            K(m, :, :) = (alpha)*Y(m, :, :)/R_t(m);
        end
    
        %%% V. Fixed Point to Solve for Firm Beliefs%%%
    
        [cp, resnorm, residual, exitflag, output] = nfp(t, M_t, Q, States, ...
            A_hist, PiV, S, nbs);
        cp_mat = reshape(cp, [M_t, Q]);
        At_hist = A_hist{t, nbs}(:);
        EPi = exp_profit(cp_mat, States, At_hist, PiV, M_t, Q, S);
        Beta = zeros(4, NS);
        for ns = 1:NS
            %%% VI. Entry and Ex-Post Outcomes  %%%
            M_J = zeros(E_t, 1);
            S_M = zeros(M_t, Q);
            for j = 1:E_t
                q = Aq_idx{t, nbs}(j);
                shocks = reshape(epsilon{t, ns}(:, j), [M_t, 1]);
                Pi_j = EPi(:, q) + shocks ;
                [Pi_mj, m_j] = max(Pi_j);
                M_J(j) = m_j;
                S_M(m_j, q) = S_M(m_j, q) + 1;
            end
            
            State_M = zeros(M_t, 1);
            data_c = zeros(1, ncol);
            for m = 1:M_t
                s_m = S_M(m, :);
                n_m = sum(s_m);
                [~, s] = ismember(s_m, States, 'rows');
                if s == 0
                    disp('Warning: eqm state is outside state space');
                    disp(s_m);
                else
                    State_M(m) = s;
                    for q = 1:Q
                        if s_m(q) > 0
                            for n = 1:s_m(q)
                                l = L(m, s, q);
                                k = K(m, s, q);
                                data_c = [data_c; t m n_m W_t(m) R_t(m) X_t(m) l k PiV(m, s, q) + W_t(m)*l + R_t(m)*k];
                            end
                        end
                    end
                end
            end
            data_c(1, :) = [];
            data_sim = [data_sim; data_c];
        end
    end
    data_sim(1, :) = [];
end

function A_draw = draw_firms(mu, sigma, NBS, E, Q, Qv, T)
    
    A_draw = cell(2, 1);
    A = cell(T, NBS);
    A_hist = cell(T, NBS);
    Aq_idx = cell(T, NBS);
    
    for t = 1:T
        for bs = 1:NBS
            A{t, bs} = random('Lognormal', mu, sigma, [E(t), 1]);
            A_hist{t, bs} = zeros(Q, 1);
            Aq_idx{t, bs} = zeros(E(t), 1);
        end
    end
    
    for t = 1:T
        for nbs = 1:NBS
            Aq = zeros(E(t), 1);
            for j = 1:E(t)
                q = 1;
                qj = 0;
                while qj == 0
                    if q < Q
                        A_j = A{t, nbs}(j);
                        if (A_j - Qv(q) > 0) & (A_j - Qv(q + 1) <= 0)
                            d1 = A_j - Qv(q);
                            d2 = Qv(q) - A_j;
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
                Aq_idx{t, nbs}(j) = qj;
                Aq(j) = Qv(qj);
            end
            for q = 1:Q
                flags = (Aq_idx{t, nbs} == q);
                A_hist{t, nbs}(q) = sum(flags);
            end
        end
    end

    A_draw{1} = A_hist;
    A_draw{2} = Aq_idx;
end

function [cp, resnorm, residual, exitflag, output] = nfp(t, M, Q, States, ...
        A_hist, PiV, S, nbs)
    cp_mat = ones(M, Q)/(M*Q);
    cp = cp_mat(:);
    lb = zeros(M, Q);
    ub = ones(M, Q);
    sparse = zeros(M*Q, M*Q);
    k = 1;
    for q = 1:Q
        for m = 1:M
            for i = 1:Q
                for l = 1:M
                    if (m == l) | (q == i)
                        sparse(k) = 1;
                    end
                    k = k + 1;
                end
            end
        end
    end
    options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective', ...
            'JacobPattern', sparse, 'UseParallel', true);
    At_hist = A_hist{t, nbs}(:);
    f = @(z) fpe(z, M, Q, States, At_hist, PiV, S);
    [cp, resnorm, residual, exitflag, output] = lsqnonlin(f, cp, ...
        lb, ub, options);
end

function z = fpe(cp, M, Q, States, At_hist, PiV, S)
    cp_mat = reshape(cp, [M, Q]);
    z = choice_prob(exp_profit(cp_mat, States, At_hist, PiV, M, Q, S), ...
        M, Q) - cp_mat;
end

function cp = choice_prob(EPi, M, Q)
    cp = zeros(M, Q);
    exp_EPi = exp(EPi);
    for q = 1:Q
        denom = sum(exp_EPi(:, q));
        cp(:, q) = exp_EPi(:, q)/denom;
    end
end

function EPi = exp_profit(p, States, At_hist, PiV, M, Q, S) 
    EPi = zeros(M, Q);
    for q = 1:Q
        for m = 1:M
            p_s = zeros(S, 1);
            for s = 1:S
                p_s(s) = pr_state(p(m, :), s, States, At_hist, Q);
            end
            EPi(m, q) = PiV(m, :, q)*p_s;
        end
    end
end

function p_s = pr_state(pr, s, States, At_hist, Q) 
    p_s = 1;
    impossible = 0;
    for q = 1:Q
        n = At_hist(q);
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
