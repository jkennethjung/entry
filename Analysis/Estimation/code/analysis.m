clc;
clear;
close all;
echo off;

diary ../output/analysis.log
diary on;
rng(1);

%%% I. Initialize  %%%

parpool(24);
load_as = '../temp/small.csv';
save_as = '../output/estimates.csv';
NS = 100;
mu = 20;
alpha = 0.4;
gamma = 40;
eta = 1;
Nq_max = 9;

mu_eps = 0;
beta_eps = 1;
mean_eps = mu_eps + beta_eps*0.57721;
Q = 4;
Qv = [10; 20; 30; 40];
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
    t_nozero = (data(:, 3) > 0);
    E(t) = sum(t_rows & t_nozero);
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

A_hists = {};
for t = 1:T
    A_hists{end+1} = nsumk(Q, E(t));
end

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
n_obs = size(data, 1);
Z = [data(:, 3:6) ones(n_obs, 1)];
y = data(:, 7);

data_mt = collapse_data(data, T, M);
nm_obs = size(data_mt, 1);
nm = data_mt(:, 3);
Zm = [data_mt(:, 4:6) ones(nm_obs, 1)];

Beta_1 = [mean(y); var(y)];
Beta_2 = (Z.'*Z)^(-1)*(Z.'*y);
Beta_3 = (Zm.'*Zm)^(-1)*(Zm.'*nm);
Beta_0 = [Beta_1; Beta_2; Beta_3];
disp(Beta_0);

theta = [gamma; eta; mu];
c = clock;
fix(c)
aux = @(theta) auxiliary(theta, Beta_0, NS, alpha, A_hists, M, E, W, R, ...
        X, Q, Qv, T, S, States, epsilon, ncol);
[theta, dBeta, exitflag, output] = patternsearch(aux, theta, [], [], [], [], ...
    theta*1e-2, theta*1e2);
disp('Outer loop optimization finished');
disp(theta);
c = clock;
fix(c)
writematrix(theta, save_as);

function dBeta = auxiliary(theta, Beta_0, NS, alpha, A_hists, M, E, W, R, X, Q, ...
        Qv, T, S, States, epsilon, ncol)
    disp('theta:')
    disp(theta)
    gamma = theta(1);
    eta = theta(2);
    mu = theta(3);
    pr_hist = {};
    for t = 1:T
        pr_hist{end+1} = pr_firms(A_hists{t}, E(t), mu, Qv, Q);
    end
    data_sim = simulate(theta, A_hists, pr_hist, NS, alpha, M, E, W, R, X, Q, Qv, T, S, ...
        States, epsilon, ncol);
    Betas = [];
    for ns = 1:NS
        data_ns = data_sim{ns};
        n_obs = size(data_ns, 1);
        Z = [data_ns(:, 3:6) ones(n_obs, 1)];
        y = data_ns(:, 7);

        data_mt = collapse_data(data_ns, T, M);
        nm_obs = size(data_mt, 1);
        nm = data_mt(:, 3);
        Zm = [data_mt(:, 4:6) ones(nm_obs, 1)];

        Beta_1 = [mean(y); var(y)];
        Beta_2 = (Z.'*Z)^(-1)*(Z.'*y);
        Beta_3 = (Zm.'*Zm)^(-1)*(Zm.'*nm);
        Beta = [Beta_1; Beta_2; Beta_3];
        Betas = [Betas Beta];
    end
    Beta = mean(Betas, 2);
    disp('Auxiliary model:');
    disp(Beta);
    dBeta = (Beta_0 - Beta).'*(Beta_0 - Beta);
    disp('Objective function');
    disp(dBeta);
end

function data_sim = simulate(theta, A_hists, pr_hist, NS, alpha, M, E, W, R, X, Q, Qv, T, S, ...
        States, epsilon, ncol)
    gamma = theta(1);
    eta = theta(2);
    mu = theta(3);

    A_draw = draw_firms(mu, NS, E, Q, Qv, T);
    Aq_idx = A_draw{2};
    data_sim = cell(NS, 1);
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
                            a = Qv(q);
                            if Iq(q) == 1
                                Y(m, s, q) = y(j);
                                [lab, cap] =  factor_demand(y(j), a, w, r, alpha);
                                L(m, s, q) = lab;
                                K(m, s, q) = cap;
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
           
        %%% V. Fixed Point to Solve for Firm Beliefs%%%
    
        At_hist = A_hists{t}; 
        prt_hist = pr_hist{t};
        [cp, resnorm, residual, exitflag, output] = nfp(t, M_t, Q, States, ...
            At_hist, prt_hist, PiV, S);
        cp_mat = reshape(cp, [M_t, Q]);
        EPi = exp_profit(cp_mat, States, At_hist, prt_hist, PiV, M_t, Q, S);
        for ns = 1:NS
            %%% VI. Entry and Ex-Post Outcomes  %%%
            M_J = zeros(E_t, 1);
            S_M = zeros(M_t, Q);
            for j = 1:E_t
                q = Aq_idx{t, ns}(j);
                shocks = reshape(epsilon{t, ns}(:, j), [M_t, 1]);
                Pi_j = EPi(:, q) + shocks ;
                [Pi_mj, m_j] = max(Pi_j);
                M_J(j) = m_j;
                S_M(m_j, q) = S_M(m_j, q) + 1;
            end
            
            State_M = zeros(M_t, 1);
            data_c = [];
            for m = 1:M_t
                s_m = S_M(m, :);
                n_m = sum(s_m);
                [~, s] = ismember(s_m, States, 'rows');
                if s == 0
                    disp('Warning: outside state space');
                elseif s == 1
                    data_c = [data_c; t m 0 W_t(m) R_t(m) X_t(m) 0 0 0];
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
            data_sim{ns} = [data_sim{ns}; data_c];
        end
    end
end

function A_draw = draw_firms(mu, NS, E, Q, Qv, T)
    
    A_draw = cell(2, 1);
    A = cell(T, NS);
    A_hist = cell(T, NS);
    Aq_idx = cell(T, NS);
    
    for t = 1:T
        for bs = 1:NS
            A{t, bs} = random('Exponential', mu, [E(t), 1]);
            A_hist{t, bs} = zeros(Q, 1);
            Aq_idx{t, bs} = zeros(E(t), 1);
        end
    end
    
    for t = 1:T
        for ns = 1:NS
            Aq = zeros(E(t), 1);
            for j = 1:E(t)
                q = 1;
                qj = 0;
                while qj == 0
                    if q < Q
                        A_j = A{t, ns}(j);
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
                Aq_idx{t, ns}(j) = qj;
                Aq(j) = Qv(qj);
            end
            for q = 1:Q
                flags = (Aq_idx{t, ns} == q);
                A_hist{t, ns}(q) = sum(flags);
            end
        end
    end

    A_draw{1} = A_hist;
    A_draw{2} = Aq_idx;
end

function [cp, resnorm, residual, exitflag, output] = nfp(t, M, Q, States, ...
        At_hist, prt_hist, PiV, S)
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
            'JacobPattern', sparse, 'UseParallel', false);
    f = @(z) fpe(z, M, Q, States, At_hist, prt_hist, PiV, S);
    [cp, resnorm, residual, exitflag, output] = lsqnonlin(f, cp, ...
        lb, ub, options);
end

function z = fpe(cp, M, Q, States, At_hist, prt_hist, PiV, S)
    cp_mat = reshape(cp, [M, Q]);
    z = choice_prob(exp_profit(cp_mat, States, At_hist, prt_hist, PiV, M, Q, S), ...
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

function EPi = exp_profit(p, States, At_hist, prt_hist, PiV, M, Q, S) 
    EPi = zeros(M, Q);
    len = size(At_hist, 1);
    p_sm = zeros(S, M);
    parfor s = 1:S
        for m = 1:M
            pr_states = zeros(len, 1);
            for i = 1:len 
                hist = At_hist(i, :);
                pr_states(i) = pr_state(p(m, :), s, States, hist, Q);
            end
            p_s(s, m) = prt_hist.' * pr_states;
        end
    end
    for m = 1:M
        p_s = p_sm(:, m);
        for q = 1:Q
            EPi(m, q) = PiV(m, :, q)*p_s;
        end
    end
end

function p_s = pr_state(pr, s, States, At_hist, Q) 
    p_s = 1;
    state = States(s, :);
    impossible = (state > At_hist);
    if sum(impossible) > 0
        p_s = 0;
    else
        for q = 1:Q
            n = At_hist(q);
            k = States(s, q);
            pr_q = pr(q);
            p_s = p_s*nchoosek(n, k)*pr_q^k*(1-pr_q)^(n-k); 
        end
    end
end

function mc = marginal_cost(a, w, r, alpha)
    [l, k] = factor_demand(1, a, w, r, alpha);
    mc = w*l + r*k;
end 

function [l, k] = factor_demand(y, a, w, r, alpha)
    l = y*((1-alpha)*r/(alpha*w))^alpha/a;
    k = y*(alpha*w/((1-alpha)*r))^(1-alpha)/a;
end 

function pr_part = pr_firms(At_hists, e, mu, Qv, Q)
    probs = discrete_probs(mu, Qv, Q);
    len = size(At_hists, 1);
    pr_part = [];
    for j = 1:len
        pr_part = [pr_part; prod(probs.^At_hists(j, :))];
    end
end

function probs = discrete_probs(mu, Qv, Q)
    probs = [cdf('Exponential', Qv(1), mu)];
    for q = 2:(Q-1)
        probs = [probs, cdf('Exponential', Qv(q), mu) - cdf('Exponential', Qv(q-1), mu)];
    end
    probs = [probs, 1 - cdf('Exponential', Qv(Q-1), mu)];
end

function x = nsumk(n, k)
    m = nchoosek(k+n-1,n-1);
    dividers = [zeros(m,1),nchoosek((1:(k+n-1))',n-1),ones(m,1)*(k+n)];
    x = diff(dividers,1,2)-1;
end

function collapsed = collapse_data(data, T, M)
    collapsed = [];
    for t = 1:T
        for m = 1:M(t)
            mt_rows  = ((data(:, 2) == m) & (data(:, 1) == t));
            data_mt = data(mt_rows, :);
            collapsed = [collapsed; data_mt(1, :)];
        end
    end
end
