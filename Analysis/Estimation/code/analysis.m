clc;
clear;
close all;
echo off;

diary ../output/analysis.log
diary on;
rng(1);

%%% I. Initialize  %%%

%parpool(feature('20'));
load_as = '../temp/data.csv';
data = readmatrix(load_as);
T = 5;
E = 30;
M = 6;
Nq_max = 4;
alpha = 0.6;
mu_eps = 0;
beta_eps = 1;
mean_eps = mu_eps + beta_eps*0.57721;

mu = 1;
sigma = 1;
gamma = 1;
eta = 1;

Q = 6;
Q_band = 2;
Qv = zeros(Q, 1);
for q = 1:Q
    Qv(q + 1) = Q_band*q;
end

theta = [mu; sigma; gamma; eta];
llh = loglikelihood(theta, data, T, E, M, Q, Qv, Nq_max, alpha, mu_eps, beta_eps, mean_eps);
disp(llh)

function llh = loglikelihood(theta, data, T, E, M, Q, Qv, Nq_max, alpha, mu_eps, beta_eps, mean_eps)
    mu = theta(1);
    sigma = theta(2);
    gamma = theta(3);
    eta = theta(4);

    llh = 0;
    for cl = 1:T
        trows = (data(:, 1) == cl);
        data_t = data(trows, :);
        N_t = sum(trows);
        [~, idx] = unique(data_t(:, 2));
        data_mt = data_t(idx, :);
        W = data_t(:, 4);
        R = data_t(:, 5);
        X = data_t(:, 6);
        N_m = data_t(:, 3);
        N_mt = zeros(M, 1);
        for m = 1:M
            m_idx = (data_t(:,2) == m);
            N_mt(m) = sum(m_idx);
        end 
        disp('Checks');
        assert(sum(N_mt) == N_t);
        N_m(idx)
        N_mt
        assert(sum(N_m(idx) ~= N_mt) == 0);

        %%% II. Draw Characteristics %%%
        A = zeros(N_t, 1);
        A_cdf = round(N_t*cdf('Lognormal', Qv, mu, sigma)); 
        A_hist = zeros(Q, 1);
        A_hist(1) = A_cdf(1);
        for q = 2:Q
            A_hist(q) = A_cdf(q) - A_cdf(q-1);
        end
        disp(A_hist)
        d = sum(A_hist) - N_t; % review later--eventually want E entrants
        disp(d)
        disp(E)
        if d > 0
            [~, idx] = max(A_hist, [], 'linear');
            A_hist(idx) = A_hist(idx) - d;
        elseif d < 0
            [~, idx] = min(A_hist, [], 'linear');
            A_hist(idx) = A_hist(idx) - d;
        end 
        disp('A_hist')
        disp(A_hist)
        Aq = zeros(E, 1);
        Aq_idx = zeros(E, 1);
        
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

               
        %%% IV. Calculate Cournot Payoffs in Each State %%%
        
        P = zeros(M, S);
        Y = zeros(M, S, Q);
        L = zeros(M, S, Q);
        K = zeros(M, S, Q);
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
        
        for m = 1:M
            L(m, :, :) = (1-alpha)*Y(m, :, :)/W(m);
            K(m, :, :) = (alpha)*Y(m, :, :)/R(m);
        end

        VStates = cell(M, 1);
        VS_sz = zeros(1, M); 
        for m = 1:M
            k = 1;
            for s = 1:S
                if (N_mt(m) == sum(States(s, :))) & ...
                  (sum(States(s,:) <= A_hist') == Q)
                    VStates{m}(k) = s;
                    k = k + 1;
                end
            end
            VS_sz(m) = k - 1;
        end
        disp(N_mt(m));
        disp('Size of state space');
        disp(VS_sz);
        disp('Number of firms');
        disp(N_mt);
        disp('A_hist')
        disp(A_hist)

        disp('Begin diagnostics');
        tstate = ones(1, M);
        tstate(M) = 0; % do we need this?
        TStates = [];
        NTStates = [];
        maxed = zeros(1, M);
        St = 1;
        disp('Loop started');
        count = 0;
        big = [0 0];
        while sum(maxed) < M
            maxed = (tstate == VS_sz);
            % More aggressive skipping strategy
            if M > 4
                pairs = nchoosek(1:(M-3), 2);
                npairs = size(pairs, 1);
                row = 1;
                while (row <= npairs) & sum(big) == 0
                    pair = pairs(row, :);
                    s1 = VStates{pair(1)}(tstate(pair(1)));
                    n1 = States(s1, :);
                    s2 = VStates{pair(2)}(tstate(pair(2)));
                    n2 = States(s2, :);
                    if sum(n1 + n2 > A_hist') > 0
                        big = pairs(row, :);
                    end
                    row = row + 1;
                end
            end
            if sum(big) > 0
                disp('Bypassing extraneous cluster states ');
                disp(tstate);
                if maxed(big(2)) == 0
                    tstate(big(2)) = tstate(big(2)) + 1;
                else
                    l = big(2);
                    l0 = 0;
                    while l0 == 0 & l > 0
                        if maxed(l) == 1
                            l = l - 1;
                        else
                            l0 = l;
                        end
                    end
                    if l0 == 0 
                        tstate = VS_sz;
                    else
                        tstate(l0) = tstate(l0) + 1;
                        for l = (l0 + 1):M
                            tstate(l) = 1;
                        end
                    end
                end
                if tstate(M) == 0 
                    tstate(M) = 1;
                end
                disp(tstate);
            else
                if maxed(M) == 0
                    tstate(M) = tstate(M) + 1;
                else
                    l = M;
                    l0 = 0;
                    while l0 == 0 & l > 0
                        if maxed(l) == 1
                            l = l - 1;
                        else
                            l0 = l;
                        end
                    end
                    tstate(l0) = tstate(l0) + 1;
                    for l = (l0 + 1):M
                        tstate(l) = 1;
                    end
                end
            end
            ts = zeros(M, 1);
            nf = zeros(M, Q);
            for m = 1:M
                s = VStates{m}(tstate(m));
                ts(m) = s;
                nf(m, :) = States(s, :);
            end
            N_st = sum(nf, 1);
            match_N = (N_st == A_hist');
            count = count + 1;
            if mod(count, 100) == 0 
                disp('State example') 
                disp(nf) 
                disp(N_st) 
                disp(A_hist') 
                disp(match_N) 
            end
            if sum(match_N) == Q
                disp('Viable cluster state found');
                disp(nf);
                TStates(St, :) = ts';
                NTStates(St, :, :) = nf;
                St = St + 1;
            end
        end
        disp('Cluster state space successfully computed');
        disp(NTStates(1, :, :));

        %%% V. Fixed Point to Solve for Firm Beliefs%%%
        
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
            'JacobPattern', sparse, 'UseParallel', true);
        f = @(z) fpe(z, M, Q, States, A_hist, PiV, S);
        [cp, resnorm, residual, exitflag, output] = lsqnonlin(f, cp, ...
            lb, ub, options);
        
        %%% VI. Entry and Ex-Post Outcomes  %%%
        cp_mat = reshape(cp, [M, Q]);
        EPi = exp_profit(cp_mat, States, A_hist, PiV, M, Q, S);
        M_J = zeros(E, 1);
        S_M = zeros(M, Q);
        for j = 1:E
            q = Aq_idx(j);
            Pi_j = EPi(:, q) + epsilon(:, j);
            [Pi_mj, m_j] = max(Pi_j);
            M_J(j) = m_j;
            S_M(m_j, q) = S_M(m_j, q) + 1;
        end
        
        State_M = zeros(M, 1);
        for m = 1:M
            s_m = S_M(m, :);
            [~, s] = ismember(s_m, States, 'rows');
            if s == 0
                disp('Eqm state is outside state space');
                disp(s_m);
            else
                State_M(m) = s;
                for q = 1:Q
                    if s_m(q) > 0
                        for n = 1:s_m(q)
                            llh = llh + log(pr_state(cp(m, :), s, States, A_hist, Q));
                            %l = L(m, s, q);
                            %k = K(m, s, q);
                        end
                    end
                end
            end
        end
    end
end

function z = fpe(cp, M, Q, States, A_hist, PiV, S)
    cp_mat = reshape(cp, [M, Q]);
    z = choice_prob(exp_profit(cp_mat, States, A_hist, PiV, M, Q, S), ...
        M, Q) - cp_mat;
    disp(z);
end

function cp = choice_prob(EPi, M, Q)
    cp = zeros(M, Q);
    exp_EPi = exp(EPi);
    for q = 1:Q
        denom = sum(exp_EPi(:, q));
        cp(:, q) = exp_EPi(:, q)/denom;
    end
end

function EPi = exp_profit(p, States, A_hist, PiV, M, Q, S) 
    EPi = zeros(M, Q);
    for q = 1:Q
        for m = 1:M
            p_s = zeros(S, 1);
            for s = 1:S
                p_s(s) = pr_state(p(m, :), s, States, A_hist, Q);
            end
            EPi(m, q) = PiV(m, :, q)*p_s;
        end
    end
end

function p_s = pr_state(pr, s, States, A_hist, Q) 
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
