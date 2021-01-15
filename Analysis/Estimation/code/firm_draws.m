
Q = 4;
Qv = [1.5 3 4.5 6];
mu = 1;
sigma = 1;

parts = firms(5, mu, sigma, Qv, Q);
parts{1}
parts{2}

function parts = firms(e, mu, sigma, Qv, Q)
    probs = discrete_probs(mu, sigma, Qv, Q)
    A_hists = partition(e, Q);
    len = size(A_hists, 1);
    pr_part = [];
    for j = 1:len
        pr_part = [pr_part; prod(probs.^A_hists(j, :))];
    end
    parts = {A_hists, pr_part}
end

function probs = discrete_probs(mu, sigma, Qv, Q)
    probs = [cdf('Lognormal', Qv(1), mu, sigma)];
    for q = 2:(Q-1)
        probs = [probs, cdf('Lognormal', Qv(q), mu, sigma) - cdf('Lognormal', Qv(q-1), mu, sigma)];
    end
    probs = [probs, 1 - cdf('Lognormal', Qv(Q-1), mu, sigma)];
end

function list = partition(e, Q)
    x = nsumk(Q,e);
    len = size(x, 1);
    list = [];
    for j = 1:len
        perm = perms(x(j, :));
        perm = unique(perm, 'rows');
        list = [list; perm];
    end
end

function x = nsumk(n, k)
    m = nchoosek(k+n-1,n-1);
    dividers = [zeros(m,1),nchoosek((1:(k+n-1))',n-1),ones(m,1)*(k+n)];
    x = diff(dividers,1,2)-1;
end

