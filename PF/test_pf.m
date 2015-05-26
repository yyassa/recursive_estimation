


% number of runs:
n_runs = 8;
result = zeros(1,n_runs);
time   = zeros(1,n_runs);

for i=1:n_runs
rng(i);
[result(i),time(i)] = run(0,0);
countdown = n_runs - i
end

error_mean = mean(result)
error_var = cov(result)
