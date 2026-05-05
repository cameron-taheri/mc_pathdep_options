% Monte Carlo Methods for Path-Dependent Options Pricing

% compute the expected price of a simulated stock in six months

clear; clc;

% ------------------------------------------------------------------------
% Set Parameters
% ------------------------------------------------------------------------
S0 = 31;                % current stock price
K = 29;                 % strike price
B = 27;                 % barrier value
r = 0.03;               % (constant) risk-free rate
mu = 0.1537;            % (objective) drift
sigma = 0.1809;         % volatility
num_sims = 100000;      % number of simulations
N = 180;                % length of period in days (6 months, 30-day months)
dt = 1/360;             % time step in years
T = N * dt;             % total time in years (0.5 years)

% ------------------------------------------------------------------------
% Simulate Terminal Stock Prices (100,000 sims.)
% ------------------------------------------------------------------------
% simulation of terminal stock prices using the exact solution to gbm sde
% S(t) = S(0) * exp[(mu - 0.5 * sigma^2) * t + sigma * sqrt(t) * e]
% where e is a standard normal random variable
% dW(t) is discretized as 

rng(0); % seed for reproducibility

% initialize the brownian motion vector
dWT = sqrt(T) * randn(num_sims, 1); % generate standard normal random variables

% Svec will be vector of daily prices during life of the option
Svec = S0 * exp((mu - 0.5 * sigma^2) * T + sigma * dWT);

% calculate the expected value at t = T
ES6m = mean(Svec);

% output result
disp(['Expec. value of S at the end of 6-month period: ', num2str(ES6m)])

%% compute the value of the barrier option

% ------------------------------------------------------------------------
% Estimate Value of Barrier Option
% ------------------------------------------------------------------------
% 6M knock-in barrier option (Euro down-and-in put option) with strike = 29
% option is negotiated today (t=0) when stock price is 31
% option's life begins at t = 1 (number of days = 180)
% if S(t) <= 27 for any t, option becomes exercisable at expiration (t=T)
% payout = max(K - S(T),0) if barrier condition satisfied, or 0

rng(0); % set seed for reproducibility

% initialize vector of stock prices and payouts
payouts = zeros(num_sims, 1);

for i = 1:num_sims
    % generate random (Gaussian) shocks
    dW = randn(N, 1) * sqrt(dt);
    
    % tomorrow's stock price (first value of Svec) with risk-neutral drift
    Svec = zeros(N, 1);
    Svec(1) = S0 * exp((r - 0.5 * sigma^2) * dt + sigma * dW(1));
    
    % for each simulation, loop over each day and compute stock price path
    for t = 2:N
        Svec(t) = Svec(t-1) * exp((r - 0.5 * sigma^2) * dt + sigma * dW(t));
    end

    % payout of the barrier option
    if any(Svec <= B)
        payouts(i) = max(K - Svec(end), 0); % if barrier is breached
    else
        payouts(i) = 0; % no payout if barrier is not breached
    end

end


% ------------------------------------------------------------------------
% Vbarrier: compute today's (t=0) price of the barrier option & 95% CIs
% ------------------------------------------------------------------------
discounted_payouts = exp(-r * T) * payouts;
SE = std(discounted_payouts) / sqrt(num_sims);
Vbarrier = mean(discounted_payouts);

% compute the 95% confidence interval (critical value ~1.96)
VLB = Vbarrier - 1.96 * SE; % lower bound
VUB = Vbarrier + 1.96 * SE; % upper bound

% ------------------------------------------------------------------------
% Display results
% ------------------------------------------------------------------------
disp(['Price of this barrier option: ', num2str(Vbarrier)]);
CI = [VLB, VUB];
disp(['95 pct confidence interval: ', num2str(CI)]);
disp(['length of 95 pct CI: ', num2str(VUB-VLB)])

% plot the price path, barrier price, and strik price
% Generate time vector for plotting
time = (0:N-1) * dt;

% Plot the stock price path for the last simulation
figure;
plot(time, Svec, 'b', 'LineWidth', 1);
hold on;
yline(B, 'r--', 'Barrier', 'LineWidth', 1.2);
yline(K, 'g--', 'Strike', 'LineWidth', 1.2);
xlabel('Time (Years)');
ylabel('Stock Price');
title('Stock Price Path with Barrier and Strike Price');
legend('Stock Price', 'Barrier', 'Strike Price');
grid on;
hold off;

%% value of the barrier option using antithetic variable technique

% ------------------------------------------------------------------------
% Estimate Value of Barrier Option
% ------------------------------------------------------------------------
% 6M knock-in barrier option (Euro down-and-in put option) with strike = 29
% option is negotiated today (t=0) when stock price is 31
% option's life begins at t = 1 (number of days = 180)
% if S(t) <= 27 for any t, option becomes exercisable at expiration (t=T)
% payout = max(K - S(T),0) if barrier condition satisfied, or 0

rng(0); % set seed for reproducibility

% initialize vector of payouts
payouts = zeros(num_sims, 1);

for i = 1:num_sims
    
    % generate random (Gaussian) shocks
    dW = randn(N, 1) * sqrt(dt);
    
    % initialize stock paths
    Svecpos = zeros(N, 1); 
    Svecneg = zeros(N, 1); 
    
    % tomorrow's stock price (first value of Svec) with risk-neutral drift
    % positive shock
    Svecpos(1) = S0 * exp((r - 0.5 * sigma^2) * dt + sigma * dW(1));
    
    % negative shock
    Svecneg(1) = S0 * exp((r - 0.5 * sigma^2) * dt - sigma * dW(1));

    % for each simulation, loop over each day and compute stock price path
    for t = 2:N
        
        % positive path
        Svecpos(t) = Svecpos(t-1) * exp((r - 0.5 * sigma^2) * dt + sigma * dW(t));
        
        % negative path
        Svecneg(t) = Svecneg(t-1) * exp((r - 0.5 * sigma^2) * dt - sigma * dW(t));

    end

    % payout of the barrier option (positive path)
    if any(Svecpos <= B)
        payoutspos = max(K - Svecpos(end), 0); % if barrier is breached
    else
        payoutspos = 0; % no payout if barrier is not breached
    end
    
    % payout of the barrier option (negative path)
    if any(Svecneg <= B)
        payoutsneg = max(K - Svecneg(end), 0); % if barrier is breached
    else
        payoutsneg = 0; % no payout if barrier is not breached
    end

    % mean of payouts (under antithetic technique)
    payouts(i) = (payoutsneg + payoutspos) / 2;

end

% ------------------------------------------------------------------------
% Vbarrier: compute today's (t=0) price of the barrier option & 95% CIs
% ------------------------------------------------------------------------
discounted_payouts = exp(-r * T) * payouts;
SE = std(discounted_payouts) / sqrt(num_sims);
Vbarrieranti = mean(discounted_payouts);

% compute the 95% confidence interval (critical value ~1.96)
VLBanti = Vbarrieranti - 1.96 * SE; % lower bound
VUBanti = Vbarrieranti + 1.96 * SE; % upper bound

% ------------------------------------------------------------------------
% Display results
% ------------------------------------------------------------------------
disp(['Price of this barrier option (antithetic): ', num2str(Vbarrieranti)]);
CIanti = [VLBanti, VUBanti];
disp(['95 pct confidence interval (antithetic): ', num2str(CIanti)]);
disp(['length of 95 pct CI (antithetic): ', num2str(VUBanti-VLBanti)])

% plot the price path, barrier price, and strik price
% Generate time vector for plotting
time = (0:N-1) * dt;

% Plot the stock price paths (positive and negative) for the last simulation
figure;
plot(time, Svecpos, 'b', 'DisplayName', 'Positive Path');
hold on;
plot(time, Svecneg, 'r', 'DisplayName', 'Negative Path');
yline(B, 'k--', 'Barrier', 'LineWidth', 1.5);
yline(K, 'g--', 'Strike Price', 'LineWidth', 1.5);
xlabel('Time (Years)');
ylabel('Stock Price');
title('Stock Price Paths for Barrier Option');
subtitle('Antithetic Variable Technique');
legend show;
grid on;

%% risk-free short-term interest rate r(t) under Vasicek

clear; clc;

% ------------------------------------------------------------------------
% Set Parameters
% ------------------------------------------------------------------------
r0 = 0.01;              % current ST RF Rate
k = 0.5;                % speed of mean reversion
r_ss = 0.03;            % steady-state ST RF Rate
sigma_r = 0.01;         % volatility
T = 10;                 % number of years
dt = 1/360;             % time step in years
N = T / dt;             % total time in years (3600 days)

% ------------------------------------------------------------------------
% Euler-Maruyama Discretization
% ------------------------------------------------------------------------
% continuous Vasicek: dr(t) = k(r_ss - r(t)*dt + sigma*dW(t)
% discrete (E-M) Vasicek: r(t+dt) = r(t) + k(r_ss - r(t))*dt + sigma*dW
% where dW is defined as before to be dW = sqrt(dt) * randn
% simulating one path for the ST RF Rate

r = zeros(1, N);  % initialize vector of interest rates
r(1) = r0;        % set initial value of r

rng(123); % set seed

% loop over all periods
for i = 1:N-1
    % brownian motion increment
    dW = sqrt(dt) * randn;

    % Euler-Maruyama Vasicek model
    r(i+1) = r(i) + k * (r_ss - r(i)) * dt + sigma_r * dW;
end

% plot r along with r_ss and a horizontal line at r = 0
% plot the interest rate path
time_vec = (0:N-1) * dt;

figure;
plot(time_vec, r, 'b', 'LineWidth', 1);
hold on;
yline(r_ss, 'r--', 'Steady-State Rate'); % add steady-state line
yline(0, 'g--');
xlabel('Time (Years)');
ylabel('Interest Rate');
title('Simulation of Short-Term Risk-Free Rate');
subtitle('Vasicek Model (Euler-Maruyama Discr.)');
grid on;
hold off;

%% EV of risk-free short-term interest rate r(t) under Vasicek in 2yrs
clear; clc;

% ------------------------------------------------------------------------
% Set Parameters
% ------------------------------------------------------------------------
r0 = 0.01;              % initial short-term risk-free rate
k = 0.5;                % speed of mean reversion
r_ss = 0.03;            % steady-state (long-run mean) ST RF rate
sigma_r = 0.01;         % volatility
T = 2;                  % simulation horizon in years
dt = 1/360;             % time step in years
N = T / dt;             % number of time steps (720)
num_sims = 100000;      % number of simulations

% ------------------------------------------------------------------------
% Monte Carlo with antithetic variables technique
% ------------------------------------------------------------------------
rng(123); % set seed

% store terminal values of r for each simulation pair
r_terminal = zeros(num_sims, 1);

for i = 1:num_sims
    
    % generate Brownian increments for one path
    dW = sqrt(dt) * randn(N, 1);
    
    % initialize both paths at r0
    rpos = zeros(N, 1);
    rneg = zeros(N, 1);

    % set the first realization to r0
    rpos(1) = r0;
    rneg(1) = r0;
    
    % Euler-Maruyama with positive and negative paths (antithetic)
    for t = 1:N-1
        
        rpos(t+1) = rpos(t) + k * (r_ss - rpos(t)) * dt + sigma_r * dW(t);
        rneg(t+1) = rneg(t) + k * (r_ss - rneg(t)) * dt - sigma_r * dW(t);
    
    end
    
    % average terminal values from both paths
    r_terminal(i) = (rpos(end) + rneg(end)) / 2;
    
end

% ------------------------------------------------------------------------
% Compute EV of r 2yrs ahead
% ------------------------------------------------------------------------
Er2y = mean(r_terminal);

% ------------------------------------------------------------------------
% Display results
% ------------------------------------------------------------------------
disp(['The expected value r two years from now is: ', num2str(Er2y)])
