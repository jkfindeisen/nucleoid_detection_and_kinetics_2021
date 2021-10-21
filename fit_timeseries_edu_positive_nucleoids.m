function fit_timeseries_edu_positive_nucleoids()
% This script fits the measured time series of EdU-positive nucleoids
% with different kinetic models and computes important model constants and
% their errors. Reproduces values given Supplementary Table 1 and Fig. 1f in the
% manuscript.
%
% Required that the nucleoids' positions are determined before and combined
% to give EdU active fraction. Here simply stored in a csv file in the data
% folder.
%
% Time internally is in hours [0-120], displayed also in days
%
% Part of "The TFAM to mtDNA ratio defines inner-cellular nucleoid
% populations with distinct activity levels"
%
% Jan Keller-Findeisen, Dep. NanoBiophotonics, MPI Biophysical Chemsitry,
% Göttingen, Germany

%% parameters
rng(0);
Nbootstrp = 400; % number of repetitions of the fit for bootstraping

% fit options
options = optimoptions('fmincon', 'Display', 'off', 'OptimalityTolerance', 1e-8);

% load experimental data
m = readmatrix('data/edu_positive_nucleoids_timeseries.csv');
data.x = m(:, 1);
data.y = m(:, 2);
data.err = m(:, 3);

x = (0:120).'; % for display a fine grid for the first 120 hours
x2 = (0:0.01:30)*24; % compute model for a long time (30 days) to be able to calculate 95% positive EdU time interval

ci_level = 0.05; % 95% quantiles for the confidence interval level

% average effective growth rate
d = log(2) / (7*24); % division rate set to 7 days
fprintf('division rate: %.2f days\n',d*24);

%% Model A - One subpopulation (fit only first 18h)
fprintf('Model A - one subpopulation\n');

% free model parameter: growth rate
model = @(p, t) 1 - exp(-2*p(1)*t);
ix = data.x <= 18;
minimizer = @(x) sum((model(x, data.x(ix)) - data.y(ix)).^2); % least squares estimation
p0 = 0.05; % initial value
lb = d; % lower bounds (growth rates must at least be the division rate)
ub = Inf;
% minimization
p = fmincon(minimizer, p0, [], [], [], [], lb, ub, [], options);
yA = model(p, x);

%% Model B - Independency model fit
fprintf('Model B - Independency model\n\n');

% free model parameters: fraction in slow pop. (gamma), slow growth rate
% (alpha_s), fast growth rate (alpha_f) corresponding to Eq. 7 in the
% manuscript
model = @(p, t) 1 - p(1)*exp(-2*p(2)*t) -(1-p(1))*exp(-2*p(3)*t);
minimizer = @(x) sum((model(x, data.x) - data.y).^2); % least squares estimation
p0 = [0.5, 0.005, 0.05]; % initial values
lb = [0, d, d]; % lower bounds (both growth rates must at least be the division rate)
ub = [1, Inf, Inf];
% minimization
p = fmincon(minimizer, p0, [], [], [], [], lb, ub, [], options);
yB = model(p, x);

% now the error part
ex = zeros(Nbootstrp, numel(p));
ey = zeros(Nbootstrp, numel(data.y)); % values of all the random models at the data points
ey2 = zeros(Nbootstrp, numel(x2)); % values of all the random models for the first 30 hours
for i = 1 : Nbootstrp
    ye = draw_pseudo_data(data.y, data.err);
    minimizer = @(x) sum((model(x, data.x) - ye).^2); % LSE
    pi = fmincon(minimizer, p0, [], [], [], [], lb, ub, [], options);
    ex(i, :) = pi;
    ey(i, :) = model(pi, data.x);
    ey2(i, :) = model(pi, x2);
end

ci = quantile(ex, [ci_level/2,1-ci_level/2]);
ytt = quantile(ey2, [ci_level/2, 0.5, 1-ci_level/2]);

fprintf('Fast replication (active) population\n');
fprintf(' Size %.0f%% (%.0f%% - %.0f%%)\n', (1-p(1))*100, (1-ci(2, 1))*100, (1-ci(1, 1))*100);
fprintf(' Replication rate (per nucleoid, per day) %.2f (%.2f - %.2f)\n', p(3) * 24, ci(1, 3) * 24, ci(2, 3) * 24);
fprintf(' Degradation rate (per nucleoid, per day) %.2f (%.2f - %.2f)\n', (p(3) - d) * 24, (ci(1, 3) - d) * 24, (ci(2, 3) - d) * 24);

fprintf('\nSlow replication (inactive) population\n');
fprintf(' Size %.0f%% (%.0f%% - %.0f%%)\n', p(1)*100, ci(1, 1)*100, ci(2, 1)*100);
fprintf(' Replication rate (per nucleoid, per day) %.2f (%.2f - %.2f)\n', p(2) * 24, ci(1, 2) * 24, ci(2, 2) * 24);
fprintf(' Degradation rate (per nucleoid, per day) %.2f (%.2f - %.2f)\n', (p(2) - d) * 24, (ci(1, 2) - d) * 24, (ci(2, 2) - d) * 24);
fprintf('\n\n');

% time until 50%/90% EdU positive
tis = [0.5, 0.9];
for i = 1 : numel(tis)
    ti = tis(i);
    
    ti1 = log(1 - ti) / (-2 * p(3) * 24);
    ti2 = log(1 - ti) / (-2 * ci(2, 3) * 24);
    ti3 = log(1 - ti) / (-2 * ci(1, 3) * 24);
    fprintf('Time until %.1f%% of the fast replicating (active) population is EdU-labeled: %.2f days (%.2f - %.2f days)\n', ti*100, ti1, ti2, ti3);
    
    ti1 = log(1 - ti) / (-2 * p(2) * 24);
    ti2 = log(1 - ti) / (-2 * ci(2, 2) * 24);
    ti3 = log(1 - ti) / (-2 * ci(1, 2) * 24);
    fprintf('Time until %.1f%% of the slow replicating (inactive) population is EdU-labeled: %.2f days (%.2f - %.2f days)\n', ti*100, ti1, ti2, ti3);
    
    [~, idx] = min(abs(ytt(2, :) - ti)); ti1 = x2(idx) / 24;
    [~, idx] = min(abs(ytt(3, :) - ti)); ti2 = x2(idx) / 24;
    [~, idx] = min(abs(ytt(1, :) - ti)); ti3 = x2(idx) / 24;
    fprintf('Time until %.1f%% of all nucleoids are EdU-labeled: %.2f days  (%.2f - %.2f days)\n', ti*100, ti1, ti2, ti3);
end

%% Model C - Inactivation model fit
fprintf('\n\nModel C - Inactivation model\n\n');

% equivalent to the complete model with alpha_s = 0, delta_s = 0 (Eq. 8)
% two free parameters Delta_f (=Alpha_f-Beta_f), Alpha_f
model = @(p, t) 1 - ((p(1)-d)*2*p(2)*exp(-d*t)+d*(2*p(2)-p(1))*exp(-2*p(2)*t))/(p(1)*(2*p(2)-d));
minimizer = @(x) sum((model(x, data.x) - data.y).^2); % least squares estimation
p0 = [0.005, 0.03]; % initial values
lb = [d, d]; % both must be at least the division rate
ub = [Inf, Inf];
p = fmincon(minimizer, p0, [], [], [], [], lb, ub, [], options);
yC = model(p, x);
% derived parameter values
alpha_f = p(2);
beta_f = p(2) - p(1);
tau = p(1) - d; % transfer rate fast to slow is Delta_f - d
gamma = tau / p(1); % fast growing fraction

% now the error part
ex = zeros(Nbootstrp, numel(p));
exf = zeros(Nbootstrp, 4);
ey = zeros(Nbootstrp, numel(data.y));
ey2 = zeros(Nbootstrp, numel(x2));
eeysp = zeros(Nbootstrp, numel(x2));
model_sp = @(x, t) 1 - ((x(4) * (x(1) + x(2)) + x(3))*exp(-d*t)-x(3)*(1-x(4))*exp(-(x(1)+x(2)+x(3)+d)*t))/x(4)/(x(1)+x(2)+x(3)); % model for the slow replicating ones only with parameter x = alpha_f, beta_f, tau, gamma
for i = 1 : Nbootstrp
    ye = draw_pseudo_data(data.y, data.err);
    minimizer = @(x) sum((model(x, data.x) - ye).^2); % LSE
    pi = fmincon(minimizer, p0, [], [], [], [], lb, ub, [], options);
    ex(i, :) = pi;
    xif = [pi(2), pi(2)-pi(1), pi(1) - d, (pi(1) - d) / pi(1)]; % alpha_f, beta_f, tau, gamma
    exf(i, :) = xif;
    ey(i, :) = model(pi, data.x);
    ey2(i, :) = model(pi, x2);
    eeysp(i, :) = model_sp(xif, x2);
end

cif = quantile(exf, [ci_level/2,1-ci_level/2]);
ytt = quantile(ey2, [ci_level/2, 0.5, 1-ci_level/2]);
yttsp = quantile(eeysp, [ci_level/2, 0.5, 1-ci_level/2]);

fprintf('Fast replication (active) population\n');
fprintf(' Size %.0f%% (%.0f%% - %.0f%%)\n', gamma*100, cif(1, 4)*100,  cif(2, 4)*100);
fprintf(' Replication rate (per nucleoid, per day) %.2f (%.2f - %.2f)\n', alpha_f * 24, cif(1, 1) * 24, cif(2, 1) * 24);
fprintf(' Degradation rate (per nucleoid, per day) %.2f (%.2f - %.2f)\n', beta_f * 24, cif(1, 2) * 24, cif(2, 2) * 24);

fprintf('\nSlow replication (inactive) population\n');
fprintf(' Size %.0f%% (%.0f%% - %.0f%%)\n', (1-gamma)*100, (1-cif(2, 4))*100, (1-cif(1, 4))*100);

fprintf('\nTransfer from fast to slow\n');
fprintf(' (per nucleoid per day) %.2f (%.2f - %.2f)\n', tau * 24, cif(1, 3) * 24, cif(2, 3) * 24);
fprintf('\n\n');

% fprintf('alpha_f %.2f/d (%.2f - %.2f) / day\n', alpha_f * 24, cif(1, 1) * 24, cif(2, 1) * 24);
% fprintf('beta_f %.2f/d (%.2f - %.2f) / day\n', beta_f * 24, cif(1, 2) * 24, cif(2, 2) * 24);
% fprintf('tau %.2f/d (%.2f - %.2f) / day\n', tau * 24, cif(1, 3) * 24, cif(2, 3) * 24);
% fprintf('gamma %.2f/d (%.2f - %.2f)\n', gamma, cif(1, 3), cif(2, 4));
% fprintf('1-gamma %.2f/d (%.2f - %.2f)\n', 1-gamma, 1-cif(1, 4), 1-cif(2, 3));

% time until 50% labeled, 90% labeled
tx = [0.5, 0.9];
for i = 1 : 2
    ti = tx(i);
    
    ti1 = log(1 - ti) / (-2 * alpha_f * 24);
    ti2 = log(1 - ti) / (-2 * cif(2, 1) * 24);
    ti3 = log(1 - ti) / (-2 * cif(1, 1) * 24);
    fprintf('Time until %.1f%% of the fast replicating (active) population is EdU-labeled: %.2f days (%.2f - %.2f days)\n', ti*100, ti1, ti2, ti3);
    
    [~, idx] = min(abs(yttsp(2, :) - ti)); ti1 = x2(idx) / 24;
    [~, idx] = min(abs(yttsp(3, :) - ti)); ti2 = x2(idx) / 24;
    [~, idx] = min(abs(yttsp(1, :) - ti)); ti3 = x2(idx) / 24;
    fprintf('Time until %.1f%% of the slow replicating (inactive) population is EdU-labeled: %.2f days (%.2f - %.2f days)\n', ti*100, ti1, ti2, ti3);
    
    [~, idx] = min(abs(ytt(2, :) - ti)); ti1 = x2(idx) / 24;
    [~, idx] = min(abs(ytt(3, :) - ti)); ti2 = x2(idx) / 24;
    [~, idx] = min(abs(ytt(1, :) - ti)); ti3 = x2(idx) / 24;
    fprintf('Time until %.1f%% of all nucleoids are EdU-labeled: %.2f days  (%.2f - %.2f days)\n', ti*100, ti1, ti2, ti3);
end

%% recreate Figure 1f from the manuscript
lw = 1.5;
fig = figure();
hold on;
plot(x, yA*100, 'k--', 'LineWidth', lw);
plot(x, yB*100, 'r-', 'LineWidth', lw);
plot(x, yC*100, 'b-', 'LineWidth', lw);
errorbar(data.x, data.y*100, data.err*100, 'ko');
xlabel('EdU-Incubation [h]');
xlim([0, 120]);
xticks(0:20:120);
ylabel('EdU-positive nucleoids (%)');
ylim([0, 100]);
yticks(0:20:100);
box on;
grid on;
title('reproduction of Fig. 1f');


end

function ye = draw_pseudo_data(y, e)
% draws random pseudo data given a mean (measured mean value) and an error
% assuming that the measured data would follow a normal distribution with
% this mean and standard deviation (used for artificial data generation in
% the bootstraping procedure)

ye = y + randn(size(e)) .* e;
% ye = max(ye, 0); % not necessary because LSE

end