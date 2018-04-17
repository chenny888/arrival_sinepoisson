%% The arrival data is input as a vector "obs", with obs(i) being the time stamp of i-th arrivals
% The time starts from 0, and the maximal time is T
freq_grid = (0:(1400*6))/1400;  % the frequency grid to compute the periodogram
a = length(obs)/T;    % average arrival, to centralize the periodogram
periodogram_window = center_periodogram(T, obs, freq_grid, a);  % windowed periodogram
tau=tau_simulate(max(periodogram_window), T, length(obs)/T, freq_grid); % compute the data-driven threshold
tau_574 = tau+(0.0574-0.0181)*max(periodogram_window);   % the threshold of the periodogram above which the frequencies are selected
% fitted_freq is the set of selected frequencies, a is the constant, c is the coefficient of the cosine, d is the coefficient of the sine
% the arrival rate is rate(t,fitted_freq,a,c,d)
[ fitted_freq, a, c, d ] = lse_time_cont( obs, periodogram_window, freq_grid, tau, T);

%% plot
% plot the periodogram
plot(freq_grid,periodogram_window);hold on;xlim([0.0,3.2]);ylim([0,4]);xlabel('\nu');ylabel('|H_c(\nu)|');
xL=get(gca,'XLim');line(xL,[tau_574,tau_574],'Color','r','LineStyle','--','LineWidth',1);
yL = get(gca,'YLim');
for jj=1:length(fitted_freq)
line([fitted_freq(jj) fitted_freq(jj)],yL,'Color','r','LineStyle','--','LineWidth',1);
end
