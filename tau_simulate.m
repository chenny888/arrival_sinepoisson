function [ threshold ] = tau_simulate( max_h, T, lambda, freq_grid )
% compute the threshold tau by simulating an NPP
obs = generate_data(T, 0, 0, 0, lambda);
p = center_periodogram(T, obs, freq_grid, lambda);
threshold = 0.0181*max_h + 1.06*max(p); %0.018 corresponds to a separation of 3/T
end
