function y = rate(t, freq, constant, cos_coef, sin_coef) 
% generate the rate at time t, with frequency, phase, and magnitude generated in the .m file
% vectorize the t argument
    y = constant+(cos_coef*cos(2*pi*freq'*t)+sin_coef*sin(2*pi*freq'*t));
end
