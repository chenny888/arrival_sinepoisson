function [ obs ] = generate_data(T, freq, phase, mag, const_term ) % using phase and mag, not cos and sin
% generate the arrival process given freq, phase, mag, const_term
%   for the formula refers to the rate function. basically the rate is
%   mag*cos(freq*t+phase)
%   NOT the complex version
ub=const_term+sum(abs(mag)); % upper bound of the rate
cos_coef = cos(phase).*mag;
sin_coef = -sin(phase).*mag;
tt = 0;
obs = [];
while tt<T
    tt = tt+exprnd(1/ub);
    %     if unifrnd(0,1) < 1.3*exp(rate(tt,freq,const_term,cos_coef,sin_coef))/ub
    if unifrnd(0,1) < (rate(tt,freq,const_term,cos_coef,sin_coef))/ub
        obs = [obs, tt];
    end
end
obs = obs(1:(length(obs)-1));
end

