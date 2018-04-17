function [ fitted_freq, a, c, d ] = lse_time_cont(obs, periodogram, freq_grid, tau, T)
% estimate the frequencies and magnitudes using LSE in the time domain
% a is the constant term, c and d are the cos and sin coefficients.
% need to bucket the arrivals
% hann = @(t) sin(pi*t/T).^4; % another window
% hann_ft = @(v) (3*T/2)*exp(-1i*pi*T*v).*sinc(T*v)./(4-5*(T*v).^2+(T*v).^4);
% hann = @(t) (1-cos(2*pi*t/T))/2;
% hann_ft = @(v) T/2*exp(-1i*pi*T*v).*sinc(T*v)./(1-(T*v.*(abs(v)~=1/T)).^2).*(abs(v)~=1/T)-T/4.*(abs(v)==1/T);
% find frequencies first
ind = find((periodogram>[periodogram(2:end),0]) & (periodogram>[0,periodogram(1:(end-1))]) & freq_grid>3/T & periodogram>tau);
[B,I]=sort(periodogram(ind),'descend');
I = ind(I);ind = [];
while ~isempty(I)
    ind=[ind,I(1)];
    I = I(abs(freq_grid(I)-freq_grid(I(1)))>3/T);
end
fitted_freq = freq_grid(ind);
freq_num = length(fitted_freq);
% below adding extra frequencies around the selected frequencies for the least squares fit.
% I believe that this will improve the prediction performance in the setting where \eps(T) \notto 0
% extras=[];
% for j=1:freq_num
%     extras=[extras,(fitted_freq(j)-1.5/T):(3/(2*log(T)*T)):(fitted_freq(j)+1.5/T)]
% end
% fitted_freq=[fitted_freq,extras]
% freq_num=length(fitted_freq)

% LSE in time domain continuously
% freq_double = [0, fitted_freq, -fitted_freq];
% xty=zeros([length(freq_double),1]); % x^T*y vector
% for j=1:length(freq_double)
%     xty(j) = sum(exp(freq_double(j)*2*pi*obs*1i));
% end
%
% xtx=zeros([length(freq_double), length(freq_double)]); % X^T*X
% for i1 = 1:length(freq_double)
%     for i2 = (i1):length(freq_double)
%         xtx(i1,i2) = (exp((freq_double(i1)+freq_double(i2))*2*pi*T*1i)-1)/2/pi/(freq_double(i1)+freq_double(i2))/1i;
%         xtx(i2,i1) = xtx(i1,i2);
%     end
% end
% xtx(1,1) = T;
% for i1 = 2:(freq_num+1)
%     xtx(i1, i1+(freq_num)) = T;
%     xtx(i1+(freq_num), i1) = T;
% end

% complex version
freq_double = [0, fitted_freq, -fitted_freq];
xty=zeros([length(freq_double),1]); % x^T*y vector
for j=1:length(freq_double)
    xty(j) = sum(exp(-freq_double(j)*2*pi*obs*1i));
end
xtx=zeros([length(freq_double), length(freq_double)]); % X^T*X
for i1 = 1:length(freq_double)
    for i2=1:length(freq_double)
        myeta=freq_double(i1)-freq_double(i2);
        xtx(i1,i2)=T*exp(-1i*pi*T*myeta)*sinc(T*myeta);
    end
end
% above complex version

coef_est = xtx\xty; %%sahand comment: change this to matlab version of solve
% a = abs(coef_est(1)); c=(2*real(coef_est(2:(1+freq_num))))'; d=-2*(imag(coef_est(2:(1+freq_num))))';
a = abs(coef_est(1)); c=(2*real(coef_est(2:(1+freq_num))))'; d=-2*(imag(coef_est(2:(1+freq_num))))';
end

