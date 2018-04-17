function [ periodogram ] = center_periodogram( T, obs, freq_grid, a )
% compute the centralized (by a) periodogram using hann window
hann = @(t) (1-cos(2*pi*t/T))/2;
hann_ft = @(v) T/2*exp(-1i*pi*T*v).*sinc(T*v)./(1-(T*v.*(abs(v)~=1/T)).^2).*(abs(v)~=1/T)-T/4.*(abs(v)==1/T);
if length(obs)*length(freq_grid)<1000000000
    periodogram = 1/T*abs(transpose(sum(repmat(hann(obs),[length(freq_grid),1])...
        .*exp(-1i*2*pi*freq_grid'*obs),2))-a*hann_ft(freq_grid));
else
    periodogram=zeros([1,length(freq_grid)]);
    for jj=1:length(freq_grid)
        periodogram(jj)=1/T*abs(sum(hann(obs).*exp(-1i*2*pi*freq_grid(jj)*obs))-a*hann_ft(freq_grid(jj)));
    end
end

end

