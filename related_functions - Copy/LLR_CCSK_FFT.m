function L  = LLR_CCSK_FFT(y, etaq, N, sigm)
L = zeros(size(y));
for n = 1 : N
    L(:,n) = ifft(fft(etaq).*conj(fft(y(:,n))));
    L(:,n)=L(:,n)-min(L(:,n));
    L(:,n)=L(:,n)*(2/sigm^2);
end
end

