function SNR = compute_SNR(NSizeGrid, Nfft, N0, TxPower, pathLoss)
    fftOccupancy = 12*NSizeGrid/Nfft;
    SNRdB = (TxPower-30) - pathLoss - 10*log10(fftOccupancy) - 10*log10(2*N0^2);
    SNR = 10^(SNRdB/10);
end