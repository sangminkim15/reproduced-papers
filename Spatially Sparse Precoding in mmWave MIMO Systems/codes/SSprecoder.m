function [FRF, FBB] = SSprecoder(Ns, NtRF, H, At)

% Spatially Sparse Precoding (SS Precoding)

% I/O
% Nt    # of transmit antennas
% NtRF  # of transmit RF chains
% FRF   RF Precoder (Nt * NtRF)
% FBB   Baseband Precoder (NtRF * Ns)

FRFnext = zeros(0,0);
Fo = Fopt(Ns, H);
Fres = Fo;

for i = 1 : NtRF
    FRFprev = FRFnext;
    PSI = At' * Fres;
    k = diag(PSI * PSI') == max(diag(PSI * PSI'));
    
    FRFnext = horzcat(FRFprev, At(:,k));
    FBB = (FRFnext' * FRFnext) \ (FRFnext' * Fo);
    
    Fres = (Fo - FRFnext * FBB) / norm(Fo - FRFnext * FBB, 'fro');
end

FRF = FRFnext;
FBB = sqrt(Ns) * FBB / norm(FRF * FBB, 'fro');

end