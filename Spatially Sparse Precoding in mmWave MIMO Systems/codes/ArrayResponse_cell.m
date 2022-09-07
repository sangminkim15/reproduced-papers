function [Acell] = ArrayResponse_cell(N, Ncl, Nray, std)

% Generates Tx and Rx Array Response Vector Matrix

% I/O
% N     square root of # of transmit antennas (Nt = N^2)
% Ncl   # of clusters
% Nray  # of rays per cluster

Acell = cell(1, Ncl);
PHI = zeros(1, Nray);
THETA = zeros(1, Nray);

for i = 1 : Ncl
    j = 1;
    k = 1;
    
    phimean = (2*pi) * rand(1,1);
    thetamean = (pi) * rand(1,1);
    
    while j <= Nray
        phi = phimean + (std) * randl(1,1);             
        
        if phi >= phimean - pi/6 && phi <= phimean + pi/6          % azimuth angle range = 60 deg
            PHI(j) = phi;
            j = j + 1;
        end
    end
    
    while k <= Nray
        theta = thetamean + (std) * randl(1,1);    
        
        if theta >= thetamean - pi/18 && theta <= thetamean + pi/18   % elevation angle range = 20 deg
            THETA(k) = theta;
            k = k + 1;
        end
    end
    
    Amat = zeros(N.^2, Nray);
    
    for l = 1 : Nray
        Amat(:,l) = ArrayResponse_vec(N, PHI(l), THETA(l));
    end
    
    Acell{1,i} = Amat;
    
end

end