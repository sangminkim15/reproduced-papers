function x = randl(m, n)

% Generates random numbers that follows Laplacian distribution

% I/O
% m   number of matrix rows
% n   number of matrix columns
% x   matrix with Laplacian distributed random numbers 
%     with mean mu = 0 and std sigma = 1 (columnwise)

u1 = rand(m, n);
u2 = rand(m, n);

x = (1 / sqrt(2)) * log(u1./u2);

end