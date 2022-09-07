function [] = Convergence ()

u1 = ones(1,4);
u2 = ones(1,8);
Etx = 1;

iteration = zeros(1,31);
WSR1 = zeros(1000,31);
WSR2 = zeros(1000,31);

for q = 1 : 1000
    H1 = channel(4, 1, 4, 10);

    [K, ~] = size(H1);

    B = cell(1,K);
    tr = 0;
    for i = 1 : K
        B{1,i} = H1{i,1}';
        tr = tr + trace(B{1,i} * B{1,i}');
    end

    for i = 1 : K
        B{1,i} = (Etx / tr).^0.5 * B{1,i};
    end

    [~,Q] = size(B{1,1});

    A = cell(K,K);
    W = cell(K,K);
    E = cell(1,K);
    R = zeros(1,K);
    
    n = 0;
    while n <= 30
        for i = 1 : K
            [E{1,i}, ~] = WMMSE(u1(i), H1{i,1}, B, i);
            R(i) = log(det(E{1,i}^(-1))) / log(2);
        end
        iteration(n+1) = n;
        WSR1(q,n+1) = sum(u1.*R);
        
        for j = 1 : K
            for l = 1 : K
                if j == l
                    A{j,l} = RxFilter(H1{j,1}, B, j);
                    [~, W{j,l}] = WMMSE(u1(j), H1{j,1}, B, j);
                else
                    A{j,l} = zeros(Q,Q);
                    W{j,l} = zeros(Q,Q);
                end
            end
        end
    
        Hmatrix = cell2mat(H1);
        Amatrix = cell2mat(A);
        Wmatrix = cell2mat(W);
    
        Bmatrix = Beamforming(Hmatrix, Amatrix, Wmatrix, Etx);
    
        B = matrixdiv(Bmatrix, K);
        
        n = n + 1;
    end
end

for q = 1 : 1000
    H2 = channel(4, 1, 8, 10);

    [K, ~] = size(H2);

    B = cell(1,K);
    tr = 0;
    for i = 1 : K
        B{1,i} = H2{i,1}';
        tr = tr + trace(B{1,i} * B{1,i}');
    end

    for i = 1 : K
        B{1,i} = (Etx / tr).^0.5 * B{1,i};
    end

    [~,Q] = size(B{1,1});

    A = cell(K,K);
    W = cell(K,K);
    E = cell(1,K);
    R = zeros(1,K);
    
    n = 0;
    while n <= 30
        for i = 1 : K
            [E{1,i}, ~] = WMMSE(u2(i), H2{i,1}, B, i);
            R(i) = log(det(E{1,i}^(-1))) / log(2);
        end
        iteration(n+1) = n;
        WSR2(q,n+1) = sum(u2.*R);
        
        for j = 1 : K
            for l = 1 : K
                if j == l
                    A{j,l} = RxFilter(H2{j,1}, B, j);
                    [~, W{j,l}] = WMMSE(u2(j), H2{j,1}, B, j);
                else
                    A{j,l} = zeros(Q,Q);
                    W{j,l} = zeros(Q,Q);
                end
            end
        end
    
        Hmatrix = cell2mat(H2);
        Amatrix = cell2mat(A);
        Wmatrix = cell2mat(W);
    
        Bmatrix = Beamforming(Hmatrix, Amatrix, Wmatrix, Etx);
    
        B = matrixdiv(Bmatrix, K);
        
        n = n + 1;
    end
end

WSRmean1 = mean(WSR1);
WSRmean2 = mean(WSR2);

plot(iteration, WSRmean1, iteration, WSRmean2, 'Marker', '*', 'LineWidth', 2);
title('Convergence Properties');
xlabel('Iteration Number');
ylabel('Sum-Rate');
ylim([0 15]);
grid on;

end