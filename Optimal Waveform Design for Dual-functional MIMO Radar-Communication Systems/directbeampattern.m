function [Rd] = directbeampattern (p)

cvx_solver sedumi
cvx_begin quiet

variable b
variable Rd(p.N, p.N) hermitian semidefinite

expression u(length(p.theta))

for idx = 1 : length(p.theta)
    a = zeros(p.N, 1);
    
    for jdx = 1 : p.N
        a(jdx, 1) = exp(1i * pi * (jdx - ceil(p.N / 2)) * sin(p.theta(idx)));
    end
    
    u(idx) = (b * p.Pd_theta(idx) - a' * Rd * a);
end

minimize norm(u, 2)

subject to
trace(Rd) == p.Pt;
cvx_end
    
end