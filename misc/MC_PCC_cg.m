%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the Mass and Coriolis + Gravity Matrices
% given the current state and geometric parameters of the pushpuppet robot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M, C] = MC_4_cg(q,qd,m,r,L0,d)
% q = state
% qd = velocity
% qdd = acceleration
% m = mass of a module;
% r = radius of the hex plate;
% L0 = Length of a module;
% d = distance to cable;
% N = number of links (number of modules plus number of motors)

N = 2;
a_grav = [0;0;-9.81;0;0;0];

I = zeros(6,6,N);
inertia = [1/4*m*r^2; 1/4*m*r^2; 1/2*m*r^2];
I(:,:,1) = diag([m;m;m;inertia]);
I(:,:,2) = I(:,:,1);


Smod = zeros(6,3,N);
Xup = zeros(6,6,N);
v = zeros(6,N);
a = zeros(6,N);
f = zeros(6,N);
C = zeros(2*N,1);

Xtree = eye(6);
i = 1;
qi = 1:3;
[XJ, Smod(:,:,i),~,dJ] = PCC_jacobian(q(qi),d,L0,qd(qi));
Xup(:,:,1) = XJ*Xtree;
Xtree = eye(6);
vJ = Smod(:,:,i)*qd(qi);
v(:,i) = vJ;
a(:,i) = Xup(:,:,1)*(-a_grav) + spatial_cross(v(:,i))*vJ + dJ * qd(qi);
f(:,i) = I(:,:,i)*a(:,i) + -spatial_cross(v(:,i)).'*I(:,:,i)*v(:,i);

% Recursive Newton Euler to Calculate C+G
for i = 2:N
    qi = (i-1)*3+1:3*i;
    [XJ, Smod(:,:,i),~,dJ] = PCC_jacobian(q(qi),d,L0,qd(qi));
    Xup(:,:,i) = XJ*Xtree;
    Xtree = eye(6);
    vJ = Smod(:,:,i)*qd(qi);
    v(:,i) = Xup(:,:,i)*v(:,i-1) + vJ;
    a(:,i) = Xup(:,:,i)*a(:,i-1) + dJ * qd(qi) + spatial_cross(v(:,i))*vJ;
    f(:,i) = I(:,:,i)*a(:,i) + -spatial_cross(v(:,i)).'*I(:,:,i)*v(:,i);
end

for i = N:-1:1
    qi = (i-1)*3+1:3*i;
    C(qi,1) = 2 * Smod(:,:,i)' * f(:,i);
    if i ~= 1
        f(:,i-1) = f(:,i-1) + Xup(:,:,i)'*f(:,i);
    end
end

% Composite Rigid Body Algorithm to calculate M
IC = I;				% composite inertia calculation

for i = N:-1:1
    if i ~= 1
        IC(:,:,i-1) = IC(:,:,i-1) + Xup(:,:,i)'*IC(:,:,i)*Xup(:,:,i);
    end
end

H = zeros(N*3);

for i = 1:N
    qi = (i-1)*3+1:3*i;
    fh3 = IC(:,:,i) * Smod(:,:,i);
    H(qi,qi) = Smod(:,:,i)' * fh3;
    j = i;
    while j < N+1 && j > 1
        fh3 = Xup(:,:,j)' * fh3;
        j = j - 1;
        qj = (j-1)*3+1:3*j;
        H(qi,qj) = fh3' * Smod(:,:,j); %(S{j}' * fh).';
        H(qj,qi) = Smod(:,:,j)' * fh3;
    end
end

M = H;
end