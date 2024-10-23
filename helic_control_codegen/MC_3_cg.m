%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the Mass and Coriolis + Gravity Matrices
% given the current state and geometric parameters of the pushpuppet robot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M, C] = MC_3_cg(q,qd,m,r,L0,d)
% q = state
% qd = velocity
% qdd = acceleration
% m = mass of a module;
% mm = mass of a motor;
% hm = height of a motor;
% rm = "radius" a motor;
% r = radius of the hex plate;
% L0 = Length of a module;
% d = distance to cable;
% N = number of links (number of modules plus number of motors)

N = 4;
a_grav = [0;0;9.81;0;0;0];

% qi = {1, 2:4, 5, 6:8, 9, 10:12, 13, 14:16};

I = zeros(6,6,4);
inertia = [1/4*m*r^2; 1/4*m*r^2; 1/2*m*r^2];
I(:,:,1) = diag([m;m;m;inertia]);
I(:,:,2) = I(:,:,1);
I(:,:,3) = I(:,:,1);
I(:,:,4) = I(:,:,1);


Smod = zeros(6,3,3);
Smotz = [0;0;0;0;0;1];
Smoty = [0;0;0;0;1;0];
Xup = zeros(6,6,4);
v = zeros(6,4);
a = zeros(6,4);
f = zeros(6,4);
C = zeros(10,1);

Xtree = eye(6); %[diag(ones(3,1)) [0;0;0]; 0 0 0 1];
i = 1;
qi = 1;
g = [cos(q(qi)) -sin(q(qi)) 0 0; sin(q(qi)) cos(q(qi)) 0 0; 0 0 1 0; 0 0 0 1];
XJ = adj_calc(g,1);
Xup(:,:,1) = XJ*Xtree;
vJ = Smotz*qd(qi);
v(:,i) = vJ;
a(:,i) = Xup(:,:,1)*(-a_grav) + spatial_cross(v(:,i))*vJ;
f(:,i) = I(:,:,i)*a(:,i) + -spatial_cross(v(:,i)).'*I(:,:,i)*v(:,i);

% Recursive Newton Euler to Calculate C+G
for i = 2:N
    qi = (i-2)*3 + 2:(i-2)*3 + 4;
    [XJ, Smod(:,:,i-1),~,dJ] = PCC_jacobian(q(qi),d,L0,qd(qi));
    Xup(:,:,i) = XJ*Xtree;
    Xtree = eye(6);
    vJ = Smod(:,:,i-1)*qd(qi);

    v(:,i) = Xup(:,:,i)*v(:,i-1) + vJ;
    a(:,i) = Xup(:,:,i)*a(:,i-1) + dJ * qd(qi) + spatial_cross(v(:,i))*vJ;
    f(:,i) = I(:,:,i)*a(:,i) + -spatial_cross(v(:,i)).'*I(:,:,i)*v(:,i);
end

for i = N:-1:1
    if i == 1
        qi = 1;
        C(qi,1) = 2 * Smotz' * f(:,i);
    else
        qi = (i-2)*3 + 2:(i-2)*3 + 4;
        C(qi,1) = 2 * Smod(:,:,i-1)' * f(:,i);
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

H = zeros(10);
% fh3 = zeros(6,3);
% fh1 = zeros(6,1);

for i = 1:N
    if i == 1
        fh1 = IC(:,:,i) * Smotz;
        H(i,i) = Smotz' * fh1;
    else
        qi = (i-2)*3 + 2:(i-2)*3 + 4;
        fh3 = IC(:,:,i) * Smod(:,:,i-1);
        H(qi,qi) = Smod(:,:,i-1)' * fh3;
        j = i;
        while j < N+1 && j > 1
            fh3 = Xup(:,:,j)' * fh3;
            j = j - 1;
            if j == 1
                H(qi,j) = fh3' * Smotz; %(S{j}' * fh).';
                H(j,qi) = Smotz' * fh3;
            else
                qj = (j-2)*3 + 2:(j-2)*3 + 4;
                H(qi,qj) = fh3' * Smod(:,:,j-1); %(S{j}' * fh).';
                H(qj,qi) = Smod(:,:,j-1)' * fh3;
            end
        end
    end
end


% 
% for i = 1:N
%     if mod(i,2) == 0
%         qi = 2 + (2*i-4):4 + (2*i-4);
%         fh3 = IC(:,:,i) * Smod(:,:,i/2);
%         H(qi,qi) = Smod(:,:,i/2)' * fh3;
%         j = i;
%         while j < N+1 && j > 1
%             fh3 = Xup(:,:,j)' * fh3;
%             j = j - 1;
%             if mod(j,2) == 0
%                 qj = 2 + (2*j-4):4 + (2*j-4);
%                 H(qi,qj) = fh3' * Smod(:,:,j/2); %(S{j}' * fh).';
%                 H(qj,qi) = Smod(:,:,j/2)' * fh3;
%             elseif j == 1
%                 qj = 1 + (2*j-2);
%                 H(qi,qj) = fh3' * Smotz; %(S{j}' * fh).';
%                 H(qj,qi) = Smotz' * fh3;
%             else
%                 qj = 1 + (2*j-2);
%                 H(qi,qj) = fh3' * Smoty; %(S{j}' * fh).';
%                 H(qj,qi) = Smoty' * fh3;
%             end
%         end
%     else
%         qi = 1 + (2*i-2);
%         if i == 1
%             fh1 = IC(:,:,i) * Smotz;
%             H(qi,qi) = Smotz' * fh1;
%         else
%             fh1 = IC(:,:,i) * Smoty;
%             H(qi,qi) = Smoty' * fh1;
%         end
% 
%         j = i;
%         while j < N+1 && j > 1
%             fh1 = Xup(:,:,j)' * fh1;
%             j = j - 1;
%             if mod(j,2) == 0
%                 qj = 2 + (2*j-4):4 + (2*j-4);
%                 H(qi,qj) = fh1' * Smod(:,:,j/2); %(S{j}' * fh).';
%                 H(qj,qi) = Smod(:,:,j/2)' * fh1;
%             elseif j == 1
%                 qj = 1 + (2*j-2);
%                 H(qi,qj) = fh1' * Smotz; %(S{j}' * fh).';
%                 H(qj,qi) = Smotz' * fh1;
%             else
%                 qj = 1 + (2*j-2);
%                 H(qi,qj) = fh1' * Smoty; %(S{j}' * fh).';
%                 H(qj,qi) = Smoty' * fh1;
%             end
%         end
%     end
% 
% 
% end

M = H;
end