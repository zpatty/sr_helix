function [tau, tau_r, x, M, C, A, cq] = helix_controller(q,dq,qd,dqd,ddqd,d,m,r,kb,ks,bb,bs,bm,L0,Kp,KD,Kpx, KDx, xd, dxd, dxr, conv_pcc,conv_motor) %#codegen


qp = zeros(3,3);
qp(:,1) = q(2:4);
qp(:,2) = q(5:7);
qp(:,3) = q(8:10);

[M, C] = MC_3_cg(q,dq,m,r,L0,d);

K = diag([0; repmat([kb;kb;ks],3,1)]);
D = diag([bm; repmat([bb;bb;bs],3,1)]);

% Fci = cell(1,4);
% Ai = cell(1,4);
Ai = zeros(3,3,3);
cq = zeros(3,1);
for i = 1:length(qp)
    dx = qp(1,i);
    dy = qp(2,i);
    dL = qp(3,i);

    del = sqrt(dx^2+dy^2);
    theta = del/d;

    Dq = del - sin(del);
    L = dL + L0;

    if dx < 10e-6 && dy < 10e-6
        c = (L0+dL)/3;
        Aq = [0, -1, 0; 1, 0, 0; 0, 0, 1];
    else
        c = 2*((L0 + dL)/(theta) - d)*sin(theta/6);

        Aq = [dx * dy * Dq / del^3, (- dx^2 * del - dy^2 * sin(del)) / del^3, dx * Dq * L / del^3;
            (dy^2 * del + dx^2 * sin(del)) / del^3, - dx * dy * Dq / del^3, dy * Dq * L / del^3;
            0, 0, sin(del) / del];
    end
    cq(i) = c;

    Al = [d*cosd(60) d*cosd(60) -d;
    -d*cosd(30) d*cosd(30) 0;
    1 1 1];

    At = [1,0,0;0,1,0;0,0,1]*-1;
    Ai(:,:,i) = Aq * At * Al;
end
A = blkdiag(1,Ai(:,:,1),Ai(:,:,2),Ai(:,:,3));
conversion = diag([conv_motor, repmat([conv_pcc, conv_pcc, conv_pcc],1,3)]);
tau = conversion*(A\(C + M*ddqd + K*qd + D*dqd + Kp*(qd-q) + KD*(dqd - dq)));

Kpr = Kp;
Kpr(1,1) = 0;
KDr = KD;
KDr(1,1) = 0;
[J, x] = J_r(q,L0,d);
Lam = eye(3)/(J * (M \ J'));
Jbar = M \ J' * Lam;
tau_r = conversion*(A\(C + J'*Jbar'*(K*q + D*dq) + J'*Lam*(Kpx*(xd-x) + KDx*(dxd - dxr))) + (eye(10) - J'*Jbar')*(-Kpr * q - KDr * dq));
end

