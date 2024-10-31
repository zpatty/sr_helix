%% Pushpuppets Multi Segment Lagrangian Cleaned Up
clear

syms d L0 hm hb dc
syms dx1 dy1 dL1
syms ddx1 ddy1 ddL1
syms dx2 dy2 dL2 
syms ddx2 ddy2 ddL2 
syms dx3 dy3 dL3 
syms ddx3 ddy3 ddL3
syms th
syms dth

qp1 = [dx1 dy1 dL1].';
qp2 = [dx2 dy2 dL2].';
qp3 = [dx3 dy3 dL3].';


dqp1 = [ddx1 ddy1 ddL1].';
dqp2 = [ddx2 ddy2 ddL2].';
dqp3 = [ddx3 ddy3 ddL3].';


q = [th; qp1; qp2; qp3];
dq = [dth; dqp1; dqp2; dqp3];


Rq = @(dx,dy,d,del) [1 + dx^2/del^2*(cos(del/d) - 1), dx*dy/del^2*(cos(del/d) - 1), -dx/del*sin(del/d);
    dx*dy/del^2*(cos(del/d) - 1), 1+dy^2/del^2*(cos(del/d) - 1), -dy/del*sin(del/d);
    dx/del*sin(del/d), dy/del*sin(del/d), cos(del/d)].';
tq = @(dx,dy,dL,del,L0,d) d*(L0 + dL)/del^2*[dx*(1 - cos(del/d)); dy*(1 - cos(del/d)); del*sin(del/d)];

gpcc = @(dx,dy,dL,L0,d) [Rq(dx,dy,d,sqrt(dx^2 + dy^2)) tq(dx,dy,dL,sqrt(dx^2 + dy^2),L0,d); 0 0 0 1];


Ry = @(x) [cos(x), 0, sin(x);
        0, 1, 0
        -sin(x), 0, cos(x)];

Rz = @(x) [cos(x), -sin(x), 0;
    sin(x), cos(x), 0;
    0, 0, 1];

% del = sqrt(dx^2+dy^2);
% syms del
g = {};
g{1} = [Rz(th), [0; 0; 0]; 0 0 0 1];
g12 = gpcc(qp1(1),qp1(2),qp1(3),-L0,d);
g23 = gpcc(qp2(1),qp2(2),qp2(3),-L0,d);
g3e = gpcc(qp3(1),qp3(2),qp3(3),-L0,d);
g{2} = simplify(g{1}*gpcc(qp1(1),qp1(2),qp1(3),-L0/2,d));
g{3} = simplify(g{1}*g12*gpcc(qp2(1),qp2(2),qp2(3),-L0/2,d));
g{4} = simplify(g{1}*g12*g23*gpcc(qp3(1),qp3(2),qp3(3),-L0/2,d));
g_ee = simplify(g{1}*g12*g23*g3e*[eye(3) Rz(0)*[0.0; 0; 0]; 0 0 0 1]);

% gc = sym(zeros(6,1));
% syms phi
% syms dc
% temp = g{1}*[eye(3) Rz(0)*[dc; 0; 0]; 0 0 0 1];
% gc = simplify(temp(3,4));
% 
% 
% J = simplify(jacobian(gc,[dx;dy;dL]));
% dJ = simplify([jacobian(J(1),qp1)*dqp1, jacobian(J(2),qp1)*dqp1, jacobian(J(3),qp1)*dqp1]);

% temp = simplify(g{2}*[eye(3) Rz(0)*[dc; 0; 0]; 0 0 0 1]);
g_x = simplify(g_ee(1:3,4));

Jr = simplify(jacobian(g_x,q));
vpa(subs(Jr, [d; L0; q; dc], [0.01; 0.1; 0.000000000001*ones(10,1); 0.001]))
% dJr = simplify([jacobian(J_r(1),q)*dq, jacobian(J_r(2),q)*dq, jacobian(J_r(3),q)*dq]);
% matlabFunction(Jr,"File","J_r_new");
% matlabFunction(dJr,"File","dJr");

% limit(limit(J_r,dx1,0),dy1,0)
% limit(limit(dJ,dx,0),dy,0)
% 
% limit(limit(J,dx,0),dy,0)
% limit(limit(dJ,dx,0),dy,0)

%%
ginvf = @(g) [g(1:3,1:3).' -g(1:3,1:3).'*g(1:3,4); zeros(1,3) 1];
ginv = {};
for i=1:length(g)
    ginv{i} = ginvf(g{i});
end

%%
dg = {};
for i=1:length(g)
    for j=1:length(q)
        dg{i,j} = simplify(diff(g{i},q(j)));
    end
end

%%
J = {};
for i=1:length(g)
    Jg = [];
    for j=1:length(q)
        Jg = [Jg (vee(ginv{i}*dg{i,j}))];
    end
    J{i} = (Jg);
end


%%

% J0 = limit(limit(subs(Jxq,del,sqrt(dx^2+dy^2)),dx,0),dy,0);
% g0 = limit(limit(subs(T_q,del,sqrt(dx^2+dy^2)),dx,0),dy,0);

%% Mass, Coriolis, and Gravity

syms I1 [3,1] real
syms I2 [3,1] real
syms m [2,1] real

I = {I1, I2};

Mi = {};
M = zeros(length(q));
for i = 1:length(g)
    if i == 1
        mi = m(1);
        Ii = I1;
    else
        mi = m(2);
        Ii = I2;
    end
    Mi{i} = diag([mi; mi; mi; Ii]);
    M = M + (J{i}.'*Mi{i}*J{i});
end
%%
dM = sym(zeros(length(q)));
for i = 1:length(q)
    for j = 1:length(q)
        dM(i,j) = jacobian(M(i,j),q)*dq;
    end
end
fid = fopen('dM.txt','w');
fprintf(fid,'%s',char(dM));
fclose(fid);
%%
C = zeros(length(M));
for k = 1:length(M)
    for i = 1:length(M)
        for j = 1:length(M)
            christoffel(i,j,k) = (diff(M(i,j),q(k))+diff(M(i,k),q(j))-diff(M(k,j),q(i)));
%             christoffel0(i,j,k) = (diff(M0(i,j),q(k))+diff(M0(i,k),q(j))-diff(M0(k,j),q(i)));
        end
    end
C = C + christoffel(:,:,k)*dq(k);
% C0 = C0 + christoffel0(:,:,k)*qd(k);
end

matlabFunction(C,"File","C_calc")
% C0 = limit(limit(C,dx,0),dy,0)

% C = simplify(subs(C,sqrt(dx^2+dy^2),delta));
% M = subs(M,sqrt(dx^2+dy^2),delta);
%%
h = [g{1}(3,4); g{2}(3,4); g{3}(3,4); g{4}(3,4)];
syms gr
Vg = [m1; m2; m2; m2]'*gr*h;
Fg = jacobian(Vg,q).';
%%
mv = 0.03;
r = 0.089/2;
L0v = 0.09;

a_grav = [0;0;-9.81;0;0;0];

qv = [0.0001;0.0;0.0;0.0001;0.0;-0.0];
dqv = [0;0.1;-0.1;0.01; 0.2; 0.0];
qdd = [0;0;0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
qi = {1, 2:4, 5, 6:8, 9, 10:12};
dv = 0.04;
inertia_values = [1/4*mv*r^2; 1/4*mv*r^2; 1/2*mv*r^2; 
                1/4*mv*r^2; 1/4*mv*r^2; 1/2*mv*r^2];
G2pcc_test = double(subs(Fg,[m;L0;q;dq;d;I1;I2;gr],[mv;mv;L0v;qv;dqv;dv;inertia_values;9.81]))
[Mt, Ct, G] = MC_PCC_cg(qv,dqv,mv,r,L0v,dv)
M2pcc_test = double(subs(M,[m;L0;q;dq;d;I1;I2],[mv;mv;L0v;qv;dqv;dv;inertia_values]))
C2pcc_test = double(subs(C,[m;L0;q;dq;d;I1;I2],[mv;mv;L0v;qv;dqv;dv;inertia_values]))*dqv

%%
fid = fopen('M_2.txt','w');
fprintf(fid,'%s',char(M));
fclose(fid);
fid = fopen('C_2.txt','w');
fprintf(fid,'%s',char(C));
fclose(fid);

%%
z = tq(3);
G_x = subs(diff(z,dx),del,delta);
G_y = subs(diff(z,dy),del,delta);
G_L = subs(diff(z,dL),del,delta);

limit(limit(diff(z,dx),dx,0),dy,0)
limit(limit(diff(z,dy),dx,0),dy,0)
limit(limit(diff(z,dL),dx,0),dy,0)

%% Acuation
Aa = [-cos(phi)*sin(theta), -sin(phi)*sin(theta), 0;
    -sin(phi), cos(phi), (L0+dL)*((theta) - sin(theta)/theta^2);
    0, 0, sin(theta)/theta];
simplify(Jaq.'*subs(Aa, [phi, theta], [acos(dx/delta), delta/d]))