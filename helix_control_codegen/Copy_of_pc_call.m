%% simulate puppet segment
clear
close all

L0 = 0.19;



N = 4;
q0 = zeros(1 + (N-1)*3,1);
dq0 = zeros(1 + (N-1)*3,1);
r = 0.03;
d = 0.0395;
kb = 0;
ks = 0;
bb = 0;
bs = 0;
bm = 0;
e = L0/1000;
segment_mass = 0.071;
med_plates = 0.036;
end_plate = 0.023;
m = 2*segment_mass + 2 * med_plates + end_plate;
kps = 0;
kpb = 0;
km = 0;
Kp = diag([km, repmat([kpb, kpb, kps],1,N-1)]);
Dm = 0.0;
Ds = 0;
KD = diag([Dm, repmat([Ds, Ds, Ds],1,N-1)]);
qhist = [q0];
tauhist = [];
chist = [];
Kpx = 50.0;
KDx = 50.0;
xd = [0.1;0.0;0.2];
dxd = [0;0;0];
dxr = [0;0;0];


conv_pcc = 1;
conv_motor = 1;
[tau, tau_r, x, M, C, A, cq] = helix_controller(q0,dq0,q0,q0,q0,d,m,r,kb,ks,bb,bs,bm,L0,Kp,KD,Kpx, KDx, xd, dxd, dxr, conv_pcc,conv_motor);
C_test = C_calc(q0,dq0,1/4*m*r^2, 1/4*m*r^2, 1/2*m*r^2,L0,r,m);
% Cdq = C_test(1:4,1:4)*dq0