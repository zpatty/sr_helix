%% simulate puppet segment
clear
close all
d = 0.04;
theta = [0 0.4 0.1 -0.2;
        0 0.0 0.0 0.0;
        0 1 2 2;
        0 0 -0.2 -0.5];
phi = [0 pi/6 pi/4 pi/3;
        0 0 0 0;
        0 -pi/4 -pi/4 -pi/4;
        0 0 0 0];
dL = [-0.0 -0.02 -0.05 -0.05;
    -0.0 -0.065 -0.065 -0.065;
    -0.0 -0.0 -0.0 -0.0;
    -0.0 -0.01 -0.02 -0.02];

dx = theta.*d.*cos(phi);
dy = theta.*d.*sin(phi);

motor = [0 0 pi/4 pi/3;
        0 -0.0 -pi/4 -pi/4;
        0 pi/4 pi/3 0.0;
        0 0 -pi/6 -pi/6];

wayPoints = [motor(1,:);dx(1,:);dy(1,:);dL(1,:);motor(2,:);dx(2,:);dy(2,:);dL(2,:);motor(3,:);
    dx(3,:);dy(3,:);dL(3,:);motor(4,:);dx(4,:);dy(4,:);dL(4,:)];
timePoints = [0 3 6 9];
fr = 0.05;
tSamples = 0:fr:9;
N = length(tSamples);
[qd,dqd,ddqd,pp] = cubicpolytraj(wayPoints,timePoints,tSamples);

L0 = 0.19;

q0 = 0.0*ones(10,1);%zeros(11,1);
% q0(end-1) = 0.02;
% qd = q0;
dq0 = zeros(10,1);

N = 4;
r = 0.01075;
d = 0.0395;
kb = 0;
ks = 0;
bb = 0;
bs = 0;
e = L0/1000;
segment_mass = 0.071;
med_plates = 0.036;
end_plate = 0.023;
m = 2*segment_mass + 2 * med_plates + end_plate;
kps = 0;
kpb = 0;
km = 0;
Kp = diag([km;kpb;kpb;kps;kpb;kpb;kps;kpb;kpb;kps]);
Dm = 0.0;
Ds = 0;
KD = diag([Dm;Ds;Ds;Ds;Ds;Ds;Ds;Ds;Ds;Ds]);
qhist = [q0];
tauhist = [];
chist = [];
Kpx = 50.0;
KDx = 50.0;
xd = [0.1;0.0;0.2];
dxd = [0;0;0];
dx = [0;0;0];
inds{1} = 2:4;
inds{2} = 6:8;
inds{3} = 10:12;
inds{4} = 14:16;
% qd = [0.01;0;0];
% dqd = [0;0;0];
% ddqd = [0;0;0];
k = 1;
offset = -0.005;
steepness = 100;
contact = 1;
attenuation = 1;
% q0 = [-0.491;0.03;-0.003;-0.008;1.108;0.01;-0.006;-0.019;0.245;-0.008;0.004;-0.001;-1.158;-0.003;0.035;-0.024];
% xd = [-0.2;0.1;0.2];
if size(qd(1,:)) == 1
    qd = repmat(qd,1,length(tSamples));
    dqd = zeros(16,length(tSamples));
    ddqd = zeros(16,length(tSamples));
end
% qd(:,1) = [0;0;0;-0.2;0;-0.13;0;-0.2;0;-0.13;0;-0.05;0;0;0.13;0];
conv_pcc = 1;
conv_motor = 1;
[tau, tau_r, x, M, C, A, cq] = helix_controller(q0,dq0,q0,q0,q0,N,d,m,r,N,kb,ks,bb,bs,L0,Kp,KD,Kpx, KDx, xd, dxd, dx, conv_pcc,conv_motor);
