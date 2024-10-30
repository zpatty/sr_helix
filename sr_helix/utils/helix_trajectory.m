%% 3D Helix Arm Trajectory
d = 0.04337;
theta = [0 10 10 0;
        0 10 10 0.0;
        0 10 10 0];
phi = [0 0 pi 0;
        0 0 pi 0;
        0 0 pi 0];
dL = [-0.0 -0.04 -0.04 -0.0;
    -0.0 -0.04 -0.04 -0.0;
    -0.0 -0.04 -0.04 -0.0];

dx = theta.*d.*cos(phi);
dy = theta.*d.*sin(phi);

motor = [0 0 0 0];

wayPoints = [motor(1,:);dx(1,:);dy(1,:);dL(1,:);dx(2,:);dy(2,:);dL(2,:);
    dx(3,:);dy(3,:);dL(3,:)];
timePoints = [0 3 9 12];
tvec = 0:0.05:12;
N = length(tvec);
[qd,dqd,ddqd,pp] = cubicpolytraj(wayPoints,timePoints,tvec);

save("~/DRL/sr_helix/sr_helix/conf/qd.mat", "qd")
save("~/DRL/sr_helix/sr_helix/conf/dqd.mat", "dqd")
save("~/DRL/sr_helix/sr_helix/conf/ddqd.mat", "ddqd")
save("~/DRL/sr_helix/sr_helix/conf/tvec.mat", "tvec")



