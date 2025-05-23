clear
close all

% analytical
Nhv = [3, 4, 5];
% variables = [Nh; w; t; D; H];
for j = 1:3
    ShA = 45;
    E = 0.0981 * (56 + 7.62336 * ShA)/(0.137505 * (254 - 2.54 * ShA))*1e6;
    w = 8/1000;
    t = 4/1000;
    H = 0.145;
    D = 60/1000;
    Nh = Nhv(j);
    Np = 2;
    alpha = atan(2*H/(pi * (Np)*D))*180/pi;
    alpha_mid = atan(H/2 / (pi*(D - w)))*180/pi;
    h = pi/2 * D / Nh * tand(alpha) ;
    di = (pi * (D - w) - Nh*2*t)/Nh;
    I = 1/12 * w * t^3;
    % alpha_i = D * tand(alpha) / 2 ./ (D - 2*w) * 180/pi;
    alpha_i = atan(2*H./(pi * (Np)*(D - w - 2*t)))* 180/pi;
    si = 2*pi/Nh/2 * (D/2 - w/2) - t/2 / sind(alpha_mid);
    hi = (h - t) / 2;
    Li = sqrt(h^2 + di.^2) / 2;
    L = (1/Nh*sqrt(h^2/4 + (pi*D - 2*t*Nh)^2/4));
    L = sqrt((h/2)^2 + (pi*D/Nh/2)^2);
    Li = sqrt(hi^2 + si^2);
    % L = 0.0203
    % L = 0.0226
    di_avg = 1/w*(- 2*t*w - (w*(pi*w - pi*D))/Nh);
    Li_avg = sqrt(h^2 + di_avg.^2) / 2;
    wv = w;
    alpha_i_avg = 180/pi*atan(1/w*(H*(log(D) - log(D - 2*wv)))/(Np*pi));
    alpha_i_avg = atan(2*H./(pi * (Np)*(D - w)))* 180/pi;
    F = 0.1;
    y = 1/12 * F * (Li_avg) .^3 ./ I / E .* cosd(alpha_i_avg/2).^2; % this model works for compression
    % y = 1/12 * F * (L) .^3 ./ I / E * cosd(alpha/2)^2  % this model works for compression
    % y = 1/12 * F * (Li_avg) .^3 ./ I / E
    y = 1/12 * F * (Li) .^3 ./ I / E .* cosd(alpha_mid).^2;  % this model works for compression
    % y = 1/12 * F * (Li) .^3 ./ I / E;
    
    knh(j) = F/(y);

    M = 0.001;
    r = D/2;
    ri = r - w;
    If = pi/4 * (r^4 -ri^4);
    F = (r^3*4)/3 / If / pi;% integral average of stress over cylinder
    F = 10*t*(r^2 - ri^2)/If;
    % F = 4 * M / pi / (D/2) * (1 - cos(pi/Nh))
    F = - 4 * (1/r - 1/ri);
    F = (r^2)/If*w*3;
    alpha = atan(H/2 / (pi*(D)))*180/pi;
    so = 2*pi/Nh/2 * (D/2) - t/2 / sind(alpha);
    % F = 4 * M / pi / (D/2) * (1 - cos(pi/Nh))
    Lo = sqrt(hi^2 + so^2);
    y = 1/12 * F * (Lo) .^3 ./ I / E * cosd(alpha)^2/2;  %needs * 4/3, 2 from other side
    
    L0 = H/2;
    L = L0 - y/1000;
    p = r/(L - L0) * L0;
    r = D/2;
    ri = r - w;
    rm = (r+ri)/2;
    If = pi/4 * (rm^4);
    A = pi*(rm)^2;
    kbNh(j) = 1/(y/r) /Nh*Nh/pi;
    kbNh(j) = knh(j)*If/A*rm/H*3*3;
end

wve = linspace(2,20);
% wve = [4, 8, 12];
for j = 1:length(wve)
    ShA = 45;
    E = 0.0981 * (56 + 7.62336 * ShA)/(0.137505 * (254 - 2.54 * ShA))*1e6;
    w = wve(j)/1000;
    t = 4/1000;
    H = 0.145;
    D = 60/1000;
    Nh = 4;
    Np = 2;
    alpha = atan(2*H/(pi * (Np)*D))*180/pi;
    alpha_mid = atan(H/2 / (pi*(D - w)))*180/pi;
    h = pi/2 * D / Nh * tand(alpha) ;
    di = (pi * (D - w) - Nh*2*t)/Nh;
    I = 1/12 * w * t^3;
    % alpha_i = D * tand(alpha) / 2 ./ (D - 2*w) * 180/pi;
    alpha_i = atan(2*H./(pi * (Np)*(D - w - 2*t)))* 180/pi;
    si = 2*pi/Nh/2 * (D/2 - w/2) - t/2 / sind(alpha_mid);
    hi = (h - t) / 2;
    Li = sqrt(h^2 + di.^2) / 2;
    L = (1/Nh*sqrt(h^2/4 + (pi*D - 2*t*Nh)^2/4));
    L = sqrt((h/2)^2 + (pi*D/Nh/2)^2);
    Li = sqrt(hi^2 + si^2);
    % L = 0.0203
    % L = 0.0226
    di_avg = 1/w*(- 2*t*w - (w*(pi*w - pi*D))/Nh);
    Li_avg = sqrt(h^2 + di_avg.^2) / 2;
    wv = w;
    alpha_i_avg = 180/pi*atan(1/w*(H*(log(D) - log(D - 2*wv)))/(Np*pi));
    alpha_i_avg = atan(2*H./(pi * (Np)*(D - w)))* 180/pi;
    F = 0.1;
    y = 1/12 * F * (Li_avg) .^3 ./ I / E .* cosd(alpha_i_avg/2).^2; % this model works for compression
    % y = 1/12 * F * (L) .^3 ./ I / E * cosd(alpha/2)^2  % this model works for compression
    % y = 1/12 * F * (Li_avg) .^3 ./ I / E
    y = 1/12 * F * (Li) .^3 ./ I / E .* cosd(alpha_mid).^2;  % this model works for compression
    % y = 1/12 * F * (Li) .^3 ./ I / E;
    
    kw(j) = F/(y);

    M = 0.001;
    r = D/2;
    ri = r - w;
    If = pi/4 * (r^4 -ri^4);
    F = ((r^3 - ri^3))/3 / If * pi;% integral average of stress over cylinder
    F = 10*t*(r^2 - ri^2)/If;
    F = 3/w*(r^2 - ri^2)/If;
    rm = (r+ri)/2;
    If = pi/4 * (rm^4);
    F = (r^3*4)/3 / If / pi;% integral average of stress over cylinder
    ri = r - w;
    % F = 1 / r^(3);
    % F = 4*w/r^2;
    % F = 4*w/r^2*0.65;
    % F = - 4 * (1/r - 1/ri);
    % F = ((r^3 - ri^3)*2)/3 / If * pi;
    % F = 8/3/pi * (r^3 - ri^3)/(r^4 - ri^4);
    F = (r^2 - ri^2)/If*w*3;
    F = (rm^2)/If*w*3/2;
    % F = 4 * (1/r);
    % F = 4 * M / pi / (D/2) * (1 - cos(pi/Nh))
    alpha = atan(H/2 / (pi*(D)))*180/pi;
    so = 2*pi/Nh/2 * (D/2) - t/2 / sind(alpha);
    % F = 4 * M / pi / (D/2) * (1 - cos(pi/Nh))
    Lo = sqrt(hi^2 + so^2);
    y = 1/12 * F * (Lo) .^3 ./ I / E * cosd(alpha)^2/2;  %needs * 4/3, 2 from other side
    
    L0 = H/2;
    L = L0 - y/1000;
    p = r/(L - L0) * L0;
    
    % A = pi*(rm)^2;
    % A = pi*(r^2 - ri^2);
    If = pi/4 * (rm^4);
    A = pi*(rm)^2;
    kbw(j) = 1/(y/r) /Nh*Nh/pi;
    kbw(j) = kw(j)*If/A*rm/H*3*3;
end

tv = linspace(1,8);
% tv = [2,4,6];
for j = 1:length(wve)
    ShA = 45;
    E = 0.0981 * (56 + 7.62336 * ShA)/(0.137505 * (254 - 2.54 * ShA))*1e6;
    w = 8/1000;
    t = tv(j)/1000;
    H = 0.145;
    D = 60/1000;
    Nh = 4;
    Np = 2;
    alpha = atan(2*H/(pi * (Np)*D))*180/pi;
    alpha_mid = atan(H/2 / (pi*(D - w)))*180/pi;
    h = pi/2 * D / Nh * tand(alpha) ;
    di = (pi * (D - w) - Nh*2*t)/Nh;
    I = 1/12 * w * t^3;
    % alpha_i = D * tand(alpha) / 2 ./ (D - 2*w) * 180/pi;
    alpha_i = atan(2*H./(pi * (Np)*(D - w - 2*t)))* 180/pi;
    si = 2*pi/Nh/2 * (D/2 - w/2) - t/2 / sind(alpha_mid);
    hi = (h - t) / 2;
    Li = sqrt(h^2 + di.^2) / 2;
    L = (1/Nh*sqrt(h^2/4 + (pi*D - 2*t*Nh)^2/4));
    L = sqrt((h/2)^2 + (pi*D/Nh/2)^2);
    Li = sqrt(hi^2 + si^2);
    % L = 0.0203
    % L = 0.0226
    di_avg = 1/w*(- 2*t*w - (w*(pi*w - pi*D))/Nh);
    Li_avg = sqrt(h^2 + di_avg.^2) / 2;
    wv = w;
    alpha_i_avg = 180/pi*atan(1/w*(H*(log(D) - log(D - 2*wv)))/(Np*pi));
    alpha_i_avg = atan(2*H./(pi * (Np)*(D - w)))* 180/pi;
    F = 0.1;
    y = 1/12 * F * (Li_avg) .^3 ./ I / E .* cosd(alpha_i_avg/2).^2; % this model works for compression
    % y = 1/12 * F * (L) .^3 ./ I / E * cosd(alpha/2)^2  % this model works for compression
    % y = 1/12 * F * (Li_avg) .^3 ./ I / E
    y = 1/12 * F * (Li) .^3 ./ I / E .* cosd(alpha_mid).^2;  % this model works for compression
    % y = 1/12 * F * (Li) .^3 ./ I / E;
    
    kt(j) = F/(y);

    M = 0.001;
    r = D/2;
    ri = r - w;
    If = pi/4 * (r^4 -ri^4);
    F = (r^3*4)/3 / If / pi;% integral average of stress over cylinder
    ri = r - w;
    % F = 1 / r^(3);
    % F = 4/w*log(r/(r-w));
    F = - 4 * (1/r - 1/ri);
    F = 10*t*(r^2 - ri^2)/If;
    F = (r^2)/If*w*3;
    % F = - 4 * (1/r - 1/ri);
% F = 8/3/pi * (r^3 - ri^3)/(r^4 - ri^4)
    alpha = atan(H/2 / (pi*(D)))*180/pi;
    so = 2*pi/Nh/2 * (D/2) - t/2 / sind(alpha);
    % F = 4 * M / pi / (D/2) * (1 - cos(pi/Nh))
    Lo = sqrt(hi^2 + so^2);
    y = 1/12 * F * (Lo) .^3 ./ I / E * cosd(alpha)^2/2;  %needs * 4/3, 2 from other side
    r = D/2;
    ri = r - w;
    rm = (r+ri)/2;
    If = pi/4 * (rm^4);
    L0 = H/2;
    L = L0 - y/1000;
    p = r/(L - L0) * L0;

    A = pi*(rm)^2;
    kbt(j) = 1/(y/r) /Nh*Nh/pi;
    kbt(j) = kt(j)*If/A*rm/H*3*3;
end

Dv = linspace(20,100);
Dv = [40,60,80];
for j = 1:length(Dv)
    ShA = 45;
    E = 0.0981 * (56 + 7.62336 * ShA)/(0.137505 * (254 - 2.54 * ShA))*1e6;
    w = 8/1000;
    t = 4/1000;
    H = 0.145;
    D = Dv(j)/1000;
    Nh = 4;
    Np = 2;
    alpha = atan(2*H/(pi * (Np)*D))*180/pi;
    alpha_mid = atan(H/2 / (pi*(D - w)))*180/pi;
    h = pi/2 * D / Nh * tand(alpha) ;
    di = (pi * (D - w) - Nh*2*t)/Nh;
    I = 1/12 * w * t^3;
    % alpha_i = D * tand(alpha) / 2 ./ (D - 2*w) * 180/pi;
    alpha_i = atan(2*H./(pi * (Np)*(D - w - 2*t)))* 180/pi;
    si = 2*pi/Nh/2 * (D/2 - w/2) - t/2 / sind(alpha_mid);
    hi = (h - t) / 2;
    Li = sqrt(h^2 + di.^2) / 2;
    L = (1/Nh*sqrt(h^2/4 + (pi*D - 2*t*Nh)^2/4));
    L = sqrt((h/2)^2 + (pi*D/Nh/2)^2);
    Li = sqrt(hi^2 + si^2);
    di_avg = 1/w*(- 2*t*w - (w*(pi*w - pi*D))/Nh);
    Li_avg = sqrt(h^2 + di_avg.^2) / 2;
    wv = w;
    alpha_i_avg = 180/pi*atan(1/w*(H*(log(D) - log(D - 2*wv)))/(Np*pi));
    alpha_i_avg = atan(2*H./(pi * (Np)*(D - w)))* 180/pi;
    F = 0.1;
    y = 1/12 * F * (Li_avg) .^3 ./ I / E .* cosd(alpha_i_avg/2).^2; % this model works for compression
    % y = 1/12 * F * (L) .^3 ./ I / E * cosd(alpha/2)^2  % this model works for compression
    % y = 1/12 * F * (Li_avg) .^3 ./ I / E
    y = 1/12 * F * (Li) .^3 ./ I / E .* cosd(alpha_mid).^2;  % this model works for compression
    Di = D - 2*w;
    r = D/2;
    ri = r - w;
    alpha = atan(H/2 / (pi*(D)))*180/pi;
    y = 1/12 * F * (Li) .^3 ./ I / E .* cosd(alpha_mid).^2;
    kD(j) = F/(y);

    
    Dm = (D + Di)/2;
    y = Dm^3/E/w/t^3/(Nh*2)^4*Nh*Di/D*3.88;
    
    mu = w/Dm;
    G = E/2/(1.5);
    S = w*t^3/16*(16/3 - 3.36*t/w*(1-t^4/12/w^4));
    lam = E*1/12*w*t^3/G/S;

    % kD(j) = E*t^3/Dm^2 * 4*mu*Nh * ((1+lam)*pi/Nh + (lam-1)*sin(pi/Nh))/(3*lam*(pi/Nh*((1+lam)*pi/Nh + (lam-1)*sin(pi/Nh)) - 4*lam*(1-cos(pi/Nh))));

    M = 0.001;
    r = D/2;
    ri = r - w;
    rm = (r+ri)/2;
    If = pi/4 * (r^4);
    
    F = (r^3*4)/3 / If / pi;% integral average of stress over cylinder
    ri = r - w;
    % F = 1 / r^(3);
    % F = 4*w/r^2;
    % F = 4*w/r^2*0.65;
    % F = - 4 * (1/r - 1/ri);
    % F = ((r^3 - ri^3)*2)/3 / If * pi;
    % F = 8/3/pi * (r^3 - ri^3)/(r^4 - ri^4);
    F = (r^2 - ri^2)/If*w*3;
    F = (rm^2)/If*w*3/2;
    alpha = atan(H/2 / (pi*(D)))*180/pi;
    th = t/2/sind(alpha)/r;
    Iann(j) = (r^4 - ri^4)*((th + sin(th))/8 - 8*sin(th/2)^2/9/th) + (8*sin(th/2)^2/9/th/(r-ri)^2 * (ri^4*r^2-r^4*ri^2));
    % F = -pi*(r^3-ri^3)/3/Iann/5;
    alpha = atan(H/2 / (pi*(D)))*180/pi;
    so = 2*pi/Nh/2 * (D/2) - t/2 / sind(alpha);
    % F = 4 * M / pi / (D/2) * (1 - cos(pi/Nh))
    Lo = sqrt(hi^2 + so^2);
    y = 1/12 * F * (Lo) .^3 ./ I / E * cosd(alpha)^2;  %needs * 4/3, 2 from other side
    
    L0 = H/2;
    L = L0 - y/1000;
    p = r/(L - L0) * L0;
    
    theta = (L0 - y)./(r*L0/y - r);
    
    kbD(j) = 1/(y/r) /Nh*Nh/pi;
    % kbD(j) = 2/pi*kD(j)*Dm^2*4;
    % A = pi*(r^2 - ri^2);
    If = pi/4 * (rm^4);
    A = pi*(rm)^2;
    y = 1/12 * (Lo) .^3 ./ I / E * cosd(alpha)^2;  %needs * 4/3, 2 from other side
    kbD(j) = kD(j)*If/A*r/ri/2;
    kbD(j) = kD(j)*If/A*rm/H*3*3;
end

Hv = linspace(0.05,0.2);
Hv = [0.12, 0.145, 0.17];
for j = 1:length(Hv)
    ShA = 45;
    E = 0.0981 * (56 + 7.62336 * ShA)/(0.137505 * (254 - 2.54 * ShA))*1e6;
    w = 8/1000;
    t = 4/1000;
    H = Hv(j);
    D = 60/1000;
    Nh = 4;
    Np = 2;
    alpha = atan(2*H/(pi * (Np)*D))*180/pi;
    alpha_mid = atan(H/2 / (pi*(D - w)))*180/pi;
    h = pi/2 * D / Nh * tand(alpha) ;
    di = (pi * (D - w) - Nh*2*t)/Nh;
    I = 1/12 * w * t^3;
    % alpha_i = D * tand(alpha) / 2 ./ (D - 2*w) * 180/pi;
    alpha_i = atan(2*H./(pi * (Np)*(D - w - 2*t)))* 180/pi;
    si = 2*pi/Nh/2 * (D/2 - w/2) - t/2 / sind(alpha_mid);
    hi = (h - t) / 2;
    Li = sqrt(h^2 + di.^2) / 2;
    L = (1/Nh*sqrt(h^2/4 + (pi*D - 2*t*Nh)^2/4));
    L = sqrt((h/2)^2 + (pi*D/Nh/2)^2);
    Li = sqrt(hi^2 + si^2);
    % L = 0.0203
    % L = 0.0226
    di_avg = 1/w*(- 2*t*w - (w*(pi*w - pi*D))/Nh);
    Li_avg = sqrt(h^2 + di_avg.^2) / 2;
    wv = w;
    alpha_i_avg = 180/pi*atan(1/w*(H*(log(D) - log(D - 2*wv)))/(Np*pi));
    alpha_i_avg = atan(2*H./(pi * (Np)*(D - w)))* 180/pi;
    F = 1;
    y = 1/12 * F * (Li_avg) .^3 ./ I / E .* cosd(alpha_i_avg/2).^2; % this model works for compression
    % y = 1/12 * F * (L) .^3 ./ I / E * cosd(alpha/2)^2  % this model works for compression
    % y = 1/12 * F * (Li_avg) .^3 ./ I / E
    y = 1/12 * F * (Li) .^3 ./ I / E .* cosd(alpha_mid).^2;  % this model works for compression
    % y = 1/12 * F * (Li) .^3 ./ I / E ;
    kH(j) = 1/y;
    theta = 0;
    G = E/2/(1.5);
    I = 1/12*w*t^3;
    R = D/2;
    K = w*t^3/16*(16/3 - 3.36*t/w*(1-t^4/12/w^4));
    beta = E*I/G/K;
    alpha = atan(H/2 / (pi*(D)))*180/pi;
    alpha_mid = atan(H/2 / (pi*(D - w)))*180/pi;
    phi = H/R*tand(alpha_mid)/Nh;
    phi = Li/R;
    
    
    Ca1 = (1 + beta)/2 * ( (phi - theta)*sin(phi - theta) - beta*(1 - cos(phi - theta)) );
    Ca2 = (1 + beta)/2 * ( (phi - theta)*cos(phi - theta) - sin(phi - theta) );
    Ca3 = -beta*( (phi - theta) - sin(phi - theta) ) - Ca2;
    Ca4 = (1 + beta)/2 * ( (phi - theta)*cos(phi - theta) ) + (1 - beta)/2 * sin(phi - theta);
    Ca5 = - (1 + beta)/2 * (phi - theta) * sin(phi - theta);
    Ca6 = Ca1;
    Ca7 = Ca5;
    Ca8 = (1 - beta)/2 * sin(phi - theta) - (1 + beta)/2 * ( (phi - theta)*cos(phi - theta) );
    Ca9 = Ca2;
    
    
    
    
    F = 1;
    M = -(F*R*(Ca5*Ca9 - Ca6*Ca8))/(Ca4*Ca8 - Ca5*Ca7);
    T = (F*R*(Ca4*Ca9 - Ca6*Ca7))/(Ca4*Ca8 - Ca5*Ca7);
    F = 1*cosd(alpha_mid);
    s = sin(phi);
    c = cos(phi);
    y = T*R^2/E/I * (Ca5*s - Ca8*(1-c) - Ca2) + M*R^2/E/I * (Ca4*s - Ca7*(1-c) - Ca1) - F*R^3/E/I * (Ca6*s - Ca9*(1-c) - Ca3);

    % kH(j) = 1/(-y);

    M = 0.001;
    r = D/2;
    ri = r - w;
    If = pi/4 * (r^4 -ri^4);
    F = (r^3*4)/3 / If / pi;% integral average of stress over cylinder
    F = 10*t*(r^2 - ri^2)/If;
    F = (r^2)/If*w*3;
    so = 2*pi/Nh/2 * (D/2) - t/2 / sind(alpha);
    % F = 4 * M / pi / (D/2) * (1 - cos(pi/Nh))
    Lo = sqrt(hi^2 + so^2);
    y = 1/12 * F * (Lo) .^3 ./ I / E * cosd(alpha)^2;  %needs * 4/3, 2 from other side
    
    L0 = H/2;
    L = L0 - y/1000;
    p = r/(L - L0) * L0;
    r = D/2;
    ri = r - w;
    rm = (r+ri)/2;
    If = pi/4 * (rm^4);
    A = pi*(rm)^2;
    theta = (L0 - y)./(r*L0/y - r);
    kbH(j) = 1/(y/r) /Nh*Nh/pi;
    kbH(j) = kH(j)*If/A*rm/H*3*3;
end

fkNh = 0.1./[1.2, 0.44, 0.1887]*1000;
fkw = 0.1./[1.224, 0.44, 0.2062]*1000;
fkt = 0.1./[1.8, 0.44, 0.1]*1000;
fkD = 0.1./[0.2076, 0.44, 0.8277]*1000;
fkH = 0.1./[0.2884, 0.44, 0.5878]*1000;

fkbNh = 0.001./[2.19, 0.79, 0.34]/(pi/180);
fkbw = 0.001./[1.8, 0.79, 0.48]/(pi/180);
fkbt = 0.001./[6.55, 0.79, 0.19]/(pi/180);
fkbD = 0.001./[1.0, 0.79, 0.78]/(pi/180);
fkbH = 0.001./[0.53, 0.79, 1.05]/(pi/180);


savef = 0;
% Compression
figure
yyaxis left
plot(Nhv,knh, "LineWidth",2.0)
hold on
scatter(Nhv,fkNh, "LineWidth",2.0)
ylabel("Axial Stiffness (N/m)", 'Interpreter', 'latex', 'FontSize', 28)

yyaxis right
plot(Nhv,kbNh, "LineWidth",2.0)
hold on
scatter(Nhv,fkbNh, "LineWidth",2.0)
ylabel("Bending Stiffness (Nm/rad)", 'Interpreter', 'latex', 'FontSize', 28)

legend(["Analytical", "FEM"])
xlabel("$N_h$", 'Interpreter', 'latex', 'FontSize', 28)
set(gca,'FontSize',28)
ax = gca;
set(ax, 'TickLabelInterpreter', 'latex')
set(gcf, 'Position', [303 495 1120 630])
% Modify the figure size
% set(gcf, 'Position', [303 495 1120 630])
% Get the legend and set its interpreter
l = findobj(gcf, 'Type', 'legend');
set(l, 'Interpreter', 'latex')
% Turn on grid and box
grid on
box on 
fig = gcf;
if savef == 1
    print(fig, "Nh", '-dpdf')
end

%%%% w %%%%
figure
yyaxis left
plot(wve,kw,'LineWidth',2)
hold on
scatter([4, 8, 12],fkw, "LineWidth",2.0)
ylabel("Axial Stiffness (N/m)", 'Interpreter', 'latex', 'FontSize', 28)

yyaxis right
plot(wve,kbw,'LineWidth',2)
hold on
scatter([4, 8, 12],fkbw, "LineWidth",2.0)
ylabel("Bending Stiffness (Nm/rad)", 'Interpreter', 'latex', 'FontSize', 28)

legend(["Analytical", "FEM"])
xlabel("$w$ (m)", 'Interpreter', 'latex', 'FontSize', 28)
set(gca,'FontSize',28)
ax = gca;
set(ax, 'TickLabelInterpreter', 'latex')
set(gcf, 'Position', [303 495 1120 630])
% Modify the figure size
% set(gcf, 'Position', [303 495 1120 630])
% Get the legend and set its interpreter
l = findobj(gcf, 'Type', 'legend');
set(l, 'Interpreter', 'latex')
% Turn on grid and box
grid on
box on 
fig = gcf;
if savef == 1
    print(fig, "w", '-dpdf')
end

%%%% t %%%%
figure
yyaxis left
plot(tv,kt,'LineWidth',2)
hold on
scatter([2, 4, 6],fkt, "LineWidth",2.0)
ylabel("Axial Stiffness (N/m)", 'Interpreter', 'latex', 'FontSize', 28)

yyaxis right
plot(tv,kbt,'LineWidth',2)
hold on
scatter([2, 4, 6],fkbt, "LineWidth",2.0)
ylabel("Bending Stiffness (Nm/rad)", 'Interpreter', 'latex', 'FontSize', 28)

legend(["Analytical", "FEM"])
xlabel("$t$ (m)", 'Interpreter', 'latex', 'FontSize', 28)
set(gca,'FontSize',28)
ax = gca;
set(ax, 'TickLabelInterpreter', 'latex')
set(gcf, 'Position', [303 495 1120 630])
% Modify the figure size
% set(gcf, 'Position', [303 495 1120 630])
% Get the legend and set its interpreter
l = findobj(gcf, 'Type', 'legend');
set(l, 'Interpreter', 'latex')
% Turn on grid and box
grid on
box on 
fig = gcf;
if savef == 1
    print(fig, "t", '-dpdf')
end

%%%% D %%%%
figure
yyaxis left
plot(Dv,kD,'LineWidth',2)
hold on
scatter([40,60,80],fkD, "LineWidth",2.0)
ylabel("Axial Stiffness (N/m)", 'Interpreter', 'latex', 'FontSize', 28)

yyaxis right
plot(Dv,kbD,'LineWidth',2)
hold on
scatter([40,60,80],fkbD, "LineWidth",2.0)
ylabel("Bending Stiffness (Nm/rad)", 'Interpreter', 'latex', 'FontSize', 28)

legend(["Analytical", "FEM"])
xlabel("$D$ (m)", 'Interpreter', 'latex', 'FontSize', 28)
set(gca,'FontSize',28)
ax = gca;
set(ax, 'TickLabelInterpreter', 'latex')
set(gcf, 'Position', [303 495 1120 630])
% Modify the figure size
% set(gcf, 'Position', [303 495 1120 630])
% Get the legend and set its interpreter
l = findobj(gcf, 'Type', 'legend');
set(l, 'Interpreter', 'latex')
% Turn on grid and box
grid on
box on 
fig = gcf;
if savef == 1
    print(fig, "D", '-dpdf')
end

%%%%% H %%%%%
figure
yyaxis left
plot(Hv,kH,'LineWidth',2)
hold on
scatter([0.12,0.145,0.17],fkH, "LineWidth",2.0)
ylabel("Axial Stiffness (N/m)", 'Interpreter', 'latex', 'FontSize', 28)

yyaxis right
plot(Hv,kbH,'LineWidth',2)
hold on
scatter([0.12,0.145,0.17],fkbH, "LineWidth",2.0)
ylabel("Bending Stiffness (Nm/rad)", 'Interpreter', 'latex', 'FontSize', 28)

legend(["Analytical", "FEM"])
xlabel("$H$ (m)", 'Interpreter', 'latex', 'FontSize', 28)
set(gca,'FontSize',28)
ax = gca;
set(ax, 'TickLabelInterpreter', 'latex')
set(gcf, 'Position', [303 495 1120 630])
% Modify the figure size
% set(gcf, 'Position', [303 495 1120 630])
% Get the legend and set its interpreter
l = findobj(gcf, 'Type', 'legend');
set(l, 'Interpreter', 'latex')
% Turn on grid and box
grid on
box on 
fig = gcf;
if savef == 1
    print(fig, "H", '-dpdf')
end