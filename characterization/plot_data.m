%% Axial
clear
close all
segments = ["a", "b", "4", "5", "6", "7"];
% segments = ["1", "2", "3", "4", "5", "6"];

% segment = "a"
axial_array = zeros(6,5);
for s = 1:6
    segment = segments(s);
    for i = 1:5
        fname = sprintf('sept_24/mod_' + segment + '_compression/mod_' + segment + '_compression_test_trial_%i.is_comp_RawData/Specimen_RawData_1.csv', i);
        data_array = readtable(fname);
        N = size(data_array, 1);
        time = data_array.Var1;
        extension = data_array.Var2;
        load = data_array.Var3;
        x = extension(7:N); 
        y = load(7:N);
        slope = (y(2:end) - y(1:end-1)) ./ (x(2:end) - x(1:end-1));
        c = polyfit(x(4000:end-7000),y(4000:end-7000),1);
        y_est = polyval(c,x(4000:end-7000));
        % figure
        % plot(x, y);
        % hold on 
        % plot(x(4000:end-7000),y_est,'r--','LineWidth',2)
        slope = (y_est(end) - y_est(1))/(x(end-7000) - x(4000))*1000;
        axial_array(s,i) = slope;
        % xlabel('extension (mm)'); 
        % ylabel('load (N)'); 
        % legend({'mod 5', 'linear fit'})
        % 
    end
end

[Y,M] = std(axial_array,0,2)


%% Bending

close all
segments = ["a", "b", "4", "5", "6", "7"];
% segment = "a"
bending_array = zeros(6,5);
for s = 1:6
    segment = segments(s);
    for i = 1:5
        fname = 'bend_test_sept_26/mod_' + segment + '_bend_test/mod_' + segment + '_bend_test_trial_' + i + '.is_tens_RawData/Specimen_RawData_1.csv';
        data_array = readtable(fname);
        N = size(data_array, 1);
        time = data_array.Var1;
        extension = data_array.Var2;
        load = data_array.Var3;
        x = extension(7:N); 
        y = load(7:N);
        figure
        plot(x, y, 'LineWidth', 2); 
        hold on 
        % print smoothed out signal 
        windowSize = 5;
        y_filtered = movmean(y, windowSize);
        plot(x, y_filtered, 'r-', 'LineWidth', 2, 'DisplayName', 'Filtered Data');
        
        
        if s == 6
           % define the range (delta x) we will look at 
            a = 130;     % starting index
            b = 140;    % ending index
        else
            a = 100;
            b = 110;
        end

     
        % range of extension we'll look at
        c = polyfit(x(a:b),y_filtered(a:b),1);
        y_est = polyval(c,x(a:b));
        plot(x(a:b), y_filtered(a:b), 'g-', 'LineWidth',2);
    
        % calculate bending stiffness
        L0 = 60; % in mm 
        l1 = L0;
        l2 = L0;
        l3 = L0 - x;
        dl1 = 0;
        dl2 = 0;
        dl3 = x(b)-x(a);        % in mm 
        d = 30/1000; % based off of CAD radius measurement (in mm)
        % calculating curvature
        k = 2 * sqrt(l1^2 + l2^2 + l3.^2 - l1*l2 - l2*l3 - l1*l3)./(d*(l1 + l2+ l3));
        delta_k = 2 * sqrt(dl1^2 + dl2^2 + dl3.^2 - dl1*dl2 - dl2*dl3 - dl1*dl3)./(d*(dl1 + dl2+ dl3));
        dy = y_filtered(b) - y_filtered(a);
        bending_stiffness = dy * d^2 / (dl3/1000)
        bending_array(s,i) = bending_stiffness;
        % dx = k(b) - k(a);
        % stiffness = (y_filtered(b)-y_filtered(a))/(k(b) - k(a))
        title_name = "bending stiffness: " + bending_stiffness;
        xlabel('extension (mm)'); 
        ylabel('load (N)'); 
        legend({'segment b', 'smooth segment b', 'for calculating stiffness'})
        title(title_name);
        % Specify folder path and file name
        % folder = '/home/ranger/Desktop/robosoft_2024_characterization/esolo/segment_' + segment;  % Replace with your desired folder
        % filename = 'segment_' + segment + '_trial_' + i + '.png';  % File name with extension
        % 
        % % Create folder if it doesn't exist
        % if ~exist(folder, 'dir')
        %    mkdir(folder);
        % end
        % 
        % % Save the figure to the specified folder
        % saveas(gcf, fullfile(folder, filename));
    end
end
close all
[Y,M] = std(bending_array,0,2)
errorbar(M,Y)
% plot(slope, x(2:end))

%% axial compression

ShA = 45;
E = 0.0981 * (56 + 7.62336 * ShA)/(0.137505 * (254 - 2.54 * ShA))*1e6;
w = 8/1000;
t = 4/1000;
H = 0.12;
alpha = 32.48/2;
D = 60/1000;
Nh = 3;
Np = 2;
alpha = atan(2*H/(pi * (Np)*D))*180/pi;
h = pi/2 * D / Nh * tand(alpha) - t;
di = (pi * (D - w) - Nh*2*t)/Nh;
I = 1/12 * w * t^3;
% alpha_i = D * tand(alpha) / 2 ./ (D - 2*w) * 180/pi;
alpha_i = atan(2*H./(pi * (Np)*(D - w - 2*t)))* 180/pi;
Li = sqrt(h^2 + di.^2) / 2;
di_avg = 1/w*(- 2*t*w - (w*(pi*w - pi*D))/Nh);
Li_avg = sqrt(h^2 + di_avg.^2) / 2;
wv = w;
alpha_i_avg = 180/pi*atan(1/w*(H*(log(D) - log(D - 2*wv)))/(Np*pi));
alpha_i_avg = atan(2*H./(pi * (Np)*(D - w)))* 180/pi;
F = 0.5;
y = 1/12 * F * (Li_avg) .^3 ./ I / E .* cosd(alpha_i_avg).^2 % this model works for compression
y = 1/12 * F * (Li) .^3 ./ I / E .* cosd(alpha_i).^2 % this model works for compression

k = F/(y)
% wv = linspace(0.0000001,w);
% y = 0;
% for i=1:length(wv)
%     w= wv(i);
%     alpha = atan(2*H/(pi * (Np)*D))*180/pi;
%     h = pi/2 * D / Nh * tand(alpha);
%     di = (pi * (D - 2*w) - Nh*2*t)/Nh;
%     % alpha_i = D * tand(alpha) / 2 ./ (D - 2*w) * 180/pi;
%     alpha_i = atan(2*H./(pi * (Np)*(D - 2*w)))* 180/pi;
%     Li = sqrt(h^2 + di.^2) / 2;
% 
%     F = 0.5;
%     y = y + 1/12 * F * (Li) .^3 ./ I / E .* cosd(alpha_i).^2;
% 
% end
% yavg = y/100
% y = (F * (D - w).^3./I / E .* cosd(alpha_i).^2*(pi/4 - 1/pi)/2);
% y = (1/12 *F/E* (Li) .^3 ./ I.* cosd(alpha_i).^2 + 1/12 *F/E* sqrt(h^2/4 + D^2/4) .^3 ./ I .* cosd(alpha))/2;


%% bending
M = 0.001;
r = D/2;
ri = r - w;
If = pi/4 * (r^4 -ri^4);
F = (r^3*4)/3 / If / pi% integral average of stress over cylinder
% F = 1 / r
% F = - 8*M/pi * (1/r - 1/ri)
% F = 8/3/pi * (r^3 - ri^3)/(r^4 - ri^4)
Nh = 3;
% F = r^3 /If * 4/3*sin(pi/Nh)

alpha = 32.48/2;
h = pi/2 * D * tand(alpha*2);
% F = 4 * M / pi / (D/2) * (1 - cos(pi/Nh))

y = 1/12 * F * (1/Nh*sqrt(h^2/4 + (pi*D - 2*t*Nh)^2/4)) .^3 ./ I / E * cosd(alpha)^2  %needs * 4/3, 2 from other side
% th_e = 1/8 * F * cos(alpha) * sqrt(h^2/4 + (pi*D/Nh)^2/4) .^2 ./ I / E / (Nh/pi);
L0 = 0.06;
% y = 1/(1/(Nh*y) + 2/(Nh*y/2))
% y = Nh*y/4
% y = 1/(1/(Nh*y) + 2/(Nh*y/2))/2
L = L0 - y/1000;
p = r/(L - L0) * L0;
% th = L0/p*180/pi
% k = M/th*180/pi

k = 1/(y/r)
k = k/Nh*Nh/pi*2