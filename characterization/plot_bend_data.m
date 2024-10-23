
%% plotting time :))
clear
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


