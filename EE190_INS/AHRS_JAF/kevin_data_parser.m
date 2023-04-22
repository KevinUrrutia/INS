clear
lines_str = readlines('kevin_data.csv');
data_length = length(lines_str);
g = 9.8; % gravity

% Get first IMU measurement as reference of time
line0_str = lines_str(1);
line_0 = parse_line(line0_str);

% pre-allocation
elps_t = zeros(data_length,1);
ax = zeros(data_length,1);
ay = zeros(data_length,1);
az = zeros(data_length,1);
wx = zeros(data_length,1);
wy = zeros(data_length,1);
wz = zeros(data_length,1);

% parse measurements into arrays
for i = 1:length(lines_str)
   line_i = parse_line(lines_str(i));
   elps_t(i) = (line_i.hrs - line_0.hrs) * 3600 ...
                  + (line_i.min - line_0.min) * 60  ...
                  + (line_i.sec - line_0.sec);
   ax(i) = line_i.ax; % Confirm unit with Kevin, mG or m/s/s
   ay(i) = line_i.ay;
   az(i) = line_i.az;
   wx(i) = line_i.wx; % deg/s ? 
   wy(i) = line_i.wy;
   wz(i) = line_i.wz;
end


%% Plot
figure(1)
plot(elps_t, ax,'.')
xlabel('time, s')
ylabel('ax, ?')
grid on
title("Accel. x-axis")

figure(2)
plot(elps_t, ay,'.')
xlabel('time, s')
ylabel('ay, ?')
grid on
title("Accel. y-axis")

figure(3)
plot(elps_t, az,'.')
xlabel('time, s')
ylabel('az, ?')
grid on
title("Accel. z-axis")

figure(4)
plot(elps_t, wx,'.')
xlabel('time, s')
ylabel('wx, deg/s')
grid on
title("Gyro x-axis")

figure(5)
plot(elps_t, wy,'.')
xlabel('time, s')
ylabel('wy, deg/s')
grid on
title("Gyro y-axis")

figure(6)
plot(elps_t, wz,'.')
xlabel('time, s')
ylabel('wz, deg/s')
grid on
title("Gyro z-axis")


%% Utility function

% Parsing one line of IMU data string
function [line] = parse_line(line_str)
    line_str = convertStringsToChars(line_str);

    % naive check
    if(line_str(3)~= ':' || line_str(6)~= ':')
        error("ERROR: Check code or data for problem: The format is not suitable to use this parser.")
    end
    

    % parsing
    splt_lines = split(line_str,','); % split line to phrases
    hms_str = splt_lines{1}; % time phrase string
    colons_idx = strfind(hms_str,':');
    
    line.hrs = str2double(hms_str(1:colons_idx(1)-1)); % hours
    line.min = str2double(hms_str(colons_idx(1)+1:colons_idx(2)-1)); % minutes
    line.sec =  str2double(hms_str(colons_idx(2)+1:end)); % seconds
    line.cnt = str2double(splt_lines(2)); % counts
    line.ax = str2double(splt_lines(3)); % acceleration x
    line.ay = str2double(splt_lines(4)); % acceleration y
    line.az = str2double(splt_lines(5)); % acceleration z
    line.wx = str2double(splt_lines(6)); % gyro x
    line.wy = str2double(splt_lines(7)); % gyro y
    line.wz = str2double(splt_lines(8)); % gyro z  
end