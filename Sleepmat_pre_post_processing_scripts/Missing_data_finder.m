

Monitor_file='18';

%EV_path=['C:\Users\shijusis\OneDrive - Michigan Medicine\Desktop\Shiju_sisobhan\GitHub_folder\Allada_Lab_University_of_Michigan\sleepmat_Windows\Example','\Monitor',Monitor_file,'.txt'];

EV_path=['C:\Users\shijusis\OneDrive - Michigan Medicine\Desktop\Re_ Sleep architecture_rebound for paper\GW\M30','\Monitor',Monitor_file,'.txt'];


T_Environment=importdata(EV_path);

% Convert Date and Time to datetime format
try
datetime_combined = datetime(T_Environment.textdata(:,2), 'InputFormat', 'dd MMM yy') + timeofday(datetime(T_Environment.textdata(:,3)));
catch
    try
    datetime_combined = datetime(T_Environment.textdata(:,2)) + timeofday(datetime(T_Environment.textdata(:,3)));
    catch
disp('ERROR! Date and time format is not standard, check Monitor file')
      set(handles.M_box,'String',[oldmsgs;{'ERROR! Date and time format is not standard, check Monitor file'}] );drawnow
       return
    end
end


% Create a complete timeline with 1-minute intervals
start_time = datetime_combined(1);
end_time = datetime_combined(end);
complete_timeline = (start_time:minutes(1):end_time)';

% Find missing times
missing_times = setdiff(complete_timeline, datetime_combined);

if isempty(missing_times)
    fprintf('No missing values found')
else
    writematrix(missing_times, "Missing_times.xls")
    disp('ERROR! Missing data found')
      %set(handles.M_box,'String',[oldmsgs;{'ERROR! Missing data found'}] );drawnow
       return
end

% Display the missing times
disp('Missing times:');
disp(missing_times);