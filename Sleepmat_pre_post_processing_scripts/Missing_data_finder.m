

Monitor_file='18';

%EV_path=['C:\Users\shijusis\OneDrive - Michigan Medicine\Desktop\Shiju_sisobhan\GitHub_folder\Allada_Lab_University_of_Michigan\sleepmat_Windows\Example','\Monitor',Monitor_file,'.txt'];

EV_path=['C:\Users\shijusis\OneDrive - Michigan Medicine\Desktop\Re_ Sleep architecture_rebound for paper\GW\M30','\Monitor',Monitor_file,'.txt'];


X_fly_Raw=importdata(EV_path);

%time_col = datetime(X_fly_Raw.textdata(:,3), 'Format', 'HH:mm:ss');          % Time column (col3)
time_col = datetime(X_fly_Raw.textdata(:,3), 'Format', 'HH:mm');          % Time column (col3)

time_col.Second=0; % set all the seconds value to zero (Sometimes it was not automatically)


try
datetime_combined = datetime(X_fly_Raw.textdata(:,2), 'InputFormat', 'dd MMM yy') + timeofday(datetime(time_col));
catch
    try
    datetime_combined = datetime(X_fly_Raw.textdata(:,2)) + timeofday(datetime(time_col));
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

  [~, idx_in_complete] = ismember(datetime_combined,complete_timeline); % change here


  new_data(idx_in_complete,:)=X_fly_Raw.data;

  new_date_col=cellstr(datestr(complete_timeline, 'dd mmm yy'));
  new_time_col = cellstr(datestr(complete_timeline, 'HH:MM:ss'));
  new_arbitrary_number= num2cell(zeros(length(complete_timeline),1));

new_textdata = [new_arbitrary_number, new_date_col, new_time_col];

X_fly_Raw.data = new_data;
X_fly_Raw.textdata = new_textdata;


end

