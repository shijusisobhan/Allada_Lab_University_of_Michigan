
% Run_number='M63';
% Environment_file='51';

Run_number='M1';
Environment_file='860';

% EV_path=['C:\Users\shijusis\OneDrive - Michigan Medicine\Desktop\Michigan Behavior Analysis BL GW\GW_new_8_15\',Run_number,'\Monitor',Environment_file,'.txt'];
EV_path=['C:\Users\shijusis\OneDrive - Michigan Medicine\Desktop\Tuhin\Environment','\Monitor',Environment_file,'.txt'];

T_Environment=readtable(EV_path);
T_Environment_cleaned=table2array(T_Environment(:, end-31:end));


Time_Environment=linspace(0,length(T_Environment_cleaned)/(60*24),length(T_Environment_cleaned));

title(['Environmental conditions Run-',Run_number])
yyaxis left
plot(Time_Environment,T_Environment_cleaned(:,4));
ylabel('Average light over min (lux)')

yyaxis right
plot(Time_Environment,T_Environment_cleaned(:,9)/10);
ylabel('Average temperature over min (degC)')
xlabel('Days')


% path='C:\Users\shijusis\OneDrive - Michigan Medicine\Desktop\Michigan Behavior Analysis BL GW\GW_new_8_15\Environment_condition_plot';

path='C:\Users\shijusis\OneDrive - Michigan Medicine\Desktop\Tuhin\Environment';


figure(1);
        temp=[path,filesep,'Environment_chart_',Run_number,'.png'];

saveas(gca,temp);
