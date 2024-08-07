
T_Environment=readtable('C:\Users\shijusis\OneDrive - Michigan Medicine\Desktop\Michigan Behavior Analysis BL GW\Monitor51.txt');
T_Environment_cleaned=table2array(T_Environment(:, end-31:end));


Time_Environment=linspace(0,length(T_Environment_cleaned)/(60*24),length(T_Environment_cleaned));

title('Environmental conditions Run-M57')
yyaxis left
plot(Time_Environment,T_Environment_cleaned(:,4));
ylabel('Average light (lux)')

yyaxis right
plot(Time_Environment,T_Environment_cleaned(:,9)/10);
ylabel('Average temperature (degC)')
xlabel('Days')


path='C:\Users\shijusis\OneDrive - Michigan Medicine\Desktop\Michigan Behavior Analysis BL GW'
folder = mkdir([path,filesep,'Environment_condition']);
path  = [path,filesep,'Environment_condition'] ;

figure(1);
        temp=[path,filesep,'M57','.png'];

saveas(gca,temp);
