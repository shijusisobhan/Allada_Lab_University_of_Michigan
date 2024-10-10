classdef sleepmat < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        figure1                matlab.ui.Figure
        uipanel17              matlab.ui.container.Panel
        text44                 matlab.ui.control.Label
        text43                 matlab.ui.control.Label
        M_box                  matlab.ui.control.TextArea
        uipanel15              matlab.ui.container.Panel
        cum_fig                matlab.ui.container.ButtonGroup
        cum_No                 matlab.ui.control.RadioButton
        cum_Yes                matlab.ui.control.RadioButton
        text48                 matlab.ui.control.Label
        SG_int                 matlab.ui.control.EditField
        SD_int                 matlab.ui.control.EditField
        text47                 matlab.ui.control.Label
        Triaged_cndn           matlab.ui.container.ButtonGroup
        Threshold3             matlab.ui.control.EditField
        Threshold2             matlab.ui.control.EditField
        Threshold4             matlab.ui.control.EditField
        SL_ch3                 matlab.ui.control.CheckBox
        SL_ch2                 matlab.ui.control.CheckBox
        SL_ch1                 matlab.ui.control.CheckBox
        Triage_date            matlab.ui.container.ButtonGroup
        T_date_No              matlab.ui.control.RadioButton
        T_date_Yes             matlab.ui.control.RadioButton
        SD_per                 matlab.ui.control.EditField
        text12                 matlab.ui.control.Label
        SD_Th                  matlab.ui.container.ButtonGroup
        SD_Th_No               matlab.ui.control.RadioButton
        SD_Th_Yes              matlab.ui.control.RadioButton
        SD_dur                 matlab.ui.control.EditField
        text11                 matlab.ui.control.Label
        SD_hr                  matlab.ui.control.EditField
        text10                 matlab.ui.control.Label
        SD_day                 matlab.ui.control.EditField
        text9                  matlab.ui.control.Label
        Sleep_dep              matlab.ui.container.ButtonGroup
        SD_No                  matlab.ui.control.RadioButton
        SD_Yes                 matlab.ui.control.RadioButton
        per_Interval           matlab.ui.control.EditField
        text21                 matlab.ui.control.Label
        analysis_per_interval  matlab.ui.container.ButtonGroup
        Int_No                 matlab.ui.control.RadioButton
        Int_Yes                matlab.ui.control.RadioButton
        N_days_ana             matlab.ui.control.EditField
        text2                  matlab.ui.control.Label
        st_day_sleep           matlab.ui.control.EditField
        text32                 matlab.ui.control.Label
        sleep_ana              matlab.ui.container.ButtonGroup
        SL_No                  matlab.ui.control.RadioButton
        SL_Yes                 matlab.ui.control.RadioButton
        uipanel16              matlab.ui.container.Panel
        uibuttongroup27        matlab.ui.container.ButtonGroup
        Threshold4_per         matlab.ui.control.EditField
        Threshold2_per         matlab.ui.control.EditField
        Threshold3_per         matlab.ui.control.EditField
        per_ch3                matlab.ui.control.CheckBox
        per_ch2                matlab.ui.control.CheckBox
        per_ch1                matlab.ui.control.CheckBox
        text27                 matlab.ui.control.Label
        text26                 matlab.ui.control.Label
        Triaged_cndn_ed        matlab.ui.container.ButtonGroup
        Threshold4_ed          matlab.ui.control.EditField
        Threshold2_ed          matlab.ui.control.EditField
        Threshold3_ed          matlab.ui.control.EditField
        ed_ch3                 matlab.ui.control.CheckBox
        ed_ch2                 matlab.ui.control.CheckBox
        ed_ch1                 matlab.ui.control.CheckBox
        Period_figure          matlab.ui.container.ButtonGroup
        PF_No                  matlab.ui.control.RadioButton
        PF_Yes                 matlab.ui.control.RadioButton
        text38                 matlab.ui.control.Label
        N_act                  matlab.ui.control.EditField
        st_act                 matlab.ui.control.EditField
        text37                 matlab.ui.control.Label
        Acto                   matlab.ui.container.ButtonGroup
        act_No                 matlab.ui.control.RadioButton
        act_Yes                matlab.ui.control.RadioButton
        N_per                  matlab.ui.control.EditField
        text31                 matlab.ui.control.Label
        Triage_date_per        matlab.ui.container.ButtonGroup
        T_date_No_per          matlab.ui.control.RadioButton
        T_date_Yes_per         matlab.ui.control.RadioButton
        max_period             matlab.ui.control.EditField
        text34                 matlab.ui.control.Label
        min_period             matlab.ui.control.EditField
        min                    matlab.ui.control.Label
        sig_val                matlab.ui.control.EditField
        sig_level              matlab.ui.control.Label
        periodogram            matlab.ui.container.ButtonGroup
        periodo_No             matlab.ui.control.RadioButton
        periodo_Yes            matlab.ui.control.RadioButton
        Triage_date_ed         matlab.ui.container.ButtonGroup
        T_date_ed_No           matlab.ui.control.RadioButton
        T_date_ed_Yes          matlab.ui.control.RadioButton
        start_per              matlab.ui.control.EditField
        text30                 matlab.ui.control.Label
        N_ed_anti              matlab.ui.control.EditField
        st_ed_anti             matlab.ui.control.EditField
        Eduction               matlab.ui.container.ButtonGroup
        Eduction_No            matlab.ui.control.RadioButton
        Eduction_Yes           matlab.ui.control.RadioButton
        Anticipation           matlab.ui.container.ButtonGroup
        Anti_No                matlab.ui.control.RadioButton
        Anti_Yes               matlab.ui.control.RadioButton
        Ed_plot                matlab.ui.container.ButtonGroup
        Ed_No                  matlab.ui.control.RadioButton
        Ed_Yes                 matlab.ui.control.RadioButton
        uipanel14              matlab.ui.container.Panel
        Out_files              matlab.ui.container.ButtonGroup
        Out_No                 matlab.ui.control.RadioButton
        Out_Yes                matlab.ui.control.RadioButton
        Max_LF                 matlab.ui.control.EditField
        LF_th                  matlab.ui.control.EditField
        lifespan               matlab.ui.container.ButtonGroup
        LF_No                  matlab.ui.control.RadioButton
        LF_Yes                 matlab.ui.control.RadioButton
        text46                 matlab.ui.control.Label
        lifespan_th            matlab.ui.control.Label
        text35                 matlab.ui.control.Label
        project_name           matlab.ui.control.EditField
        P_name                 matlab.ui.control.Label
        t_Loff1                matlab.ui.control.EditField
        text5                  matlab.ui.control.Label
        t_LON                  matlab.ui.control.EditField
        text4                  matlab.ui.control.Label
        File_Name              matlab.ui.control.EditField
        open_file              matlab.ui.control.Button
        Start                  matlab.ui.control.Button
    end

    
    methods (Access = private)
        function SD_sleep_1hr = Calculate_sleep_1hr(app, Data_SD, N_ch, Triaged_matrix)
            
            N_inval_min=60;
            
            [nr, nc]=size(Data_SD);
            
            sleep_matrix_inval=zeros(floor(nr/N_inval_min), N_ch);
            
            for i_seg=1:floor(nr/N_inval_min)
                data_seg=Data_SD(((i_seg-1)*N_inval_min+1):i_seg*N_inval_min,:);%get the data with in the interval
            
                for i_colm=1:N_ch % for first chanel of selection
                    N_count=0;
            
                    if Triaged_matrix(i_colm)== 0
            
            
                sleep_matrix_inval(:, i_colm)=NaN;
            
            
                    else
                    for i_row=1:length(data_seg)
                       try
                        if data_seg(i_row,i_colm)==0
                            N_count=N_count+1;
                        elseif data_seg(i_row,i_colm)>0 & N_count>=5 % condition for sleep detection
                            sleep_matrix_inval(i_seg, i_colm)=sleep_matrix_inval(i_seg,i_colm)+N_count;% for sleep
            
            
            
                            N_count=0;
            
                           % Check zeros in the first rows and is it a continuation from previous segment
                       elseif (data_seg(i_row,i_colm)>0 & N_count<5) & Data_SD((((i_seg-1)*N_inval_min+i_row)-5):(((i_seg-1)*N_inval_min+i_row)-1),i_colm)==0 % condition for sleep detection
                            sleep_matrix_inval(i_seg, i_colm)=sleep_matrix_inval(i_seg,i_colm)+N_count;% for bout length
                            N_count=0;
            
                        else
                            N_count=0;
                        end
                       catch
                           N_count=0;
                   end
                    end
                    try
                         % Check the zeros are in last entries OR zeros are initiate in this segment and continues in next segment
                    if (data_seg(i_row,i_colm)==0 & N_count>=5)  | (N_count>0 & N_count<5 & Data_SD((((i_seg-1)*N_inval_min+i_row)-(N_count-1)):(((i_seg-1)*N_inval_min+i_row)+(5-N_count)),i_colm)==0)  % condition for sleep if last element of chnel is zero
                            sleep_matrix_inval(i_seg, i_colm)=sleep_matrix_inval(i_seg,i_colm)+N_count;
            
            
            
                    end
                    catch
                    end
            
                end
                end
            end
            
            SD_sleep_1hr=sleep_matrix_inval';
        end
        
        function Start_hr_CreateFcn(app, hObject, eventdata, handles)
            % --- Executes during object creation, after setting all properties.
            
            % hObject    handle to Start_hr (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    empty - handles not created until after all CreateFcns called
            
            % Hint: edit controls usually have a white background on Windows.
            %       See ISPC and COMPUTER.
            if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
                set(hObject,'BackgroundColor','white');
            end
        end
        
        function Start_min_CreateFcn(app, hObject, eventdata, handles)
            % --- Executes during object creation, after setting all properties.
            
            % hObject    handle to Start_min (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    empty - handles not created until after all CreateFcns called
            
            % Hint: edit controls usually have a white background on Windows.
            %       See ISPC and COMPUTER.
            if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
                set(hObject,'BackgroundColor','white');
            end
        end
        
        function Threshold1_CreateFcn(app, hObject, eventdata, handles)
            % --- Executes during object creation, after setting all properties.
            
            % hObject    handle to Threshold1 (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    empty - handles not created until after all CreateFcns called
            
            % Hint: edit controls usually have a white background on Windows.
            %       See ISPC and COMPUTER.
            if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
                set(hObject,'BackgroundColor','white');
            end
        end
        
        function Threshold1_ed_CreateFcn(app, hObject, eventdata, handles)
            % --- Executes during object creation, after setting all properties.
            
            % hObject    handle to Threshold1_ed (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    empty - handles not created until after all CreateFcns called
            
            % Hint: edit controls usually have a white background on Windows.
            %       See ISPC and COMPUTER.
            if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
                set(hObject,'BackgroundColor','white');
            end
        end
        
        function Threshold1_per_CreateFcn(app, hObject, eventdata, handles)
            % --- Executes during object creation, after setting all properties.
            
            % hObject    handle to Threshold1_per (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    empty - handles not created until after all CreateFcns called
            
            % Hint: edit controls usually have a white background on Windows.
            %       See ISPC and COMPUTER.
            if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
                set(hObject,'BackgroundColor','white');
            end
        end
        
        function Qp = calc_Qp(app, x, y, minpoints, maxpoints, variance)
            %Qp calculation
            
             len_Qp=maxpoints-minpoints+1;
             Qp=zeros(1,len_Qp);
            
             npoints=minpoints:maxpoints;
             for k=1:length(npoints)
             nblocks=floor(length(x)/npoints(k));
             table=zeros(nblocks,npoints(k));
             for i=1:nblocks
             for j=1:npoints(k)
             pos=j+npoints(k)*(i-1);
             table(i,j)=y(pos);
             end
             end
             if (len_Qp==1)
             Qp=npoints(k)*nblocks*var(mean(table,1))/variance;
             else
             Qp(k)=npoints(k)*nblocks*var(mean(table,1))/variance;
             end
             end
        end
        
        function [sig_vector, Qp, sig_line] = chi2_periodogram(app, Time_data, Activity_data, noscip, maxperm, func, plts, geno_T, alpha_val, per_min, per_max)
            
            min_per=per_min*60; % Minimum period to evaluate
            max_per=per_max*60; % maximum period to evaluate
            
            sig_vector=per_min:0.5:per_max;
            sig_line= chi2inv(alpha_val,sig_vector*2); % for significant line with alpha=0.01
            
            per_rslt=zeros(1,noscip*3); % Final Results
            
            minpoints=length(find(Time_data<=min_per));
            maxpoints=length(find(Time_data<=max_per));
            Qp=calc_Qp(app, Time_data,Activity_data,minpoints,maxpoints,var(Activity_data));
            xaxis1=(minpoints:maxpoints)*(Time_data(2)-Time_data(1));
            xaxis=xaxis1/60; % convert minute into hrs
            
            
             plot(xaxis,Qp,'k');
             set(gca,'fontweight','b','fontsize',16);
             xlabel('Period (h)','fontweight','b','fontsize',16);
             ylabel('Qp','fontweight','b','fontsize',16);
             title(geno_T)
            
             hold on
             plot(sig_vector, sig_line)
        end
        
        function cum_SL = cum_sleep_loss(app, SDL_data1, SDL_data2, SDG_data1, SDG_data2, N_ch, GT_name, Triaged_matrix)
            
            
            SDL_data1_sleep=Calculate_sleep_1hr(app, SDL_data1,N_ch,Triaged_matrix);
            SDL_data2_sleep=Calculate_sleep_1hr(app, SDL_data2,N_ch,Triaged_matrix);
            
            SDG_data1_sleep=Calculate_sleep_1hr(app, SDG_data1,N_ch,Triaged_matrix);
            SDG_data2_sleep=Calculate_sleep_1hr(app, SDG_data2,N_ch,Triaged_matrix);
            
            duration=(length(SDL_data1)+length(SDG_data1))/60;
            
            
            
            SDL2_sleep=sum(SDL_data2_sleep,2);
            SDG2_sleep=sum(SDG_data2_sleep,2);
            
            sleep_loss=(SDL_data1_sleep)-(SDL_data2_sleep);
            sleep_gain=(SDG_data1_sleep)-(SDG_data2_sleep);
            sleep_loss_gain=[sleep_loss,sleep_gain];
            
            %cum_sleep_loss_gain_per=cumsum((sleep_loss_gain./SDL2_sleep)*100,2);
            
            cum_sleep_loss_gain_per=cumsum((sleep_loss_gain),2);
            
            
            plot(linspace(1,duration,duration),cum_sleep_loss_gain_per,'Color',[.7 .7 .7]);
            hold on;
            plot(linspace(1,duration,duration),mean(cum_sleep_loss_gain_per,"omitmissing"),'r-','linewidth',2);
            ylabel('sleep lost (min)')
            xlabel('hr')
            title(GT_name)
        end
        
        function SM_act = sleepmat_actogram(app, Data_seg_all, i_nch, N_days_actogram, geno_type, t_Lon, Monitor_number, Ch_num, path)
            
            
            % for figure saving
            
            
            
            
            % x axis labels
            x1=mod(t_Lon,24) ;x2=mod(t_Lon+8,24);x3=mod(t_Lon+16,24);x4=mod(t_Lon+24,24);
            x5=mod(t_Lon+32,24);x6=mod(t_Lon+40,24);
            
            
            Data_seg_ch=Data_seg_all(:,i_nch);
            L_ch=length(Data_seg_ch);
            Act_data1=[];
            
            
            
            
            
             for i_n=1:N_days_actogram
                 Act_start=1+(i_n-1)*(24)*60;
                 try
                 Act_data1=[Act_data1 Data_seg_ch(Act_start:Act_start+(48*60)-1,:)];
                 catch
                     try
                     A_zer_pad = [Data_seg_ch(Act_start:end,:);zeros(L_ch-Act_start+1,1)];
                         Act_data1=[Act_data1 A_zer_pad];% Zero paddding
                     catch
                         Act_data1=[Act_data1 zeros(48*60,1)];
                     end
                 end
            
             end
            
             Act_data=zeros(96,N_days_actogram);
            
                     for i_act=1:96 % 96 --> 30 min bin for 48hr
            
                             Act_data(i_act, :)=sum(Act_data1(((i_act-1)*30)+1:i_act*30,:)); %30 min bin
            
                     end
            
            
            
             t_plot=0:1/(60):48;
             [nr_act, nc_act]=size(Act_data);
            
            
             %% ploting actogram
            num_act = nc_act;
            start_Pos = fliplr(linspace(.04, .90, num_act)); % Position of each days (y axis)
            h_act = start_Pos(1) - start_Pos(2) - .005; % hight of actogram (make sure it is not overlaping)
            act_max =max(max(Act_data)); % find the maximum value of actogram
            
            
            Title_plot=[geno_type, '  ' 'Monitor= ' num2str(Monitor_number) '  '  'Ch= ' num2str(Ch_num)];
            
            
            
            
            if act_max>0% To avoide dead flies
                figure('NumberTitle', 'off', 'Name', Title_plot);
                for n_p = 1:num_act
                    axes('Position',[.05 start_Pos(n_p) .90 h_act],'Color','none',...
                        'XTick',[],'YTick',[n_p]);
                    hold on;
                    h_ag=bar(Act_data(:,n_p), 'BarWidth', 1);
                    set(gca,'YLim',[0 act_max], 'XLim', [0 96]);
                    set(h_ag,'FaceColor','k')
                end
            
                set(gca, 'XTick',0:16:96, 'XTickLabel',[x1,x2,x3,x4,x5,x6]);
            
            
            
                temp1=[path,filesep,Title_plot,'.png'];
                saveas(gca,temp1);
                close(Title_plot);
            
            end
        end
        
        function start_date_ana_CreateFcn(app, hObject, eventdata, handles)
            % --- Executes during object creation, after setting all properties.
            
            % hObject    handle to start_date_ana (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    empty - handles not created until after all CreateFcns called
            
            % Hint: edit controls usually have a white background on Windows.
            %       See ISPC and COMPUTER.
            if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
                set(hObject,'BackgroundColor','white');
            end
        end
        
        function start_date_exp_CreateFcn(app, hObject, ~, handles)
            % --- Executes during object creation, after setting all properties.
            
            % hObject    handle to start_date_exp (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    empty - handles not created until after all CreateFcns called
            
            % Hint: edit controls usually have a white background on Windows.
            %       See ISPC and COMPUTER.
            if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
                set(hObject,'BackgroundColor','white');
            end
        end
        
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function sleepmat_OpeningFcn(app, varargin)
            % --- Executes just before sleepmat is made visible.
            
            % Ensure that the app appears on screen when run
            movegui(app.figure1, 'onscreen');
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app); %#ok<ASGLU>
            
            % This function has no output args, see OutputFcn.
            % hObject    handle to figure
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            % varargin   command line arguments to sleepmat (see VARARGIN)
            
            % Choose default command line output for sleepmat
            handles.output = hObject;
            
            % Update handles structure
            guidata(hObject, handles);
        end

        % Selection changed function: Acto
        function Acto_SelectionChangedFcn(app, event)
            % --- Executes when selected object is changed in Acto.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            
            
            
            Acto_select1=get(handles.Acto, 'SelectedObject');
            Acto_select2=get(Acto_select1, 'String');
            
            
            switch Acto_select2
            
                case 'Yes'
            set(handles.st_act, 'Enable', 'on');
             set(handles.N_act, 'Enable', 'on');
            
               case 'No'
                   set(handles.st_act, 'Enable', 'off');
             set(handles.N_act, 'Enable', 'off');
            end
        end

        % Selection changed function: Anticipation
        function Anticipation_SelectionChangedFcn(app, event)
            % --- Executes when selected object is changed in Anticipation.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            
            Anticipation_select1=get(handles.Anticipation, 'SelectedObject');
            Anticipation_select2=get(Anticipation_select1, 'String');
        end

        % Selection changed function: Ed_plot
        function Ed_plot_SelectionChangeFcn(app, event)
            % --- Executes when selected object is changed in Ed_plot. (circadian
            % analysis selection)
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            
            Ed_select1=get(handles.Ed_plot, 'SelectedObject');
            Ed_select2=get(Ed_select1, 'String');
            
            
            switch Ed_select2
                case 'Yes'
                   set(handles.st_ed_anti, 'Enable', 'on');
                   set(handles.N_ed_anti, 'Enable', 'on');
            %        set(handles.Ind_Ed_No, 'Enable', 'on');
            %        set(handles.Ind_Ed_Yes, 'Enable', 'on');
                   set(handles.Anti_Yes, 'Enable', 'on');
                   set(handles.Anti_No, 'Enable', 'on');
                   set(handles.Eduction_Yes, 'Enable', 'on');
                   set(handles.Eduction_No, 'Enable', 'on');
            
            
            
                     set(handles.T_date_ed_No, 'Enable', 'on');
                     set(handles.T_date_ed_Yes, 'Enable', 'on');
                      set(handles.Threshold3_ed, 'Enable', 'on');
                     set(handles.Threshold2_ed, 'Enable', 'on');
                     set(handles.Threshold4_ed, 'Enable', 'on');
            
            
                     set(handles.act_Yes, 'Enable', 'on');
                     set(handles.act_No, 'Enable', 'on');
            
            
                   set(handles.ed_ch1, 'Enable', 'on');
                   set(handles.ed_ch2, 'Enable', 'on');
                   set(handles.ed_ch3, 'Enable', 'on')
            
            
                case 'No'
                   set(handles.st_ed_anti, 'Enable', 'off');
                   set(handles.N_ed_anti, 'Enable', 'off');
            %        set(handles.Ind_Ed_No, 'Enable', 'off');
            %        set(handles.Ind_Ed_Yes, 'Enable', 'off');
                   set(handles.Anti_Yes, 'Enable', 'off');
                   set(handles.Anti_No, 'Enable', 'off');
                   set(handles.Eduction_Yes, 'Enable', 'off');
                   set(handles.Eduction_No, 'Enable', 'off');
            
            
            
            
                     set(handles.T_date_ed_No, 'Enable', 'off');
                     set(handles.T_date_ed_Yes, 'Enable', 'off');
                     set(handles.Threshold3_ed, 'Enable', 'off');
                     set(handles.Threshold2_ed, 'Enable', 'off');
                     set(handles.Threshold4_ed, 'Enable', 'off');
                     %set(handles.Threshold4_ed, 'Enable', 'off');
                     set(handles.act_Yes, 'Enable', 'off');
                     set(handles.act_No, 'Enable', 'off');
            
                    set(handles.ed_ch1, 'Enable', 'off');
                   set(handles.ed_ch2, 'Enable', 'off');
                   set(handles.ed_ch3, 'Enable', 'off')
            
             end
        end

        % Selection changed function: Eduction
        function Eduction_SelectionChangedFcn(app, event)
            % --- Executes when selected object is changed in Eduction.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            
            Eduction_select1=get(handles.Eduction, 'SelectedObject');
            Eduction_select2=get(Eduction_select1, 'String');
        end

        % Selection changed function: Out_files
        function Out_files_SelectionChangedFcn(app, event)
            % --- Executes when selected object is changed in Out_files.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            
            Outfile_select1=get(handles.Out_files, 'SelectedObject');
            Outfile_select2=get(Outfile_select1, 'String');
        end

        % Selection changed function: Period_figure
        function Period_figure_SelectionChangedFcn(app, event)
            % --- Executes when selected object is changed in Period_figure.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            
            periodo_fig_select1=get(handles.Period_figure, 'SelectedObject');
            periodo_fig_select2=get(periodo_fig_select1, 'String');
        end

        % Selection changed function: SD_Th
        function SD_Th_SelectionChangeFcn(app, event)
            % --- Executes when selected object is changed in SD_Th.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            
            
            SD_Th_select1=get(handles.SD_Th, 'SelectedObject');
            SD_Th_select2=get(SD_Th_select1, 'String');
            
            
            switch SD_Th_select2
                case 'Yes'
                   set(handles.SD_per, 'Enable', 'on');
                case 'No'
                   set(handles.SD_per, 'Enable', 'off');
             end
        end

        % Selection changed function: Sleep_dep
        function Sleep_dep_SelectionChangeFcn(app, event)
            % --- Executes when selected object is changed in Sleep_dep.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            
            
            SD_select1=get(handles.Sleep_dep, 'SelectedObject');
            SD_select2=get(SD_select1, 'String');
            
            %handles.SD_select_temp=[handles.SD_select,{SD_select2}]
            
            
            switch SD_select2
                case 'Yes'
                   set(handles.SD_day, 'Enable', 'on');
                   set(handles.SD_hr, 'Enable', 'on');
                   set(handles.SD_dur, 'Enable', 'on');
                   set(handles.SD_per, 'Enable', 'on');
                   set(handles.SD_Th_No, 'Enable', 'on');
                   set(handles.SD_Th_Yes, 'Enable', 'on');
                   set(handles.Int_Yes, 'Enable', 'on');
                   set(handles.SD_int, 'Enable', 'on');
                   set(handles.SG_int, 'Enable', 'on');
                   set(handles.cum_Yes, 'Enable', 'on');
                   set(handles.cum_No, 'Enable', 'on');
            
            
            
                case 'No'
                   set(handles.SD_day, 'Enable', 'off');
                   set(handles.SD_hr, 'Enable', 'off');
                   set(handles.SD_dur, 'Enable', 'off');
                   set(handles.SD_per, 'Enable', 'off');
                   set(handles.SD_Th_No, 'Enable', 'off');
                   set(handles.SD_Th_Yes, 'Enable', 'off');
                   set(handles.SD_int, 'Enable', 'off');
                   set(handles.SG_int, 'Enable', 'off');
                   set(handles.cum_Yes, 'Enable', 'off');
                   set(handles.cum_No, 'Enable', 'off');
            end
        end

        % Button pushed function: Start
        function Start_Callback(app, event)
            % --- Executes on button press in Start.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            
            
             try
            
            oldmsgs = cellstr(get(handles.M_box,'String'));
            
            
            
            
            
            
            Outfile_select1=get(handles.Out_files, 'SelectedObject');
            Outfile_select2=get(Outfile_select1, 'String');
            %% ********************************************
            
            
            SL_select1=get(handles.sleep_ana, 'SelectedObject');
            SL_select2=get(SL_select1, 'String');
            
            SD_select1=get(handles.Sleep_dep, 'SelectedObject');
            SD_select2=get(SD_select1, 'String');
            
            
            
            
            
            
            % strt_time_spec_select1=get(handles.strat_time_spec, 'SelectedObject');
            % strt_time_spec_select2=get(strt_time_spec_select1, 'String');
            
            Interval_select1=get(handles.analysis_per_interval, 'SelectedObject');
            Interval_select2=get(Interval_select1, 'String');
            
            SD_Th_select1=get(handles.SD_Th, 'SelectedObject');
            SD_Th_select2=get(SD_Th_select1, 'String');
            
            cum_fig_select1=get(handles.cum_fig, 'SelectedObject');
            cum_fig_select2=get(cum_fig_select1, 'String');
            
            % anti_day_select1=get(handles.M_anti_calculation, 'SelectedObject');
            % anti_day_select=get(anti_day_select1, 'String');
            
            
            Ed_select1=get(handles.Ed_plot, 'SelectedObject');% This is for circadian analysis selection
            Ed_select2=get(Ed_select1, 'String');
            
            
            Anticipation_select1=get(handles.Anticipation, 'SelectedObject');
            Anticipation_select2=get(Anticipation_select1, 'String');
            
            Eduction_select1=get(handles.Eduction, 'SelectedObject');
            Eduction_select2=get(Eduction_select1, 'String');
            
            
            
            
            
            
            % Ind_Ed_select1=get(handles.Ind_Ed_save, 'SelectedObject');
            % Ind_Ed_select2=get(Ind_Ed_select1, 'String');
            
            
            periodo_select1=get(handles.periodogram, 'SelectedObject');
            periodo_select2=get(periodo_select1, 'String');
            
            
            T_date_select1=get(handles.Triage_date, 'SelectedObject');
            T_date_select2=get(T_date_select1, 'String');
            
            
            T_date_ed_select1=get(handles.Triage_date_ed, 'SelectedObject');
            T_date_ed_select2=get(T_date_ed_select1, 'String');
            
            
            T_date_per_select1=get(handles.Triage_date_per, 'SelectedObject');
            T_date_per_select2=get(T_date_per_select1, 'String');
            
            
            Acto_select1=get(handles.Acto, 'SelectedObject');
            Acto_select2=get(Acto_select1, 'String');
            
            periodo_fig_select1=get(handles.Period_figure, 'SelectedObject');
            periodo_fig_select2=get(periodo_fig_select1, 'String');
            
            
            Lifespan_select1=get(handles.lifespan, 'SelectedObject');
            Lifespan_select2=get(Lifespan_select1, 'String');
            
            
            
            
            %file_path1 = handles.path;
            %diary(file_path1);
            %diary on;
            
            
            
            warning('off');
            
            disp('Sleep analysis started...')
            
            set(handles.M_box,'String',[oldmsgs;{'Sleep analysis started...'}] );drawnow
            
            
            
            
            
            
            
            
            if isfield(handles,'Fullpath')
                TXT2=readtable(handles.Fullpath);
            end
            
            TXT1=table2cell(TXT2);
            
            
            
            
            geno_type1=TXT1(:,3);
            TXT_date=TXT1(:,6);
            
            try
            TXT_triage_date=TXT1(:,7);
            end
            
            try
            ref_genotype=TXT1(:,8);
            ref_genotype=ref_genotype(~cellfun('isempty',ref_genotype));
            end
            
            
            try
            Environment_data=TXT1(:,9);
            %ref_genotype=ref_genotype(~cellfun('isempty',ref_genotype));
            end
            
            
            
            
            
            t_Lon=str2num(get(handles.t_LON, 'string'));
            t_Loff=str2num(get(handles.t_Loff1, 'string'));%light off time (ih hrs)
            
            
            
            starting_hr=t_Lon;
            starting_min=0;
            
            
            
            start_day_sleep=str2num(get(handles.st_day_sleep, 'string'));
            N_days_analysis=str2num(get(handles.N_days_ana, 'string'));
            %N_days_experiment=str2num(get(handles.N_days_experiment, 'string'));
            
            
            SD_start_day=str2num(get(handles.SD_day, 'string'));
            SD_start_hr=str2num(get(handles.SD_hr, 'string'));
            SD_duration=str2num(get(handles.SD_dur, 'string'));
            SD_percentage=str2num(get(handles.SD_per, 'string'));
            
            SD_interval=str2num(get(handles.SD_int, 'string'));
            SG_interval=str2num(get(handles.SG_int, 'string'));
            
            
            Eduction_startday=str2num(get(handles.st_ed_anti, 'string'));
            N_days_Eduction=str2num(get(handles.N_ed_anti, 'string'));
            
            
            
            % Eduction_startday=1;
            % N_days_Eduction=4;
            
            start_day_anticipation=Eduction_startday;
            N_days_anti=N_days_Eduction;
            
            
            N_inval_hr=str2num(get(handles.per_Interval, 'string'));
            N_inval_min=N_inval_hr*60;
            
            Significance_level=str2num(get(handles.sig_val, 'string'));
            per_min=str2num(get(handles.min_period, 'string'));
            per_max=str2num(get(handles.max_period, 'string'));
            
            
            
            
            Project_name=(get(handles.project_name, 'string'));
            
            
            
            %Triage_th1=str2num(get(handles.Threshold1, 'string'));
            Triage_th2=str2num(get(handles.Threshold2, 'string'));
            Triage_th3=str2num(get(handles.Threshold3, 'string'));
            Triage_th4=str2num(get(handles.Threshold4, 'string'));
            
            
            %Triage_th1_ed=str2num(get(handles.Threshold1_ed, 'string'));
            Triage_th2_ed=str2num(get(handles.Threshold2_ed, 'string'));
            Triage_th3_ed=str2num(get(handles.Threshold3_ed, 'string'));
            Triage_th4_ed=str2num(get(handles.Threshold4_ed, 'string'));
            
            
            %Triage_th1_per=str2num(get(handles.Threshold1_per, 'string'));
            Triage_th2_per=str2num(get(handles.Threshold2_per, 'string'));
            Triage_th3_per=str2num(get(handles.Threshold3_per, 'string'));
            Triage_th4_per=str2num(get(handles.Threshold4_per, 'string'));
            
            
            
            st_day_per=str2num(get(handles.start_per, 'string'));
            N_day_per=str2num(get(handles.N_per, 'string'));
            
            start_day_actogram=str2num(get(handles.st_act, 'string'));
            N_days_actogram=str2num(get(handles.N_act, 'string'));
            
            %% *********************************************************
            
            
            t_start=0;
            
            
            Result_all=[];
            Result_all_SC=[];
            %ed_result=[];
            err_result=[];
            summary_all=[];
            ed_result_avg_std_all=[];
            std_error_all=[];
            ed_result_Ndays=[];
            Triage_all=[];
            Periodo_all=[];
            Periodo_all_SC=[];
            Result_all_anti=[];
            Result_all_anti_SC=[];
            screen_data=[];
            screen_data_anticipation=[];
            screen_data_period=[];
            
            
            Sleep_interval=[];
            BoutN_interval=[];
            BoutL_interval=[];
            Activity_interval=[];
            
            
            
            Sleep_interval1=[];
            BoutN_interval1=[];
            BoutL_interval1=[];
            Activity_interval1=[];
            N_active_array=[];
            N_active_array_ed=[];
            %plug_in_result=[];
            Eduction_bar=[];
            
            
            R_all=[];
            R_ed=[];
            R_ed_Ndays=[];
            R_ed_Day_by_Day=[];
            R_ind_periodo=[];
            R_all_L_D=[];
            R_all_anti=[];
            R_all_LNZ=[];
            
            
            Eduction_Ndays=[];
            Eduction_Day_by_Day=[];
            ed_result_avg_std_all=[];
            
            
            sleep_gp=[];
            Bout_N_gp=[];
            Bout_L_gp=[];
            Total_act_gp=[];
            act_count_gp=[];
            Latency_gp=[];
            Manti_gp=[];
            Eanti_gp=[];
            Manti_ph_gp=[];
            Eanti_ph_gp=[];
            Life_gp=[];
            
            
            sleep_lost_gp=[];
            Per_sleep_lost_gp=[];
            sleep_gain_gp=[];
            Per_sleep_gain_gp=[];
            Latency_SD_gp=[];
            
            R_ed_Ind=[];
            Ind_Ed_Normalized_all=[];
            
            
            periodo_data_all=[];
            Result_all_L_D=[];
            
            R_all_ed_TW=[];
            Result_all_ed_TW=[];
            
            periodogram_ploting_data_2=[];
            Environment_Cell=[];
            
            
            
            
            Run_number1=TXT1(:,1) ;
            
            
            %% ******************************************************************************************
                                 %Start of Each genotype analysis starts%
            %% ******************************************************************************************
            
            
            for i_geno=1:length(geno_type1)
            
            
            Run_number=Run_number1{i_geno}
            
            Monitor=cell2mat(TXT1(:,2));
            Monitor_number=(Monitor(i_geno))
            
            geno_type=geno_type1{i_geno}
            
            NUM_st=cell2mat(TXT1(:,4));
            NUM_end=cell2mat(TXT1(:,5));
            
            N_ch_start=NUM_st(i_geno);
            N_ch_stop=NUM_end(i_geno);
            
            N_ch=N_ch_stop-N_ch_start+1;
            
            % Monitor_file=strcat(handles.path, num2string(Monitor_number);
            %
            % X_fly=importdata(Monitor_file);
            
            %%
            
            sheet_name=['Monitor' num2str(Monitor_number) '.txt'];
            Monitor_file=strcat(handles.path, sheet_name);
            X_fly_Raw=importdata(Monitor_file);
            
            
            
            
                %% **************************** Check whether Environment monitor files are present and stores it in a cell **********************************
            
            
                try
                if ~isnan(Environment_data{i_geno})
                    % Store the value from column 2 and column 5 in the same row
                    Environment_Cell = [Environment_Cell; {Run_number, Environment_data{i_geno}}];
                end
                end
            
            %***********************************************************************************************************************************
            
            
            %% ****************   Check for data Missing  %%%%%%*******************************************************************
            
             %path = handles.path;
            folder = mkdir([handles.path,Project_name]);% create folder for all results
            %path_1  = [handles.path,Project_name] ;
            Missing_folder='Missing_data';
            folder2=mkdir([handles.path,Project_name,filesep,Missing_folder]);
            path_2=[handles.path,Project_name,filesep,Missing_folder];% create folder for missing data
            save_location_Missing_error=[path_2,filesep, strcat('Monitor', num2str(Monitor_number), '_Missing_times.xls')] ;
            
            
            
            
            
            
            time_col = datetime(X_fly_Raw.textdata(:,3), 'Format', 'HH:mm:ss');          % Time column (col3)
            time_col.Second=0; % set all the seconds value to zero (Sometimes it was not automatically)
            
            
            % Convert Date and Time to datetime format
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
                fprintf('No missing values found \n')
            else
                 writematrix(missing_times, save_location_Missing_error)
            
            
                [~, idx_in_complete] = ismember(datetime_combined,complete_timeline); % change here
            
            
              new_data(idx_in_complete,:)=X_fly_Raw.data;
            
              new_date_col=cellstr(datestr(complete_timeline, 'dd mmm yy'));
              new_time_col = cellstr(datestr(complete_timeline, 'HH:MM:ss'));
              new_arbitrary_number= num2cell(zeros(length(complete_timeline),1));
            
            new_textdata = [new_arbitrary_number, new_date_col, new_time_col];
            
            X_fly_Raw.data = new_data;
            X_fly_Raw.textdata = new_textdata;
            
            
            end
            
            % Display the missing times
            %disp('Missing times:');
            %disp(missing_times);
            
            X_plug_date=X_fly_Raw.textdata(:,2);% data start from plug in time
            X_plug_time=X_fly_Raw.textdata(:,3);
            
            %% ***********************************************************************************************
            
            start_day_index=[];
            
            starting_date_ana=TXT_date{i_geno};% date which is in GT file
            
             try
            for i_ana=1:length(X_plug_time)
            
            if strcmp(X_plug_time(i_ana), '23:59:00') &  strcmp(X_plug_date(i_ana), starting_date_ana)
               start_day_index(i_geno)=(i_ana);
                break;
            end
            end
            
                  Starting_time_min_sleep=start_day_index(i_geno)-((23*60+59)-(starting_hr*60)-1)+(start_day_sleep-1)*24*60;
                   Starting_time_lifespan=start_day_index(i_geno)-((23*60+59)-(starting_hr*60)-1);
                   Starting_time_min=Starting_time_lifespan;
             catch
            
                  disp('ERROR! Please enter valid dates (See monitor file and genotype specification file)')
                  set(handles.M_box,'String',[oldmsgs;{'ERROR! Please enter valid dates (See monitor file and genotype specification file)'}] );drawnow
                  return
            
             end
            
            
            %% ***********************************************************************************************
            
            Photo_period=t_Loff-t_Lon;
            
            if Photo_period<0
                Photo_period=24+Photo_period;
            end
            
            
            Ending_time_min_sleep=Starting_time_min_sleep +(N_days_analysis*60*24)-1;
            
            Ending_time_triaged_sleep=Ending_time_min_sleep;
            
            Starting_time_anti=Starting_time_min+(start_day_anticipation-1)*60*24;
            Ending_time_anti=Starting_time_anti +(N_days_anti*60*24)-1;
            
            SD_Starting_time_min=(Starting_time_min_sleep-(starting_hr*60)+((SD_start_day-1)*24+SD_start_hr)*60);
            SD_ending_time_min=(SD_Starting_time_min+SD_duration*60)-1;
            
            Starting_time_Eduction_min=Starting_time_min+(Eduction_startday-1)*24*60-((24-Photo_period)/2)*60;
            Ending_time_Eduction_min=Starting_time_Eduction_min +(N_days_Eduction*60*24)-1;
            
            start_time_periodo_min=(st_day_per-1)*24*60+1;
            end_time_periodo_min=start_time_periodo_min+(N_day_per*60*24)-1;
            
            
            Starting_time_actogram_min=Starting_time_min+(start_day_actogram-1)*60*24;
            Ending_time_actogram_min=Starting_time_actogram_min +(N_days_actogram*60*24)-1;
            
            
            
            
            
            
            X_fly=X_fly_Raw.data;
            
            
            
            
            
            %% Find the life span of triage
            
            
            Data_triage=X_fly_Raw.data(:, end-32+N_ch_start:end-32+N_ch_stop);
            
            
            %% *******************************************************************************************************************************
            
            
            try
            Data_all=X_fly_Raw.data(Starting_time_min_sleep:Ending_time_min_sleep, :);% starting from saturday 8AM
            catch
                disp('ERROR! start day of sleep or Number of days of sleep analysis is out of range');
                set(handles.M_box,'String',[oldmsgs;{'ERROR! start day of sleep or Number of days of sleep analysis is out of range'}] );drawnow
                return;
            end
            
            
            
            Data_seg_periodo=X_fly_Raw.data(Starting_time_min:end,end-32+N_ch_start:end-32+N_ch_stop);
            Data_seg_all=X_fly_Raw.data(Starting_time_min_sleep:Ending_time_min_sleep,end-32+N_ch_start:end-32+N_ch_stop);
            
            
            %****************************************************************************************
            
            
            %Titles';
            
            
            Title2a='Sleep analysis result/24hr';
            Title2='Number of days of analysis=';
            Title3='Light ON time=';
            Title4='Light Off time=';
            
            Title5='N=';
            
            
            %Title{1}=[Title0 ' ' string(char(geno_type)) '      ' Title1 ' ' num2str(N_days_experiment) '      ' Title2 ' ' num2str(N_days_analysis) '      ' Title3 ' ' num2str(N_Triaged)];
            
            Title{1}=[Title2a '       ' Title2 ' ' num2str(N_days_analysis) '         '  Title3 ' ' num2str(t_Lon) '       ' Title4 ' ' num2str(t_Loff)];
            
            
            
            % Write the results into excel files
            %Chanel = '';
            
            row_header=[];
            
            for k = 1 : N_ch
                ChanelID{k} = [num2str(N_ch_start+k-1,'%d')];            % Cell Array
            end
            
            
            row_header12=(ChanelID(1:N_ch))';
            row_header=str2num(char(row_header12));
            %L_NZ_str=str2num(char(L_NZ));
            
            
            
            GT=cell(length(row_header),1);
            %GT(:)={string(char(geno_type))};
            
            GT(:)={geno_type};
            
            Run_1=cell(length(row_header),1);
            Run_1(:)={Run_number};
            
            % Run_1=[Run_number*ones(length(row_header),1)];
            
            
            Monitor_1=[Monitor_number*ones(length(row_header),1)];
            
            
            
            
            
            
            
            %% New lifespan
            
            
            if strcmp(Lifespan_select2, 'Yes')
            
            
            LF_th1=str2num(get(handles.LF_th, 'string'));
            LF_days=str2num(get(handles.Max_LF, 'string'));
            
            
            Data_seg_LF=X_fly_Raw.data(Starting_time_min:end,:);
            
            
            Data_triage_date=X_fly_Raw.textdata(Starting_time_min:end,2);
            Data_triage_time=X_fly_Raw.textdata(Starting_time_min:end,3);
            
            L_NZ={};
            L_NZ_day=[];
            
            [nr_LF, nc_LF]=size(Data_seg_periodo);
            
            for i_ch=1:N_ch
            D_FF=Data_seg_periodo(:,i_ch);
            for i_LNZ=1:nr_LF
            
            try
            if sum(D_FF(i_LNZ:24*60+(i_LNZ-1)))<LF_th1
                L_NZ(i_ch,1)=Data_triage_date(i_LNZ);
                L_NZ(i_ch,2)=Data_triage_time(i_LNZ);
                L_NZ_day(i_ch)=((i_LNZ)/(60*24));
            
                if L_NZ_day(i_ch)>LF_days
                L_NZ_day(i_ch)=LF_days;
                end
            
                 break;
            
            end
            
            catch
            
                L_NZ(i_ch,1)=Data_triage_date(end);
                L_NZ(i_ch,2)=Data_triage_time(end);
            
                L_NZ_day(i_ch)=((i_LNZ)/(60*24));
                if L_NZ_day(i_ch)>LF_days
                L_NZ_day(i_ch)=LF_days;
            end
            
            end
            
            
            
            end
            end
            
            
            
            %L_NZ;
            
            R_all_LNZ=[R_all_LNZ;GT Run_1 num2cell(Monitor_1) num2cell(row_header) L_NZ num2cell(L_NZ_day')];
            
            end
            
            
            
            
            
            %% ******************************************************************************
            
            
            
            
            
            
            
            
            
            %% ***********Specific date for Triage *************************************
            
            if strcmp(SL_select2, 'Yes')
            
            
            if strcmp(T_date_select2, 'Yes')
            
            
            
                starting_time_triage=TXT_triage_date{i_geno};
            
            
            for i_trg=1:length(X_plug_time)
            
            if strcmp(X_plug_time(i_trg), '23:59:00') &  strcmp(X_plug_date(i_trg), starting_time_triage)
               start_day_triage_index(i_geno)=(i_trg);
                break;
            end
            end
                   try
                   Starting_time_triage_min=start_day_triage_index(i_geno)-((23*60+59)-(starting_hr*60));
                   catch
                       disp('Error ! check triage date')
                       set(handles.M_box,'String',[oldmsgs;{'Error ! check triage date'}] );drawnow
                return;
                       return
                   end
            
            
             %% Triage matrix
            
            
            try
            Triaged_matrix_new2=sum(X_fly(Starting_time_triage_min-(6*60):Starting_time_triage_min+(6*60), end-32+N_ch_start:end-32+N_ch_stop));
            
            
            catch
                try
              Triaged_matrix_new2=sum(X_fly(Starting_time_triage_min-(6*60):end, end-32+N_ch_start:end-32+N_ch_stop));
            
                catch
                     disp('ERROR! check triage date');
                     set(handles.M_box,'String',[oldmsgs;{'ERROR! check triage date'}] );drawnow
                     break
            
                end
            end
            
            
            
            
            
             try
            Triaged_matrix_cnd2=sum(X_fly( Starting_time_triage_min:Starting_time_triage_min+(24*60), end-32+N_ch_start:end-32+N_ch_stop));
                    catch
                       Triaged_matrix_cnd2=sum(X_fly(Starting_time_triage_min:end, end-32+N_ch_start:end-32+N_ch_stop));
                   end
            
            
            
            Triaged_matrix_new1=sum(X_fly(Starting_time_min_sleep:Ending_time_min_sleep,end-32+N_ch_start:end-32+N_ch_stop));
            
            
            
            
            
            
            
            
            
            %switch Triage_select
            
                if handles.SL_ch2.Value
            %% Triaged condition-2
            
            Triage_condn=['Triage for sleep (Last day activity between LON-6 to LON+6 < ' num2str(Triage_th2) ')' ];
            
            
            
            
            
            
            
            Triaged_matrix=Triaged_matrix_new2;
            
            for i_triaged=1:length(Triaged_matrix)
               if  Triaged_matrix_new2(i_triaged)<Triage_th2
                   Triaged_matrix(i_triaged)=0;
                   else
                   Triaged_matrix(i_triaged)=1;
               end
            end
            
                end
                %%
            
            
            
            
              if handles.SL_ch1.Value
            
            
            %% Triaged condition4 activi/30min <1
            
             Triage_condn=['Triage for sleep(Last day of activity <   ' num2str(Triage_th3), ')' ];
            
            
            
            
            
            
            
            
            
            Triaged_matrix= Triaged_matrix_cnd2;
            for i_triaged=1:length(Triaged_matrix_cnd2)
            
               if  ((Triaged_matrix_cnd2(i_triaged)<Triage_th3))
                   Triaged_matrix(i_triaged)=0;
                   else
                   Triaged_matrix(i_triaged)=1;
               end
            
            end
            
              end
            
            
              if handles.SL_ch3.Value
            %% Triaged condition-2
            
            Triage_condn=['Triage for sleep (Activity per 30 min < ' num2str(Triage_th4) ')' ];
            
            
            
            
            
            Triaged_matrix=Triaged_matrix_new1;
            
            for i_triaged=1:length(Triaged_matrix)
               if  (Triaged_matrix_new1(i_triaged)<((N_days_analysis*60*24*Triage_th4)/30))
                   Triaged_matrix(i_triaged)=0;
                   else
                   Triaged_matrix(i_triaged)=1;
               end
            end
            
            
              end
            
            %%
            if handles.SL_ch2.Value && handles.SL_ch1.Value
            
             Triage_condn=['Triage for sleep (Last day activity between LON-6 to LON+6 < ' num2str(Triage_th2) ' OR Last day of activity <   ' num2str(Triage_th3) ')' ];
            
            
             Triaged_matrix=Triaged_matrix_new2;
            
            for i_triaged=1:length(Triaged_matrix)
               if  Triaged_matrix_new2(i_triaged)<Triage_th2 | (Triaged_matrix_cnd2(i_triaged)<Triage_th3)
                   Triaged_matrix(i_triaged)=0;
                   else
                   Triaged_matrix(i_triaged)=1;
               end
            end
            
            
            end
            
            %%
            
            if handles.SL_ch2.Value & handles.SL_ch3.Value
            
            Triage_condn=['Triage for sleep (Last day activity between LON-6 to LON+6 < ' num2str(Triage_th2) ' OR activity per 30min <   ' num2str(Triage_th4) ')' ];
            
            
            Triaged_matrix=Triaged_matrix_new2;
            
            for i_triaged=1:length(Triaged_matrix)
               if  Triaged_matrix_new2(i_triaged)<Triage_th2 | (Triaged_matrix_new1(i_triaged)<((N_days_analysis*60*24*Triage_th4)/30))
                   Triaged_matrix(i_triaged)=0;
                   else
                   Triaged_matrix(i_triaged)=1;
               end
            end
            
            
            end
            
            
            %%
            
            if handles.SL_ch1.Value & handles.SL_ch3.Value
            
              Triage_condn=['Triage for sleep (Last day activity  < ' num2str(Triage_th3) ' OR activity per 30min <   ' num2str(Triage_th4) ')' ];
            
            
            Triaged_matrix=Triaged_matrix_cnd2;
            
            for i_triaged=1:length(Triaged_matrix)
               if  Triaged_matrix_cnd2(i_triaged)<Triage_th3 | (Triaged_matrix_new1(i_triaged)<((N_days_analysis*60*24*Triage_th4)/30))
                   Triaged_matrix(i_triaged)=0;
                   else
                   Triaged_matrix(i_triaged)=1;
               end
            end
            end
            
            
            
            %%
            if handles.SL_ch1.Value & handles.SL_ch2.Value & handles.SL_ch3.Value
            
              Triage_condn=['Triage for sleep (Last day activity  < ' num2str(Triage_th3) 'OR Last day activity between Lon-6 to Lon+6  < ' num2str(Triage_th2)  ' OR activity per 30min <   ' num2str(Triage_th4) ')' ];
            
            
            Triaged_matrix=Triaged_matrix_cnd2;
            
            for i_triaged=1:length(Triaged_matrix)
               if  Triaged_matrix_cnd2(i_triaged)<Triage_th3 | (Triaged_matrix_new1(i_triaged)<((N_days_analysis*60*24*Triage_th4)/30)) | Triaged_matrix_new2(i_triaged)<Triage_th2
            
                   Triaged_matrix(i_triaged)=0;
                   else
                   Triaged_matrix(i_triaged)=1;
               end
            end
            end
            
            
            
            else
            
            %Triage matrix
            
             Triaged_matrix_new1=sum(X_fly(Starting_time_min_sleep:Ending_time_min_sleep,end-32+N_ch_start:end-32+N_ch_stop));
            
            
            try
            Triaged_matrix_cnd2=sum(X_fly(Ending_time_triaged_sleep-(24*60):Ending_time_triaged_sleep, end-32+N_ch_start:end-32+N_ch_stop));
                    catch
                       Triaged_matrix_cnd2=sum(X_fly(end-(24*60):end, end-32+N_ch_start:end-32+N_ch_stop));
                   end
            
            
            
            try
            Triaged_matrix_new2=sum(X_fly(Ending_time_triaged_sleep-(6*60):Ending_time_triaged_sleep+(6*60), end-32+N_ch_start:end-32+N_ch_stop));
            
            catch
                try
             Triaged_matrix_new2=sum(X_fly(Ending_time_triaged_sleep-(6*60):end, end-32+N_ch_start:end-32+N_ch_stop));
            
                catch
                    disp('ERROR! change triage condition');
                set(handles.M_box,'String',[oldmsgs;{'ERROR! change triage condition'}] );drawnow
                end
            end
            
            
            
            
                  if handles.SL_ch2.Value
            %% Triaged condition-2
            
            Triage_condn=['Triage for sleep (Last day activity between LON-6 to LON+6 < ' num2str(Triage_th2) ')' ];
            
            
            
            Triaged_matrix=Triaged_matrix_new2;
            
            for i_triaged=1:length(Triaged_matrix)
               if  Triaged_matrix_new2(i_triaged)<Triage_th2
                   Triaged_matrix(i_triaged)=0;
                   else
                   Triaged_matrix(i_triaged)=1;
               end
            end
            
            
                  end
            
            
              if handles.SL_ch1.Value
            
            
            %% Triaged condition4 activi/30min <1
            
             Triage_condn=['Triage for sleep(Last day of activity <   ' num2str(Triage_th3), ')' ];
            %Triaged_matrix_cnd4=sum(Data_Eduction);
            
            
            %% Triaged condition-2
            
            
            Triaged_matrix= Triaged_matrix_cnd2;
            for i_triaged=1:length(Triaged_matrix_cnd2)
            
               if  ((Triaged_matrix_cnd2(i_triaged)<Triage_th3))
                   Triaged_matrix(i_triaged)=0;
                   else
                   Triaged_matrix(i_triaged)=1;
               end
            
            end
            
            
              end
            
            
            
            
            
              if handles.SL_ch3.Value
            %% Triaged condition-2
            
            Triage_condn=['Triage for sleep (Activity per 30 min < ' num2str(Triage_th4) ')' ];
            
            
            
            
            
            Triaged_matrix=Triaged_matrix_new1;
            
            for i_triaged=1:length(Triaged_matrix)
               if  (Triaged_matrix_new1(i_triaged)<((N_days_analysis*60*24*Triage_th4)/30))
                   Triaged_matrix(i_triaged)=0;
                   else
                   Triaged_matrix(i_triaged)=1;
               end
            end
            
            
              end
            
            
              %%
            if handles.SL_ch2.Value & handles.SL_ch1.Value
            
             Triage_condn=['Triage for sleep (Last day activity between LON-6 to LON+6 < ' num2str(Triage_th2) ' OR Last day of activity <   ' num2str(Triage_th3) ')' ];
            
            
             Triaged_matrix=Triaged_matrix_new2;
            
            for i_triaged=1:length(Triaged_matrix)
               if  Triaged_matrix_new2(i_triaged)<Triage_th2 | (Triaged_matrix_cnd2(i_triaged)<Triage_th3)
                   Triaged_matrix(i_triaged)=0;
                   else
                   Triaged_matrix(i_triaged)=1;
               end
            end
            
            
            end
            
            %%
            
            if handles.SL_ch2.Value & handles.SL_ch3.Value
            
            Triage_condn=['Triage for sleep (Last day activity between LON-6 to LON+6 < ' num2str(Triage_th2) ' OR activity per 30min <   ' num2str(Triage_th4) ')' ];
            
            
            Triaged_matrix=Triaged_matrix_new2;
            
            for i_triaged=1:length(Triaged_matrix)
               if  Triaged_matrix_new2(i_triaged)<Triage_th2 | (Triaged_matrix_new1(i_triaged)<((N_days_analysis*60*24*Triage_th4)/30))
                   Triaged_matrix(i_triaged)=0;
                   else
                   Triaged_matrix(i_triaged)=1;
               end
            end
            
            
            end
            
            
            %%
            
            if handles.SL_ch1.Value & handles.SL_ch3.Value
            
              Triage_condn=['Triage for sleep (Last day activity  < ' num2str(Triage_th3) ' OR activity per 30min <   ' num2str(Triage_th4) ')' ];
            
            
            Triaged_matrix=Triaged_matrix_cnd2;
            
            for i_triaged=1:length(Triaged_matrix)
               if  Triaged_matrix_cnd2(i_triaged)<Triage_th3 | (Triaged_matrix_new1(i_triaged)<((N_days_analysis*60*24*Triage_th4)/30))
                   Triaged_matrix(i_triaged)=0;
                   else
                   Triaged_matrix(i_triaged)=1;
               end
            end
            end
            
            
            
            %%
            if handles.SL_ch1.Value & handles.SL_ch2.Value & handles.SL_ch3.Value
            
              Triage_condn=['Triage for sleep (Last day activity  < ' num2str(Triage_th3) 'OR Last day activity between Lon-6 to Lon+6  < ' num2str(Triage_th2)  ' OR activity per 30min <   ' num2str(Triage_th4) ')' ];
            
            
            Triaged_matrix=Triaged_matrix_cnd2;
            
            for i_triaged=1:length(Triaged_matrix)
               if  Triaged_matrix_cnd2(i_triaged)<Triage_th3 | (Triaged_matrix_new1(i_triaged)<((N_days_analysis*60*24*Triage_th4)/30)) | Triaged_matrix_new2(i_triaged)<Triage_th2
            
                   Triaged_matrix(i_triaged)=0;
                   else
                   Triaged_matrix(i_triaged)=1;
               end
            end
            end
            
            
            
            
            end
            
            
            
            
            
            %% *************************************************************************
            
            
            N_Triaged=numel(Triaged_matrix)-nnz(Triaged_matrix);
            
            fprintf('Number of triaged flies for sleep analysis = %d \n', N_Triaged)
            
            
            
            
            
            if strcmp(SD_select2, 'Yes')
            
            
            SDL_t1=SD_Starting_time_min-((SD_interval-SD_duration)*60);
            SDL_t2=SD_ending_time_min;
            SDL_t3=SDL_t1-24*60;
            SDL_t4=SDL_t2-24*60;
            
            SDG_t1=SDL_t2;
            SDG_t2=SDG_t1+SG_interval*60-1;
            SDG_t3=SDG_t1-24*60*1;
            SDG_t4=SDG_t3+SG_interval*60-1;
            
            
            try
            SDL_data1=X_fly(SDL_t1:SDL_t2, end-32+N_ch_start:end-32+N_ch_stop);
            SDL_data2=X_fly(SDL_t3:SDL_t4, end-32+N_ch_start:end-32+N_ch_stop);
            
            SDG_data1=X_fly(SDG_t1:SDG_t2, end-32+N_ch_start:end-32+N_ch_stop);
            SDG_data2=X_fly(SDG_t3:SDG_t4, end-32+N_ch_start:end-32+N_ch_stop);
            catch
                disp('Error! check sleep deprivation entries')
                 set(handles.M_box,'String',[oldmsgs;{'Error! check sleep deprivation entries'}] );drawnow
                return
            end
            
            
            
            
            
             switch cum_fig_select2
            
                 case 'Yes'
            
            %figure;
            figure('NumberTitle', 'off', 'Name', geno_type);
            cum_sleep_loss(app, SDL_data1,SDL_data2,SDG_data1,SDG_data2,N_ch,geno_type,Triaged_matrix)
            
            %% Saving and closing cumulative sleep
            %try
            path = handles.path;
            myfolder = 'Cumulative_sleep' ;   % new folder name
            folder = mkdir([path,Project_name,filesep,myfolder]);
            path  = [path,Project_name,filesep,myfolder] ;
            
            
                    FigName   = geno_type;
                    temp1=[path,filesep,FigName,'.png'];
                            try
                    saveas(gca,temp1);
                    catch
                        disp("Error!! Please check genotype name. Don't use backslash '\' or '>'  for genotype naming if you need figures")
                                        set(handles.M_box,'String',[oldmsgs;{"Error!! Please check genotype name.  Don't use backslash '\' or '>' for genotype naming if you need figures"}] );drawnow
            
                                        return
                    end
            
              close(FigName);
            
            
            
            %end
             end
            
            
            %% ****************************************************************
            
            
            
            
            
            %% *********************************************************************
            
            
            SDL_sleep_matrix1=zeros(1, N_ch);
            SDL_sleep_matrix2=zeros(1, N_ch);
            
            SDG_sleep_matrix1=zeros(1, N_ch);
            SDG_sleep_matrix2=zeros(1, N_ch);
            
            
            
            % program for finding sleep latency after SD
            SDG_data1_Lat=X_fly(SDG_t1:end, end-32+N_ch_start:end-32+N_ch_stop);
            
            SD_Latency_matrix=zeros(1, N_ch);
            
             for Loff_colm=1:N_ch
                    SD_Count_Lat=0;
              if Triaged_matrix(Loff_colm)==0
                        SD_Latency_matrix(1,Loff_colm)=NaN;
              else
                    for Loff_row=1:length(SDG_data1_Lat)
                        if SDG_data1_Lat(Loff_row, Loff_colm)==0 & SD_Count_Lat<5
                            SD_Count_Lat=SD_Count_Lat+1;
                        elseif SD_Count_Lat==5
                                SD_Latency_matrix(1,Loff_colm)=Loff_row-7;
            
                                break;
                        else
                            SD_Count_Lat=0;
                        end
                    end
                    end
                end
            
            idx = SD_Latency_matrix(1,:) == -1; % Two condition give zero lattency
            SD_Latency_matrix(1,idx) = 0;
            SD_Latency=SD_Latency_matrix;
            
            
            
            
            
            
            
            
            for i_colm=1:N_ch % for first chanel of selection
                  N_count_SDL1=0;
                  N_count_SDL2=0;
            
                  N_count_SDG1=0;
                  N_count_SDG2=0;
            
                  if Triaged_matrix(i_colm)==0
                        SDL_sleep_matrix1(1, i_colm)=NaN;
                        SDL_sleep_matrix2(1, i_colm)=NaN;
            
                        SDG_sleep_matrix1(1, i_colm)=NaN;
                        SDG_sleep_matrix2(1, i_colm)=NaN;
                  else
                    for i_row_L=1:length(SDL_data1)
            
                 %****************************************************************************************
                        if SDL_data1(i_row_L,i_colm)==0
                            N_count_SDL1=N_count_SDL1+1;
                        elseif SDL_data1(i_row_L,i_colm)>0 & N_count_SDL1>=5 % condition for sleep detection
                             SDL_sleep_matrix1(1, i_colm)=SDL_sleep_matrix1(1,i_colm)+N_count_SDL1;
                             N_count_SDL1=0;
                        else
                            N_count_SDL1=0;
                        end
            
            
                 %***************************************************************************
                        if SDL_data2(i_row_L,i_colm)==0
                            N_count_SDL2=N_count_SDL2+1;
                        elseif SDL_data2(i_row_L,i_colm)>0 & N_count_SDL2>=5 % condition for sleep detection
                             SDL_sleep_matrix2(1, i_colm)=SDL_sleep_matrix2(1,i_colm)+N_count_SDL2;
                             N_count_SDL2=0;
                        else
                            N_count_SDL2=0;
                        end
                    end
                %***************************************************************************
            
                for i_row_G=1:length(SDG_data1)
                %****************************************************************************************
                        if SDG_data1(i_row_G,i_colm)==0
                            N_count_SDG1=N_count_SDG1+1;
                        elseif SDG_data1(i_row_G,i_colm)>0 & N_count_SDG1>=5 % condition for sleep detection
                             SDG_sleep_matrix1(1, i_colm)=SDG_sleep_matrix1(1,i_colm)+N_count_SDG1;
                             N_count_SDG1=0;
                        else
                            N_count_SDG1=0;
                        end
            
            
                 %***************************************************************************
                        if SDG_data2(i_row_G,i_colm)==0
                            N_count_SDG2=N_count_SDG2+1;
                        elseif SDG_data2(i_row_G,i_colm)>0 & N_count_SDG2>=5 % condition for sleep detection
                             SDG_sleep_matrix2(1, i_colm)=SDG_sleep_matrix2(1,i_colm)+N_count_SDG2;
                             N_count_SDG2=0;
                        else
                            N_count_SDG2=0;
                        end
                %***************************************************************************
            
            
                    end
            
                    %************************************************************************
                    if SDL_data1(i_row_L,i_colm)==0 & N_count_SDL1>=5 % condition for sleep if last element of chnel is zero
                            SDL_sleep_matrix1(1, i_colm)=SDL_sleep_matrix1(1,i_colm)+N_count_SDL1;
            
                    end
               %********************************************************************************
            
                    if SDL_data2(i_row_L,i_colm)==0 & N_count_SDL2>=5 % condition for sleep if last element of chnel is zero
                            SDL_sleep_matrix2(1, i_colm)=SDL_sleep_matrix2(1,i_colm)+N_count_SDL2;
            
                    end
              %*********************************************************************************
            
            
            
               %************************************************************************
                    if SDG_data1(i_row_G,i_colm)==0 & N_count_SDG1>=5 % condition for sleep if last element of chnel is zero
                            SDG_sleep_matrix1(1, i_colm)=SDG_sleep_matrix1(1,i_colm)+N_count_SDG1;
            
                    end
               %********************************************************************************
            
                    if SDG_data2(i_row_G,i_colm)==0 & N_count_SDG2>=5 % condition for sleep if last element of chnel is zero
                            SDG_sleep_matrix2(1, i_colm)=SDG_sleep_matrix2(1,i_colm)+N_count_SDG2;
            
                    end
              %*********************************************************************************
            
                end
              end
            
            
              sleep_lost1=SDL_sleep_matrix2- SDL_sleep_matrix1 ;
            
              sleep_gain1=SDG_sleep_matrix1- SDG_sleep_matrix2 ;
            
              sleep_lost_percentage1=(sleep_lost1./SDL_sleep_matrix2)*100;
              sleep_gain_percentage1=(sleep_gain1./SDG_sleep_matrix2)*100;
            
              for i_ch=1:N_ch
                  if sleep_lost_percentage1(i_ch)>=SD_percentage | strcmp(SD_Th_select2, 'No')
                      sleep_lost(i_ch)=sleep_lost1(i_ch);
                      sleep_gain(i_ch)=sleep_gain1(i_ch);
                      sleep_lost_percentage(i_ch)=sleep_lost_percentage1(i_ch);
                      sleep_gain_percentage(i_ch)=sleep_gain_percentage1(i_ch);
                      SD_Latency(i_ch)=SD_Latency(i_ch);
                  else
                      sleep_lost(i_ch)=NaN;
                      sleep_gain(i_ch)=NaN;
                      sleep_lost_percentage(i_ch)=NaN;
                      sleep_gain_percentage(i_ch)=NaN;
                      SD_Latency(i_ch)=NaN;
                  end
              end
            
            
            
            
             Number_SD_flys= numel(sleep_lost(~isnan(sleep_lost)));
            
             fprintf('Number of sleep deprived flies = %d \n', Number_SD_flys)
            
            
            
            
            % writing sleep deprivation results
            Column_headers1={'Sleep lost (min)', 'Percentage sleep lost', 'Sleep gain (min)', 'Percentage sleep gain', 'Latency after SD (min)'};
            Results1=[(sleep_lost(1:N_ch))' (sleep_lost_percentage(1:N_ch))' (sleep_gain(1:N_ch))' (sleep_gain_percentage(1:N_ch))', (SD_Latency(1:N_ch))'];
            
            
            %Fly_mean1=[nanmean(Results1);nanstd(Results1)];
            
            % Fly_mean1=nanmean(Results1,1);
            % Fly_std1=nanstd(Results1,0,1);
            
            end
            
            %% ******************************************************************************************
            
            %Data_all=Data_exp(Starting_time_min:Ending_time_min, :);
            
            t_latency_start=(starting_hr+((starting_min+1)/60));
            
            t_latency_array=t_latency_start:(1/60):(t_latency_start+length(Data_all)/60);
            
            Lat_matrix=[mod((t_latency_array(1:end-1))', 24) Data_all(:,end-32+N_ch_start:end-32+N_ch_stop)];% Matrix for calculating sleep latency
            
            ind_Loff=find(Lat_matrix(:,1)==t_Loff);
            ind_Lon=find(Lat_matrix(:,1)==t_Lon);
            
            
            
            
            sleep_matrix_24=zeros(1, N_ch);
            bout_matrix=zeros(1, N_ch);
            activity_matrix=zeros(1, N_ch);
            Bout_length_matrix=zeros(1, N_ch);
            
            
            
            Triaged=zeros(1, N_ch);
            
            
            seg_NZ=sum(Data_seg_all~=0,1);% find the ctivity
            total_activity_matrix=sum(Data_seg_all);% Total activity
            
            
            
            
              for i_colm=1:N_ch % for first chanel of selection
                  N_count_24=0;
                  if Triaged_matrix(i_colm)==0
                        Triaged(i_colm)=Triaged(i_colm)+1;
                        sleep_matrix_24(1,i_colm)=NaN;
                        bout_matrix(1,i_colm)=NaN;
                        activity_matrix(1,i_colm)=NaN;
                        Bout_length_matrix(1,i_colm)=NaN;
                        total_activity_matrix(1,i_colm)=NaN;
                  else
                    for i_row=1:length(Data_seg_all)
                        if Data_seg_all(i_row,i_colm)==0
                            N_count_24=N_count_24+1;
                        elseif Data_seg_all(i_row,i_colm)>0 & N_count_24>=5 % condition for sleep detection
                             sleep_matrix_24(1, i_colm)=sleep_matrix_24(1,i_colm)+N_count_24;
                             bout_matrix(1,i_colm)=bout_matrix(1,i_colm)+1;
                             Bout_length_matrix(1,i_colm)=Bout_length_matrix(1,i_colm)+N_count_24;
                            N_count_24=0;
                        else
                            N_count_24=0;
                        end
                    end
                    if Data_seg_all(i_row,i_colm)==0 & N_count_24>=5 % condition for sleep if last element of chnel is zero
                            sleep_matrix_24(1, i_colm)=sleep_matrix_24(1,i_colm)+N_count_24;
                            bout_matrix(1,i_colm)=bout_matrix(1,i_colm)+1;
                            Bout_length_matrix(1,i_colm)=Bout_length_matrix(1,i_colm)+N_count_24;
                    end
            
            
            
                end
              end
            
            
            sleep_per_24=sleep_matrix_24/N_days_analysis;% min and hrs measures
            Bout_number_per_24=bout_matrix/N_days_analysis;
            Avrage_bout_length_per_24=Bout_length_matrix./bout_matrix;
            Total_activity_per_24=total_activity_matrix/N_days_analysis;
            activity_counts_per_min_24=(total_activity_matrix./(60*24*N_days_analysis-sleep_matrix_24));
            
            
            %% Split L data and D data
            
            
             L_data=[];
             D_data=[];
             for i_n=1:N_days_analysis
                 L_start=1+(i_n-1)*(24)*60;
                 D_start=L_start+Photo_period*60;
                 L_data=[L_data;Data_seg_all(L_start:L_start+(Photo_period*60)-1,:)]; % L data
                 D_data=[D_data;Data_seg_all(D_start:(i_n*24*60),:)]; % D data
            
             end
            
            
             % for L data
            sleep_matrix_L=zeros(1, N_ch);
            bout_matrix_L=zeros(1, N_ch);
            activity_matrix_L=zeros(1, N_ch);
            Bout_length_matrix_L=zeros(1, N_ch);
            
            
            
            %Triaged=zeros(1, N_ch);
            
            
            seg_NZ_L=sum(L_data~=0,1);% find the ctivity
            total_activity_matrix_L=sum(L_data);% Total activity
            
            
            
            
              for i_L=1:N_ch % for first chanel of selection
                  N_count_L=0;
                  if Triaged_matrix(i_L)==0
                        %Triaged(i_L)=Triaged(i_L)+1;
                        sleep_matrix_L(1,i_L)=NaN;
                        bout_matrix_L(1,i_L)=NaN;
                        activity_matrix_L(1,i_L)=NaN;
                        Bout_length_matrix_L(1,i_L)=NaN;
                        total_activity_matrix_L(1,i_L)=NaN;
                  else
                    for i_row_L=1:length(L_data)
                        if L_data(i_row_L,i_L)==0
                            N_count_L=N_count_L+1;
                        elseif L_data(i_row_L,i_L)>0 & N_count_L>=5 % condition for sleep detection
                             sleep_matrix_L(1, i_L)=sleep_matrix_L(1,i_L)+N_count_L;
                             bout_matrix_L(1,i_L)=bout_matrix_L(1,i_L)+1;
                             Bout_length_matrix_L(1,i_L)=Bout_length_matrix_L(1,i_L)+N_count_L;
                            N_count_L=0;
                        else
                            N_count_L=0;
                        end
                    end
                    if L_data(i_row_L,i_L)==0 & N_count_L>=5 % condition for sleep if last element of chnel is zero
                            sleep_matrix_L(1, i_L)=sleep_matrix_L(1,i_L)+N_count_L;
                            bout_matrix_L(1,i_L)=bout_matrix_L(1,i_L)+1;
                            Bout_length_matrix_L(1,i_L)=Bout_length_matrix_L(1,i_L)+N_count_L;
                    end
            
            
            
                end
              end
            
            
            sleep_per_L=sleep_matrix_L/N_days_analysis;% min and hrs measures
            Bout_number_per_L=bout_matrix_L/N_days_analysis;
            Avrage_bout_length_per_L=Bout_length_matrix_L./bout_matrix_L;
            Total_activity_per_L=total_activity_matrix_L/N_days_analysis;
            activity_counts_per_min_L=(total_activity_matrix_L./(60*24*N_days_analysis-sleep_matrix_L));
            
            
            
            
            
            % For D data
            
            sleep_matrix_D=zeros(1, N_ch);
            bout_matrix_D=zeros(1, N_ch);
            activity_matrix_D=zeros(1, N_ch);
            Bout_length_matrix_D=zeros(1, N_ch);
            
            
            
            %Triaged=zeros(1, N_ch);
            
            
            seg_NZ_D=sum(D_data~=0,1);% find the ctivity
            total_activity_matrix_D=sum(D_data);% Total activity
            
            
            
            
              for i_D=1:N_ch % for first chanel of selection
                  N_count_D=0;
                  if Triaged_matrix(i_D)==0
                       % Triaged(i_D)=Triaged(i_D)+1;
                        sleep_matrix_D(1,i_D)=NaN;
                        bout_matrix_D(1,i_D)=NaN;
                        activity_matrix_D(1,i_D)=NaN;
                        Bout_length_matrix_D(1,i_D)=NaN;
                        total_activity_matrix_D(1,i_D)=NaN;
                  else
                    for i_row_D=1:length(D_data)
                        if D_data(i_row_D,i_D)==0
                            N_count_D=N_count_D+1;
                        elseif D_data(i_row_D,i_D)>0 & N_count_D>=5 % condition for sleep detection
                             sleep_matrix_D(1, i_D)=sleep_matrix_D(1,i_D)+N_count_D;
                             bout_matrix_D(1,i_D)=bout_matrix_D(1,i_D)+1;
                             Bout_length_matrix_D(1,i_D)=Bout_length_matrix_D(1,i_D)+N_count_D;
                            N_count_D=0;
                        else
                            N_count_D=0;
                        end
                    end
                    if D_data(i_row_D,i_D)==0 & N_count_D>=5 % condition for sleep if last element of chnel is zero
                            sleep_matrix_D(1, i_D)=sleep_matrix_D(1,i_D)+N_count_D;
                            bout_matrix_D(1,i_D)=bout_matrix_D(1,i_D)+1;
                            Bout_length_matrix_D(1,i_D)=Bout_length_matrix_D(1,i_D)+N_count_D;
                    end
            
            
            
                end
              end
            
            
            
            sleep_per_D=sleep_matrix_D/N_days_analysis;% min and hrs measures
            Bout_number_per_D=bout_matrix_D/N_days_analysis;
            Avrage_bout_length_per_D=Bout_length_matrix_D./bout_matrix_D;
            Total_activity_per_D=total_activity_matrix_D/N_days_analysis;
            activity_counts_per_min_D=(total_activity_matrix_D./(60*24*N_days_analysis-sleep_matrix_D));
            
            
            
            %% **********************************************
            
            
            
            %% program for finding sleep latency
            
            
            Latency_matrix=zeros(length(ind_Loff), N_ch);
            
            
            
            %% ******************************************=
            
            for i_off=1:length(ind_Loff)
                Lat_matrix_Loff=Lat_matrix(ind_Loff(i_off):end,:);
                for Loff_colm=1:N_ch
                    Count_Lat=0;
                    Latency_min=0;
            
                    if Triaged_matrix(Loff_colm)==0
                    Latency_matrix(i_off,Loff_colm)=NaN;
            
                    else
            
                         for Loff_row=1:length(Lat_matrix_Loff)
                        Latency_min=Latency_min+1;
            
                        if Lat_matrix_Loff(Loff_row, Loff_colm+1)==0 & Count_Lat<5
                            Count_Lat=Count_Lat+1;
                        elseif Count_Lat==5
            
                                Latency_matrix(i_off,Loff_colm)=Latency_min-7;
                                                    if Latency_matrix(i_off,Loff_colm)<0
                                                       Latency_matrix(i_off,Loff_colm)=0;
                                                    end
            
                                break;
                        else
                            Count_Lat=0;
                        end
                    end
                    end
                end
            end
            
            Latency=Latency_matrix;
            Average_Latency=nanmean(Latency_matrix,1);
            
            
            end
            
            
            
            
            
            EAP_matrix=t_Loff-5:0.5:t_Loff-1;
            MAP_matrix=t_Lon-5:0.5:t_Lon-1;
            %% ***********************************************
                    %Latency_matrix=[];
                   E_anticipation_matrix=[];
                   E_anticipation_phase_matrix=[];
                   M_anticipation_matrix=[];
                   M_anticipation_phase_matrix=[];
                   %EA_matrix=[];
                   %MA_matrix=[];
                   M_anticipation_phase=[];
                   E_anticipation_phase=[];
            
                   A_MA_ph1=[];
            
            
            
            
            %%      Eduction plot
            
            if strcmp(Ed_select2, 'Yes')
            
            try
            Data_Eduction=X_fly_Raw.data(Starting_time_Eduction_min:Ending_time_Eduction_min, end-32+N_ch_start:end-32+N_ch_stop);
            catch
                disp('ERROR! starting day or number of days of eduction/anticipation is out of range');
                set(handles.M_box,'String',[oldmsgs;{'ERROR! starting day or number of days of eduction/anticipation is out of range'}] );drawnow
                return;
            end
            
            
            
            
               % Eduction triage conditon
            
               %% ***********Specific date for Triage *************************************
            
            
            if strcmp(T_date_ed_select2, 'Yes')
            
            
            
                starting_time_triage=TXT_triage_date{i_geno};
            
            
            for i_trg=1:length(X_plug_time)
            
            if strcmp(X_plug_time(i_trg), '23:59:00') &  strcmp(X_plug_date(i_trg), starting_time_triage)
               start_day_triage_index(i_geno)=(i_trg);
                break;
            end
            end
            
                   Starting_time_triage_min=start_day_triage_index(i_geno)-((23*60+59)-(starting_hr*60));
            
            
            try
            Triaged_matrix_ed_new2=sum(X_fly(Starting_time_triage_min-(6*60):Starting_time_triage_min+(6*60), end-32+N_ch_start:end-32+N_ch_stop));
            
            catch
                try
              Triaged_matrix_ed_new2=sum(X_fly(Starting_time_triage_min-(6*60):end, end-32+N_ch_start:end-32+N_ch_stop));
            
                catch
                     disp('ERROR! check triage date');
                     set(handles.M_box,'String',[oldmsgs;{'ERROR! check triage date'}] );drawnow
                     break
            
                end
            end
            
            
            
                     try
            Triaged_matrix_ed_cnd2=sum(X_fly( Starting_time_triage_min:Starting_time_triage_min+(24*60), end-32+N_ch_start:end-32+N_ch_stop));
                    catch
                       Triaged_matrix_ed_cnd2=sum(X_fly(Starting_time_triage_min:end, end-32+N_ch_start:end-32+N_ch_stop));
                   end
            
            
             Triaged_matrix_new1=sum(X_fly_Raw.data(Starting_time_Eduction_min+(((24-Photo_period)/2)*60):Ending_time_Eduction_min+(((24-Photo_period)/2)*60), end-32+N_ch_start:end-32+N_ch_stop));
            
            
            
            
            
            
            %switch Triage_ed_select
            
                 if handles.ed_ch2.Value
            %% Triaged condition-2
            
            Triage_condn_ed=['Triage for eduction/anticipation (Last day activity between LON-6 to LON+6 < ' num2str(Triage_th2_ed) ')' ];
            
            
            
            
            Triaged_matrix_ed=Triaged_matrix_ed_new2;
            
            for i_triaged=1:length(Triaged_matrix_ed)
               if  Triaged_matrix_ed_new2(i_triaged)<Triage_th2_ed
                   Triaged_matrix_ed(i_triaged)=0;
                   else
                   Triaged_matrix_ed(i_triaged)=1;
               end
            end
            
                 end
            
            
            
            
             if handles.ed_ch1.Value
            
            
            %% Triaged condition4 activi/30min <1
            
            Triage_condn_ed=['Triage for eduction/anticipation(Last day of activity <   ' num2str(Triage_th3_ed), ')' ];
            
            
            % Triaged_matrix_ed_cnd4=sum(X_fly_Raw.data(Starting_time_Eduction_min+(((24-Photo_period)/2)*60):Ending_time_Eduction_min+(((24-Photo_period)/2)*60), end-32+N_ch_start:end-32+N_ch_stop));
            
            
            
            
            
             % Threshold for activity/30 min < Th4 : it shoud be total min/(60*24*Th4*48);  48--> 24*2 (30mins)
             %Threshod4_mod_ed= round(((Starting_time_triage_min+(24*60))-Starting_time_min)/(60*24))*Triage_th4_ed*48;
            
            
            Triaged_matrix_ed= Triaged_matrix_ed_cnd2;
            for i_triaged=1:length(Triaged_matrix_ed_cnd2)
            
               if  ((Triaged_matrix_ed_cnd2(i_triaged)<Triage_th3_ed))
                   Triaged_matrix_ed(i_triaged)=0;
                   else
                   Triaged_matrix_ed(i_triaged)=1;
               end
            
            end
            
             end
            
             if handles.ed_ch3.Value
            %% Triaged condition-2
            
            Triage_condn_ed=['Triage for eduction/anticipation (Activity per 30 min < ' num2str(Triage_th4_ed) ')' ];
            
            
            
            
            
            Triaged_matrix_ed=Triaged_matrix_new1;
            
            for i_triaged=1:length(Triaged_matrix_ed)
               if  (Triaged_matrix_new1(i_triaged)<((N_days_Eduction*60*24*Triage_th4_ed)/30))
                   Triaged_matrix_ed(i_triaged)=0;
                   else
                   Triaged_matrix_ed(i_triaged)=1;
               end
            end
            
            
             end
            
            
            %%
            
            if handles.ed_ch2.Value & handles.ed_ch1.Value
            
                Triage_condn_ed=['Triage for eduction/anticipation (Last day activity between LON-6 to LON+6 < ' num2str(Triage_th2_ed) 'OR Last day of activity <   ' num2str(Triage_th3_ed) ')' ];
            
            
            Triaged_matrix_ed=Triaged_matrix_ed_new2;
            
            for i_triaged=1:length(Triaged_matrix_ed)
               if  Triaged_matrix_ed_new2(i_triaged)<Triage_th2_ed | ((Triaged_matrix_ed_cnd2(i_triaged)<Triage_th3_ed))
                   Triaged_matrix_ed(i_triaged)=0;
                   else
                   Triaged_matrix_ed(i_triaged)=1;
               end
            end
            
            
            end
            
            %%
            
            if handles.ed_ch2.Value & handles.ed_ch3.Value
            
                Triage_condn_ed=['Triage for eduction/anticipation (Last day activity between LON-6 to LON+6 < ' num2str(Triage_th2_ed) 'OR Activity per 30 min <   ' num2str(Triage_th4_ed) ')' ];
            
            
            Triaged_matrix_ed=Triaged_matrix_ed_new2;
            
            for i_triaged=1:length(Triaged_matrix_ed)
               if  Triaged_matrix_ed_new2(i_triaged)<Triage_th2_ed | (Triaged_matrix_new1(i_triaged)<((N_days_Eduction*60*24*Triage_th4_ed)/30))
                   Triaged_matrix_ed(i_triaged)=0;
                   else
                   Triaged_matrix_ed(i_triaged)=1;
               end
            end
            
            
            end
            
            %%
            
            if handles.ed_ch1.Value & handles.ed_ch3.Value
            
                Triage_condn_ed=['Triage for eduction/anticipation (Last day activity  < ' num2str(Triage_th3_ed) 'OR Activity per 30 min <   ' num2str(Triage_th4_ed) ')' ];
            
            
            Triaged_matrix_ed=Triaged_matrix_ed_cnd2;
            
            for i_triaged=1:length(Triaged_matrix_ed)
               if  ((Triaged_matrix_ed_cnd2(i_triaged)<Triage_th3_ed)) | (Triaged_matrix_new1(i_triaged)<((N_days_Eduction*60*24*Triage_th4_ed)/30))
                   Triaged_matrix_ed(i_triaged)=0;
                   else
                   Triaged_matrix_ed(i_triaged)=1;
               end
            end
            end
            
            
            
            if handles.ed_ch1.Value & handles.ed_ch3.Value
            
                Triage_condn_ed=['Triage for eduction/anticipation (Last day activity  < ' num2str(Triage_th3_ed) 'OR Activity per 30 min <   ' num2str(Triage_th4_ed) ')' ];
            
            
            Triaged_matrix_ed=Triaged_matrix_ed_cnd2;
            
            for i_triaged=1:length(Triaged_matrix_ed)
               if  ((Triaged_matrix_ed_cnd2(i_triaged)<Triage_th3_ed)) | (Triaged_matrix_new1(i_triaged)<((N_days_Eduction*60*24*Triage_th4_ed)/30)) | (Triaged_matrix_ed_new2(i_triaged)<Triage_th2_ed)
                   Triaged_matrix_ed(i_triaged)=0;
                   else
                   Triaged_matrix_ed(i_triaged)=1;
               end
            end
            end
            
            
            
            
            else
            
            %switch Triage_ed_select
            
            
            
            Triaged_matrix_new1=sum(X_fly_Raw.data(Starting_time_Eduction_min+(((24-Photo_period)/2)*60):Ending_time_Eduction_min+(((24-Photo_period)/2)*60), end-32+N_ch_start:end-32+N_ch_stop));
            
            Triaged_matrix_ed_cnd2=sum(X_fly(Ending_time_Eduction_min+(((24-Photo_period)/2)*60)-(24*60):Ending_time_Eduction_min+(((24-Photo_period)/2)*60), end-32+N_ch_start:end-32+N_ch_stop));
            
            try
            Triaged_matrix_ed_new2=sum(X_fly((Ending_time_Eduction_min+((24-Photo_period)/2)*60)-(6*60):(Ending_time_Eduction_min+((24-Photo_period)/2)*60)+(6*60), end-32+N_ch_start:end-32+N_ch_stop));
            
            catch
                try
             Triaged_matrix_ed_new2=sum(X_fly((Ending_time_Eduction_min+((24-Photo_period)/2)*60)-(6*60):end, end-32+N_ch_start:end-32+N_ch_stop));
            
                catch
                    disp('ERROR! change triage condition');
                    set(handles.M_box,'String',[oldmsgs;{'ERROR! change triage condition'}] );drawnow
            
                end
            end
            
            
            
              if handles.ed_ch2.Value
            %% Triaged condition-2
            
            Triage_condn_ed=['Triage for eduction/anticipation (Last day activity between LON-6 to LON+6 < ' num2str(Triage_th2_ed) ')' ];
            
            
             Triaged_matrix_ed=Triaged_matrix_ed_new2;
            
            for i_triaged=1:length(Triaged_matrix_ed)
               if   Triaged_matrix_ed_new2(i_triaged)<Triage_th2_ed
                   Triaged_matrix_ed(i_triaged)=0;
               else
                  Triaged_matrix_ed(i_triaged)=1;
               end
            end
            
            
            
              end
            
             if handles.ed_ch1.Value
            
            
            %% Triaged condition4 activi/30min <1
            
             Triage_condn_ed=['Triage for eduction/anticipation (Last day of activity <   ' num2str(Triage_th3_ed), ')' ];
            
            
            
            Triaged_matrix_ed= Triaged_matrix_ed_cnd2;
            for i_triaged=1:length(Triaged_matrix_ed_cnd2)
            
               if  (Triaged_matrix_ed_cnd2(i_triaged)<Triage_th3_ed)
                   Triaged_matrix_ed(i_triaged)=0;
                    else
                  Triaged_matrix_ed(i_triaged)=1;
               end
            
            end
             end
            
            
             if handles.ed_ch3.Value
            %% Triaged condition-2
            
            Triage_condn_ed=['Triage for eduction/anticipation (Activity per 30 min < ' num2str(Triage_th4_ed) ')' ];
            
            
            
            
            
            
            Triaged_matrix_ed=Triaged_matrix_new1;
            
            for i_triaged=1:length(Triaged_matrix_ed)
               if  (Triaged_matrix_new1(i_triaged)<((N_days_Eduction*60*24*Triage_th4_ed)/30))
                   Triaged_matrix_ed(i_triaged)=0;
                   else
                   Triaged_matrix_ed(i_triaged)=1;
               end
            end
             end
            
            
            
             %%
            
            if handles.ed_ch2.Value & handles.ed_ch1.Value
            
                Triage_condn_ed=['Triage for eduction/anticipation (Last day activity between LON-6 to LON+6 < ' num2str(Triage_th2_ed) 'OR Last day of activity <   ' num2str(Triage_th3_ed) ')' ];
            
            
            Triaged_matrix_ed=Triaged_matrix_ed_new2;
            
            for i_triaged=1:length(Triaged_matrix_ed)
               if  Triaged_matrix_ed_new2(i_triaged)<Triage_th2_ed | ((Triaged_matrix_ed_cnd2(i_triaged)<Triage_th3_ed))
                   Triaged_matrix_ed(i_triaged)=0;
                   else
                   Triaged_matrix_ed(i_triaged)=1;
               end
            end
            
            
            end
            
            %%
            
            if handles.ed_ch2.Value & handles.ed_ch3.Value
            
                Triage_condn_ed=['Triage for eduction/anticipation (Last day activity between LON-6 to LON+6 < ' num2str(Triage_th2_ed) 'OR Activity per 30 min <   ' num2str(Triage_th4_ed) ')' ];
            
            
            Triaged_matrix_ed=Triaged_matrix_ed_new2;
            
            for i_triaged=1:length(Triaged_matrix_ed)
               if  Triaged_matrix_ed_new2(i_triaged)<Triage_th2_ed | (Triaged_matrix_new1(i_triaged)<((N_days_Eduction*60*24*Triage_th4_ed)/30))
                   Triaged_matrix_ed(i_triaged)=0;
                   else
                   Triaged_matrix_ed(i_triaged)=1;
               end
            end
            
            
            end
            
            %%
            
            if handles.ed_ch1.Value & handles.ed_ch3.Value
            
                Triage_condn_ed=['Triage for eduction/anticipation (Last day activity  < ' num2str(Triage_th3_ed) 'OR Activity per 30 min <   ' num2str(Triage_th4_ed) ')' ];
            
            
            Triaged_matrix_ed=Triaged_matrix_ed_cnd2;
            
            for i_triaged=1:length(Triaged_matrix_ed)
               if  ((Triaged_matrix_ed_cnd2(i_triaged)<Triage_th3_ed)) | (Triaged_matrix_new1(i_triaged)<((N_days_Eduction*60*24*Triage_th4_ed)/30))
                   Triaged_matrix_ed(i_triaged)=0;
                   else
                   Triaged_matrix_ed(i_triaged)=1;
               end
            end
            end
            
            
            
            if handles.ed_ch1.Value & handles.ed_ch3.Value
            
                Triage_condn_ed=['Triage for eduction/anticipation (Last day activity  < ' num2str(Triage_th3_ed) 'OR Activity per 30 min <   ' num2str(Triage_th4_ed) ')' ];
            
            
            Triaged_matrix_ed=Triaged_matrix_ed_cnd2;
            
            for i_triaged=1:length(Triaged_matrix_ed)
               if  ((Triaged_matrix_ed_cnd2(i_triaged)<Triage_th3_ed)) | (Triaged_matrix_new1(i_triaged)<((N_days_Eduction*60*24*Triage_th4_ed)/30)) | (Triaged_matrix_ed_new2(i_triaged)<Triage_th2_ed)
                   Triaged_matrix_ed(i_triaged)=0;
                   else
                   Triaged_matrix_ed(i_triaged)=1;
               end
            end
            end
            
            
            
            
            end
            
            
            
            %% ***********************************************
            
            
            N_Triaged_ed=numel(Triaged_matrix_ed)-nnz(Triaged_matrix_ed);
            
            fprintf('Number of triaged flies for eduction/anticipation analysis = %d \n', N_Triaged_ed)
            
            
            
            if strcmp(Anticipation_select2, 'Yes')
            
            %% ***************************************************************
            % program for finding anticipation
            
            
             Data_all_anti=X_fly_Raw.data(Starting_time_anti-(7*60):Ending_time_anti-(7*60), :);
            %t_anti_start=(1+(1/60));
            t_anti_start=(mod((t_Lon-7)+24,24)+(1/60));
            t_anti_array=t_anti_start:(1/60):(t_anti_start+length(Data_all_anti)/60);
            
            %
            
            
            Lat_matrix_anti=[mod((t_anti_array(1:end-1))', 24) Data_all_anti(:,end-32+N_ch_start:end-32+N_ch_stop)];% Matrix for calculating sleep latency
            
            
            % ind_Loff_anti=find(Lat_matrix_anti(:,1)==t_Loff);
            % ind_Lon_anti=find(Lat_matrix_anti(:,1)==t_Lon);
            
            
            ind_Loff_anti=find(abs(Lat_matrix_anti(:,1)-t_Loff)<1e-2);
            ind_Lon_anti=find(abs(Lat_matrix_anti(:,1)-t_Lon)<1e-2);
            
            %% Modified
            
            
            %% Evening anticipation
            
            %% ******************************************=
             EA_matrix_all=zeros(12, N_ch);
            for i_off=1:length(ind_Loff_anti)
                EA_matrix=[];
                    for Loff_colm=1:N_ch
            
                    if Triaged_matrix_ed(Loff_colm)==0
            
                    E_anticipation_matrix(i_off,Loff_colm)=NaN;
                    EA_matrix(:, Loff_colm)=NaN;
            
                    else
                        EE_anti=(sum(Lat_matrix_anti(((ind_Loff_anti(i_off))-(3*60)):ind_Loff_anti(i_off),Loff_colm+1))./sum(Lat_matrix_anti(((ind_Loff_anti(i_off))-(6*60)):ind_Loff_anti(i_off),Loff_colm+1)))-0.5;
                         EE_anti(isinf(EE_anti)|isnan(EE_anti)) = 0;
                        E_anticipation_matrix(i_off,Loff_colm)=EE_anti;
            
            
                   for i_E=1:12
                            EA_matrix(i_E, Loff_colm)=sum(Lat_matrix_anti(((ind_Loff_anti(i_off))-(6-(i_E-1)*0.5)*60)+1:((ind_Loff_anti(i_off))-(6-(i_E)*0.5)*60),Loff_colm+1));
                   end
            
            
             end
                    end
            
                    EA_matrix_all=EA_matrix_all+EA_matrix;  % Sum of 30 min bin for Ndays
            end
            
            %% *************************************************
            
            
             MA_matrix_all=zeros(12, N_ch);
            
            for i_on=1:length(ind_Lon_anti)
                MA_matrix=[];
              for Loff_colm=1:N_ch
            
            
                    if Triaged_matrix_ed(Loff_colm)==0
                    M_anticipation_matrix(i_on,Loff_colm)=NaN;
                    MA_matrix(:, Loff_colm)=NaN;
                    else
            
                   MM_anti=(sum(Lat_matrix_anti((ind_Lon_anti(i_on)-(3*60)):ind_Lon_anti(i_on),Loff_colm+1))./sum(Lat_matrix_anti((ind_Lon_anti(i_on)-(6*60)):ind_Lon_anti(i_on),Loff_colm+1)))-0.5;
                   MM_anti(isinf(MM_anti)|isnan(MM_anti)) = 0;
                   M_anticipation_matrix(i_on,Loff_colm)=MM_anti;
            
                  % 30 min bin from 2AM to 8 AM
                   for i_M=1:12
                            MA_matrix(i_M, Loff_colm)=sum(Lat_matrix_anti(((ind_Lon_anti(i_on))-(6-(i_M-1)*0.5)*60)+1:((ind_Lon_anti(i_on))-(6-(i_M)*0.5)*60),Loff_colm+1));
                   end
            
                    end
            
              end
             MA_matrix_all=MA_matrix_all+MA_matrix;  % Sum of 30 min bin for Ndays
            end
            
            
            %***********************************************************************
            
            MA_Dif=[];
            EA_Dif=[];
            
            for ii_M=1:9
                MA_Dif(ii_M,:)=MA_matrix_all(4+(ii_M-1),:)-MA_matrix_all(ii_M,:);% difference (3.30-4)- (2-2.30)
                EA_Dif(ii_M,:)=EA_matrix_all(4+(ii_M-1),:)-EA_matrix_all(ii_M,:);
            end
            
            
                         [MAP_max,M_ind]=max(MA_Dif); % find maximum difference
                          M_anticipation_phase_matrix=MAP_matrix(M_ind); % Midpoint
            
            
                          [EAP_max,E_ind]=max(EA_Dif);
                          E_anticipation_phase_matrix=EAP_matrix(E_ind);
            
            
            M_anticipation=nanmean(M_anticipation_matrix,1);% morning anticipation
            E_anticipation=nanmean(E_anticipation_matrix,1);% Evening anticipation
            
            
            
            % phase claculution converted to CR
            for i_ant=1:N_ch
            
                if Triaged_matrix_ed(i_ant)==0
                    M_anticipation_phase(i_ant)=NaN;
                else
                    M_anticipation_phase(i_ant)=mod((nanmean(M_anticipation_phase_matrix(:,i_ant))-t_Lon)+24,24);
                end
            
                   if Triaged_matrix_ed(i_ant)==0
                            E_anticipation_phase(i_ant)=NaN;
                    else
                          E_anticipation_phase(i_ant)=mod((nanmean(E_anticipation_phase_matrix(:,i_ant))-t_Lon)+24,24);
                end
            end
            
            
            end
            
            
            
            
            
            
            %% Data normalization, Data/(average over all data)
            
            % step1 avoide the triaged flies
            
            [nr_ed,nc_ed]=size(Data_Eduction);
            
            Triaged_ed=zeros(1, N_ch);
            
            for i_ed=1:nc_ed
                if Triaged_matrix_ed(i_ed)==0
                   Data_Eduction(:,i_ed)=NaN;
                   Triaged_ed(i_ed)=Triaged_ed(i_ed)+1;
                else
                Data_Eduction(:,i_ed)=Data_Eduction(:,i_ed);
                end
            end
            
            
            if strcmp(Eduction_select2, 'Yes')
            % step 2 average over N flies
            
            %Data_Eduction_avg_Nflies=nanmean(Data_Eduction,2);
            
            Data_Eduction_avg_Nflies=nansum(Data_Eduction,2);
            
            %%  30 min bins
            
            N_bin=length(Data_Eduction_avg_Nflies)/30;
            
            Edution_bin=[];
            for i_bin=1:N_bin
                Edution_bin(i_bin)=sum(Data_Eduction_avg_Nflies((i_bin-1)*30+1:(i_bin*30)));
            end
            
            Edution_bin_for_Ndays=Edution_bin';% 30 min for N days
            
            %Normalized_Edution_bin_for_Ndays=Edution_bin_for_Ndays/nanmean(Edution_bin_for_Ndays);% 30 min bin for N days
            
            
            
            
            
            
            % 30 bin average over N days
            
            Eduction_all=zeros(48,1);
            
            
            
            for ig=1:N_days_Eduction;
                Eduction_all=Eduction_all+Edution_bin_for_Ndays((ig-1)*48+1:(ig)*48);% 30 bin for 24 hr
            end
            
            %% ***********************************************************************
            avg_bin_N_days=Eduction_all/N_days_Eduction;% Normalized value average over N days
            %% ***********************************************************************
            
            
            
            ed_result=[avg_bin_N_days];% 30 min bins for 24 hr (avg over N days)
            %err_result=[err_result Normalized_err];
            
            %end
            ed_result_Ndays=[Edution_bin_for_Ndays];%30 min bin for N days
            
            
            
            %if strcmp(Ind_Ed_select2, 'Yes')
            
            Data_Eduction_Nflies=Data_Eduction;
            %%  30 min bins
             Edution_bin_Nflies=[];
            for i_bin=1:N_bin
                Edution_bin_Nflies(i_bin,:)=sum(Data_Eduction_Nflies((i_bin-1)*30+1:(i_bin*30),:));
            end
            
              Edution_bin_for_N_flies_Ndays=Edution_bin_Nflies';% 30 min for N days
            
            %end
            
            
            
            
            end
            
            end
            
            
            
            
            
            
            
            %% ************************************************************************
            
            
            if strcmp(Interval_select2, 'Yes')
            
            
            [nr, nc]=size(Data_seg_all);
            
            t_array_all_inval=t_start:(1/60):(t_start+length(Data_seg_all)/60);
            
            sleep_matrix_inval=zeros(floor(nr/N_inval_min), N_ch);
            bout_matrix_inval=zeros(floor(nr/N_inval_min), N_ch);
            %activity_matrix_inval=zeros(floor(nr/N_inval_min), N_ch);
            Bout_length_matrix_inval=zeros(floor(nr/N_inval_min), N_ch);
            Total_activity_matrix_inval=zeros(floor(nr/N_inval_min), N_ch);
            
            
            for i_seg=1:floor(nr/N_inval_min)
                data_seg=Data_seg_all(((i_seg-1)*N_inval_min+1):i_seg*N_inval_min,:);%get the data with in the interval
                seg_NZ=sum(data_seg~=0,1);% find the ctivity
                total_activity=sum(data_seg);% Total activity per interval
                for i_colm=1:N_ch % for first chanel of selection
                    N_count=0;
            
                    if Triaged_matrix(i_colm)== 0
            
            
                sleep_matrix_inval(:, i_colm)=NaN;
                bout_matrix_inval(:, i_colm)=NaN;
                %activity_matrix_inval(:, i_colm)=NaN;
                Bout_length_matrix_inval(:, i_colm)=NaN;
                Total_activity_matrix_inval(:, i_colm)=NaN;
            
                    else
                    for i_row=1:length(data_seg)
                       try
                        if data_seg(i_row,i_colm)==0
                            N_count=N_count+1;
                        elseif data_seg(i_row,i_colm)>0 & N_count>=5 % condition for sleep detection
                            sleep_matrix_inval(i_seg, i_colm)=sleep_matrix_inval(i_seg,i_colm)+N_count;% for sleep
            
                            if (i_row-N_count==1 & Data_seg_all((i_seg-1)*N_inval_min, i_colm)~=0) | i_row-N_count~=1
                                bout_matrix_inval(i_seg, i_colm)=bout_matrix_inval(i_seg,i_colm)+1;
                                Bout_length_matrix_inval(i_seg, i_colm)=Bout_length_matrix_inval(i_seg, i_colm)+N_count; % Bout length
            
                            end
            
                            N_count=0;
            
                           % Check zeros in the first rows and is it a continuation from previous segment
                       elseif (data_seg(i_row,i_colm)>0 & N_count<5) & Data_seg_all((((i_seg-1)*N_inval_min+i_row)-5):(((i_seg-1)*N_inval_min+i_row)-1),i_colm)==0 % condition for sleep detection
                            sleep_matrix_inval(i_seg, i_colm)=sleep_matrix_inval(i_seg,i_colm)+N_count;% for bout length
                            N_count=0;
            
                        else
                            N_count=0;
                        end
                       catch
                           N_count=0;
                   end
                    end
                    try
                         % Check the zeros are in last entries OR zeros are initiate in this segment and continues in next segment
                    if (data_seg(i_row,i_colm)==0 & N_count>=5)  | (N_count>0 & N_count<5 & Data_seg_all((((i_seg-1)*N_inval_min+i_row)-(N_count-1)):(((i_seg-1)*N_inval_min+i_row)+(5-N_count)),i_colm)==0)  % condition for sleep if last element of chnel is zero
                            sleep_matrix_inval(i_seg, i_colm)=sleep_matrix_inval(i_seg,i_colm)+N_count;
            
                            if i_row-N_count~=0 % Check whether all entries are zero
                            bout_matrix_inval(i_seg, i_colm)=bout_matrix_inval(i_seg,i_colm)+1;% for bout number
                            x=Data_seg_all((i_seg*N_inval_min-(N_count-1)):end,i_colm);
                            idx = find(x~=0, 1, 'first');
                            Bout_length_matrix_inval(i_seg, i_colm)=Bout_length_matrix_inval(i_seg, i_colm)+(idx-1);% Bout length
                            end
            
            
                    end
                    catch
                    end
                    % activity_matrix_inval(i_seg, i_colm)=total_activity(i_colm)/seg_NZ(i_colm);%Average Activty per interval
                     Total_activity_matrix_inval(i_seg, i_colm)=Total_activity_matrix_inval(i_seg,i_colm)+total_activity(i_colm);% for sleep
                end
                end
            end
            
            Avg_bout_ln=Bout_length_matrix_inval./bout_matrix_inval;% average bout length
            
            t_array_inval=t_start:N_inval_hr:((floor(nr/N_inval_min))*N_inval_hr)+t_start; % for timing interval
            
            Total_activity_per_inval=[Total_activity_matrix_inval];
            
            sleep_inval=sleep_matrix_inval';
            end
            
            %% *******************************************************************
            
            
            
            
            
            if strcmp(SL_select2, 'Yes')
            
            Column_headers={'Genotype', 'Run number','Monitor number','Channel',Triage_condn, 'Sleep/Day (min)', 'Bout number/Day','Average bout length/Day', 'Total activity/Day', 'Activity/waking min',	'Latency (min)'};
            
            
            
            try
                Column_headers_inval={'Genotype', 'Run number','Monitor number','Channel', Triage_condn};
            end
            
            
            
            
            Results=[Monitor_1 row_header Triaged' sleep_per_24' Bout_number_per_24' Avrage_bout_length_per_24' Total_activity_per_24'...
                activity_counts_per_min_24' Average_Latency'];
            
            
            end
            
            
            % anticipation results
            
            try
                Results_anti=[Monitor_1 row_header Triaged_ed' (M_anticipation)' (E_anticipation)' (M_anticipation_phase)' (E_anticipation_phase)'];
            end
            
            try
                Results_ed_TW=[Monitor_1 row_header Triaged_ed'];
            end
            
            
            head_mean={Run_number, '','Mean'};
            head_std={Run_number, '','SD', ''};
            head_SE={Run_number, '','SE', ''};
            
            %% Individual eduction
            
            try
                Column_headers_inval_ed={'Genotype', 'Run number','Monitor number','Channel', Triage_condn_ed};
            end
            
            
            try
                Result_Ind_Ed=[Monitor_1 row_header Triaged_ed' Edution_bin_for_N_flies_Ndays];
                R_ed_Ind=[R_ed_Ind;GT Run_1 num2cell(Result_Ind_Ed)];
            end
            
            
            
            %% Last activity date and time
            
            if strcmp(SL_select2, 'Yes')
            
            % try
            % R_all=[R_all;GT Run_1 num2cell(Results) num2cell(Results1) L_NZ num2cell(L_NZ_day')];
            % catch
            % R_all=[R_all;GT Run_1 num2cell(Results) L_NZ num2cell(L_NZ_day')]    ;
            % end
            
            
            try
            R_all=[R_all;GT Run_1 num2cell(Results) num2cell(Results1)];
            catch
            R_all=[R_all;GT Run_1 num2cell(Results)]    ;
            end
            
            
            %try
            
               Results_L_D=[Monitor_1 row_header Triaged' sleep_per_L' sleep_per_D' Bout_number_per_L' Bout_number_per_D'...
                   Avrage_bout_length_per_L' Avrage_bout_length_per_D' Total_activity_per_L' Total_activity_per_D' activity_counts_per_min_L'...
                   activity_counts_per_min_D'];
                R_all_L_D=[R_all_L_D;GT Run_1 num2cell(Results_L_D)]    ;
            %end
            end
            
            
            try
            R_all_anti=[R_all_anti;GT Run_1 num2cell(Results_anti)];
            end
            
            
            try
            R_all_ed_TW=[R_all_ed_TW;GT Run_1 num2cell(Results_ed_TW)];
            end
            
            
            
            
            
            
            
            
            
            try
            
            Sleep_int=[Monitor_1 row_header Triaged' sleep_inval];
            Avg_bout_ln_int=[Monitor_1 row_header Triaged' Avg_bout_ln'];
            Bout_number_int=[Monitor_1 row_header Triaged' bout_matrix_inval'];
            Total_activity_int=[Monitor_1 row_header Triaged' Total_activity_per_inval'];
            
            
            Sleep_interval1=[Sleep_interval1;GT Run_1 num2cell(Sleep_int)];
            
            BoutN_interval1=[BoutN_interval1;GT Run_1 num2cell(Bout_number_int)];
            
            BoutL_interval1=[BoutL_interval1;GT Run_1 num2cell(Avg_bout_ln_int)];
            
            Activity_interval1=[Activity_interval1; GT Run_1 num2cell(Total_activity_int)];
            
            end
            
            
            
            
            
            try
                active_ed=nnz(~(Triaged_ed));
                R_ed=[R_ed; GT{1} Run_1{1} num2cell(Monitor_1(1)) num2cell(active_ed) num2cell(ed_result')];
                R_ed_Ndays=[R_ed_Ndays; GT{1} Run_1{1} num2cell(Monitor_1(1)) num2cell(ed_result_Ndays')];
               % R_ed_Day_by_Day=[R_ed_Day_by_Day; GT{1} num2cell(Run_1(1)) num2cell(Monitor_1(1)) num2cell(Day_by_Day_Edution_bin')];
            
            end
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            %% Periodogram
            
            
            
            switch periodo_select2
                case 'Yes'
            %% Periodogram data
            %periodo_data=Data_seg_all;
            
            try
            periodo_data=Data_seg_periodo(start_time_periodo_min:end_time_periodo_min,:);
            catch
                disp('Error! start day OR number of days of periodogram is out of range')
                set(handles.M_box,'String',[oldmsgs;{'Error! start day OR number of days of periodogram is out of range'}] );drawnow
                return
            end
            
            if strcmp(T_date_per_select2, 'Yes')
            
            
            
             starting_time_triage=TXT_triage_date{i_geno};
            
            
            for i_trg=1:length(X_plug_time)
            
            if strcmp(X_plug_time(i_trg), '23:59:00') &  strcmp(X_plug_date(i_trg), starting_time_triage)
               start_day_triage_index(i_geno)=(i_trg);
                break;
            end
            end
            
                   Starting_time_triage_min=start_day_triage_index(i_geno)-((23*60+59)-(starting_hr*60));
            
            
             %% Triaged matrix
              try
            Triaged_matrix_per_new2=sum(X_fly(Starting_time_triage_min-(6*60):Starting_time_triage_min+(6*60), end-32+N_ch_start:end-32+N_ch_stop));
            
             catch
            
            
            
                try
              Triaged_matrix_per_new2=sum(X_fly(Starting_time_triage_min-(6*60):end, end-32+N_ch_start:end-32+N_ch_stop));
            
                catch
                     disp('ERROR! check triage date');
                     set(handles.M_box,'String',[oldmsgs;{'ERROR! check triage date'}] );drawnow
                     break
            
                end
            end
            
            
            
            
            
            
            
            
            try
            Triaged_matrix_per_cnd2=sum(X_fly( Starting_time_triage_min:Starting_time_triage_min+(24*60), end-32+N_ch_start:end-32+N_ch_stop));
                    catch
                       Triaged_matrix_per_cnd2=sum(X_fly(Starting_time_triage_min:end, end-32+N_ch_start:end-32+N_ch_stop));
                   end
            
            
            
            Triaged_matrix_per_new1=sum(Data_seg_periodo(start_time_periodo_min:end_time_periodo_min,:));
            
            
             %%
            
            
            %% Triage testing
            
              if handles.per_ch2.Value
            %% Triaged condition-2
            
            Triage_condn_per=['Triage for periodogram (Last day activity between LON-6 to LON+6 < ' num2str(Triage_th2_per) ')' ];
            
            Triaged_matrix_per=Triaged_matrix_per_new2;
            
            for i_triaged=1:length(Triaged_matrix_per)
               if Triaged_matrix_per_new2(i_triaged)<Triage_th2_per
                   Triaged_matrix_per(i_triaged)=0;
               else
                   Triaged_matrix_per(i_triaged)=1;
               end
            end
              end
            
            
             %%
            
            if handles.per_ch1.Value
            
            Triage_condn_per=['Triage for periodogram(Last day of activity <   ' num2str(Triage_th3_per), ')' ];
            
            Triaged_matrix_per=Triaged_matrix_per_cnd2;
            for i_triaged=1:length(Triaged_matrix_per_cnd2)
            
               if (Triaged_matrix_per_cnd2(i_triaged)<Triage_th3_per)
                   Triaged_matrix_per(i_triaged)=0;
                   else
                   Triaged_matrix_per(i_triaged)=1;
               end
            
            end
            
            end
            %%
            
            
            if handles.per_ch3.Value
            %% Triaged condition-2
            
            Triage_condn_per=['Triage for periodogram (Activity per 30 min < ' num2str(Triage_th4_per) ')' ];
            Triaged_matrix_per=Triaged_matrix_per_new1;
            
            for i_triaged=1:length(Triaged_matrix_per)
                  if  (Triaged_matrix_per_new1(i_triaged)<((N_day_per*60*24*Triage_th4_per)/30))
                   Triaged_matrix_per(i_triaged)=0;
                   else
                   Triaged_matrix_per(i_triaged)=1;
               end
            end
            end
            
            
            %%
            
                 if (handles.per_ch2.Value & handles.per_ch1.Value)
            %% Triaged condition-2
            
            Triage_condn_per=['Triage for periodogram (Last day activity between LON-6 to LON+6 < ' num2str(Triage_th2_per)  'OR Last day of activity <   ' num2str(Triage_th3_per), ')' ];
            
            Triaged_matrix_per=Triaged_matrix_per_new2;
            
            for i_triaged=1:length(Triaged_matrix_per)
               if  Triaged_matrix_per_new2(i_triaged)<Triage_th2_per| (Triaged_matrix_per_cnd2(i_triaged)<Triage_th3_per)
                   Triaged_matrix_per(i_triaged)=0;
               else
                   Triaged_matrix_per(i_triaged)=1;
               end
            end
                 end
            
            
            
            %%
            
             if (handles.per_ch2.Value & handles.per_ch3.Value)
            %% Triaged condition-2
            
            Triage_condn_per=['Triage for periodogram (Last day activity between LON-6 to LON+6 < ' num2str(Triage_th2_per)  ' OR Activity per 30 min <   ' num2str(Triage_th4_per), ')' ];
            
            Triaged_matrix_per=Triaged_matrix_per_new2;
            
            for i_triaged=1:length(Triaged_matrix_per)
               if  Triaged_matrix_per_new2(i_triaged)<Triage_th2_per| Triaged_matrix_per_new1(i_triaged)<(N_day_per*60*24*Triage_th4_per)/30
                   Triaged_matrix_per(i_triaged)=0;
               else
                   Triaged_matrix_per(i_triaged)=1;
               end
            end
            
                      end
            
            
             %%
               if (handles.per_ch1.Value & handles.per_ch3.Value)
            %% Triaged condition-2
            
            Triage_condn_per=['Triage for periodogram (Last day activity  < ' num2str(Triage_th3_per)  ' OR Activity per 30 min <   ' num2str(Triage_th4_per), ')' ];
            
            
            Triaged_matrix_per=Triaged_matrix_per_new2;
            for i_triaged=1:length(Triaged_matrix_per)
               if  Triaged_matrix_per_new1(i_triaged)<(N_day_per*60*24*Triage_th4_per)/30 | (Triaged_matrix_per_cnd2(i_triaged)<Triage_th3_per)
                   Triaged_matrix_per(i_triaged)=0;
               else
                   Triaged_matrix_per(i_triaged)=1;
               end
            end
            
                     end
            
            
            %%
              if (handles.per_ch1.Value & handles.per_ch2.Value & handles.per_ch3.Value)
            %% Triaged condition-2
            
            Triage_condn_per=['Triage for periodogram (Last day activity  < ' num2str(Triage_th3_per) ' OR Last day of activity between Lon-6 and Lon+6 <   ' num2str(Triage_th4_per)  ' OR Activity per 30 min <   ' num2str(Triage_th4_per), ')' ];
            
            Triaged_matrix_per=Triaged_matrix_per_cnd2;
            
            for i_triaged=1:length(Triaged_matrix_per)
               if  Triaged_matrix_per_new2(i_triaged)<Triage_th2_per| (Triaged_matrix_per_cnd2(i_triaged)<Triage_th3_per) | Triaged_matrix_per_new1(i_triaged)<(N_day_per*60*24*Triage_th4_per)/30
                   Triaged_matrix_per(i_triaged)=0;
               else
                   Triaged_matrix_per(i_triaged)=1;
               end
            end
            
              end
            
            %%
            
            else
            
            
            %% ***********************************************************************************
             try
            Triaged_matrix_per_new2=sum(Data_seg_periodo(end_time_periodo_min-(6*60):end_time_periodo_min+(6*60), :));
            
            catch
                try
             Triaged_matrix_per_new2=sum(Data_seg_periodo(end_time_periodo_min-(6*60):end, :));
            
                catch
                    disp('ERROR! change triage condition');
                    set(handles.M_box,'String',[oldmsgs;{'ERROR! change triage condition'}] );drawnow
            
                end
            end
            
             Triaged_matrix_per_cnd2=sum(Data_seg_periodo(end_time_periodo_min-(24*60):end_time_periodo_min, :));
            
            Triaged_matrix_per_cnd4=sum(Data_seg_periodo(start_time_periodo_min:end_time_periodo_min,:));
            
            
            %% **********************************************************************************************
            if handles.per_ch2.Value
            %% Triaged condition-2
            
            
            
            Triage_condn_per=['Triage for periodogram (Last day activity between LON-6 to LON+6 < ' num2str(Triage_th2_per) ')' ];
            
            Triaged_matrix_per=Triaged_matrix_per_new2;
            
            for i_triaged=1:length(Triaged_matrix_per)
               if  Triaged_matrix_per_new2(i_triaged)<Triage_th2_per
                   Triaged_matrix_per(i_triaged)=0;
               else
                   Triaged_matrix_per(i_triaged)=1;
               end
            end
            
                 end
            
            
            if handles.per_ch1.Value
            
            
            %% Triaged condition4 activi/30min <1
            
             Triage_condn_per=['Triage for periodogram (Last day of activity <   ' num2str(Triage_th3_per), ')' ];
            
            Triaged_matrix_per= Triaged_matrix_per_cnd2;
            for i_triaged=1:length(Triaged_matrix_per_cnd2)
            
               if  (Triaged_matrix_per_cnd2(i_triaged)<Triage_th3_per)
                   Triaged_matrix_per(i_triaged)=0;
                   else
                   Triaged_matrix_per(i_triaged)=1;
               end
            
            end
            
            end
            
            
            if handles.per_ch3.Value
            %% Triaged condition-2
            
            
            Triage_condn_per=['Triage for periodogram (Activity per 30 min < ' num2str(Triage_th4_per) ')' ];
            
            
            Triaged_matrix_per= Triaged_matrix_per_cnd4;
            for i_triaged=1:length(Triaged_matrix_per)
            
               if  Triaged_matrix_per_cnd4(i_triaged)<(N_day_per*60*24*Triage_th4_per)/30
                   Triaged_matrix_per(i_triaged)=0;
                   else
                   Triaged_matrix_per(i_triaged)=1;
               end
            
            end
            end
            
            
                 if (handles.per_ch2.Value & handles.per_ch1.Value)
            %% Triaged condition-2
            
            Triage_condn_per=['Triage for periodogram (Last day activity between LON-6 to LON+6 < ' num2str(Triage_th2_per)  'OR Last day of activity <   ' num2str(Triage_th3_per), ')' ];
            
            Triaged_matrix_per=Triaged_matrix_per_new2;
            
            for i_triaged=1:length(Triaged_matrix_per)
               if  Triaged_matrix_per_new2(i_triaged)<Triage_th2_per| (Triaged_matrix_per_cnd2(i_triaged)<Triage_th3_per)
                   Triaged_matrix_per(i_triaged)=0;
               else
                   Triaged_matrix_per(i_triaged)=1;
               end
            end
                 end
            
            
            
            
            
             if (handles.per_ch2.Value & handles.per_ch3.Value)
            %% Triaged condition-2
            
            Triage_condn_per=['Triage for periodogram (Last day activity between LON-6 to LON+6 < ' num2str(Triage_th2_per)  ' OR Activity per 30 min <   ' num2str(Triage_th4_per), ')' ];
            
            Triaged_matrix_per=Triaged_matrix_per_new2;
            
            for i_triaged=1:length(Triaged_matrix_per)
               if  Triaged_matrix_per_new2(i_triaged)<Triage_th2_per| Triaged_matrix_per_cnd4(i_triaged)<(N_day_per*60*24*Triage_th4_per)/30
                   Triaged_matrix_per(i_triaged)=0;
               else
                   Triaged_matrix_per(i_triaged)=1;
               end
            end
            
                      end
            
            
            
               if (handles.per_ch1.Value & handles.per_ch3.Value)
            %% Triaged condition-2
            
            Triage_condn_per=['Triage for periodogram (Last day activity  < ' num2str(Triage_th3_per)  ' OR Activity per 30 min <   ' num2str(Triage_th4_per), ')' ];
            
            
            Triaged_matrix_per=Triaged_matrix_per_cnd4;
            for i_triaged=1:length(Triaged_matrix_per)
               if  Triaged_matrix_per_cnd4(i_triaged)<(N_day_per*60*24*Triage_th4_per)/30 | (Triaged_matrix_per_cnd2(i_triaged)<Triage_th3_per)
                   Triaged_matrix_per(i_triaged)=0;
               else
                   Triaged_matrix_per(i_triaged)=1;
               end
            end
            
                     end
            
            
            
              if (handles.per_ch1.Value & handles.per_ch2.Value & handles.per_ch3.Value)
            %% Triaged condition-2
            
            Triage_condn_per=['Triage for periodogram (Last day activity  < ' num2str(Triage_th3_per) ' OR Last day of activity between Lon-6 and Lon+6 <   ' num2str(Triage_th4_per)  ' OR Activity per 30 min <   ' num2str(Triage_th4_per), ')' ];
            
            Triaged_matrix_per=Triaged_matrix_per_cnd2;
            
            for i_triaged=1:length(Triaged_matrix_per)
               if  Triaged_matrix_per_new2(i_triaged)<Triage_th2_per| (Triaged_matrix_per_cnd2(i_triaged)<Triage_th3_per) | Triaged_matrix_per_cnd4(i_triaged)<(N_day_per*60*24*Triage_th4_per)/30
                   Triaged_matrix_per(i_triaged)=0;
               else
                   Triaged_matrix_per(i_triaged)=1;
               end
            end
            
              end
            
            
            
            
            
            
            end
            
            %% ***********************************************
            
            
            N_Triaged_per=numel(Triaged_matrix_per)-nnz(Triaged_matrix_per);
            
            fprintf('Number of triaged flies for periodogram = %d \n', N_Triaged_per)
            
            
            
            
            
            
            
            
            periodo_data_GT=[];
            
            
            [np_r,np_c]=size(periodo_data);
            
            header_period=repmat({geno_type},1,np_c);
            
            periodo_data_GT=[header_period; num2cell(periodo_data)];
            
            periodo_data_all= [periodo_data_all periodo_data_GT];
            
            %*******************Period info for individual channel************************************
            
            
             period_ind_Ch=zeros(N_ch,6);
             Triaged_per=zeros(1, N_ch);
            for i_colm=1:N_ch
            
                    if Triaged_matrix_per(i_colm)== 0
                        period_ind_Ch(i_colm,:)=NaN;
                        Triaged_per(i_colm)=Triaged_per(i_colm)+1;
            
                    else
            
            fign=1;
            noscip=5;
            
            y_new=periodo_data(:,i_colm);
            N_bin=round(length(y_new)/30);
            ydata=[];
            
            for i_bin=1:N_bin
                ydata(i_bin)=sum(y_new((i_bin-1)*30+1:(i_bin*30)));% for 30 min bin 0.5hr resulution
            end
            
            xdata=1:30:length(y_new)-1;
            Time_data=xdata;
            Activity_data=ydata;
            variance=var(Activity_data);
            
            
            % min_per=20*60; % Minimum period to evaluate
            % max_per=28*60; % maximum period to evaluate
            %
            % sig_vector=20:0.5:28;
            
            
            
            min_per=per_min*60; % Minimum period to evaluate
            max_per=per_max*60; % maximum period to evaluate
            
            sig_vector=per_min:0.5:per_max;
            
            
            alpha_val=1-Significance_level;
            sig_line= chi2inv(alpha_val,sig_vector*2); % for significant line with alpha=0.01
            
            per_rslt=zeros(1,noscip*3); % Final Results
            
            minpoints=length(find(Time_data<=min_per));
            maxpoints=length(find(Time_data<=max_per));
            
            xaxis1=(minpoints:maxpoints)*(Time_data(2)-Time_data(1));
            xaxis=xaxis1/60; % convert minute into hrs
            %sig_line= chi2inv(0.99,sig_vector*2);
            
            
            len_Qp=maxpoints-minpoints+1;
            Qp=zeros(1,len_Qp);
            
            
             npoints=minpoints:maxpoints;
             for k=1:length(npoints)
             nblocks=floor(length(Time_data)/npoints(k));
             table1=zeros(nblocks,npoints(k));
             for i=1:nblocks
             for i_LNZ=1:npoints(k)
             pos=i_LNZ+npoints(k)*(i-1);
             table1(i,i_LNZ)=Activity_data(pos);
             end
             end
             if (len_Qp==1)
             Qp=npoints(k)*nblocks*var(mean(table1,1))/variance;
             else
             Qp(k)=npoints(k)*nblocks*var(mean(table1,1))/variance;
             end
             end
            
            
            
             Dqp1= Qp'-sig_line';
             P_0=xaxis';
            for iqp=1:length(Dqp1)
             if Dqp1(iqp)<0
                Dqp1(iqp)=0;
                P_0(iqp)=NaN;
             end
            end
            
            
            P_10=xaxis';;
            for iqq=1:length(Dqp1)
             if Dqp1(iqq)<=10
                P_10(iqq)=NaN;
             end
            end
            
            
            P_30=xaxis';
            for iqq=1:length(Dqp1)
             if Dqp1(iqq)<=30
                P_30(iqq)=NaN;
             end
            end
            
            
            
            
            
            
            periodo_result=[Qp' sig_line' Dqp1 P_0 P_10 P_30];
            
             if max(periodo_result(:,3))==0
                period_ind_Ch(i_colm,:)=[NaN NaN Dqp1(i_colm) NaN NaN NaN];
             else
            P_period=find(periodo_result(:,3)==max(periodo_result(:,3)));
            
            try
             period_ind_Ch(i_colm,:)=periodo_result(P_period,:);
            catch
                 period_ind_Ch(i_colm,:)=[NaN NaN Dqp1(i_colm) NaN NaN NaN];
            end
            
             end
            
                    end
            end
            
            R_ind_periodo1=[Monitor_1 row_header Triaged_per' period_ind_Ch];
            R_ind_periodo=[R_ind_periodo;GT Run_1 num2cell(R_ind_periodo1)];
            
            
            end
            
            
            
            %% Actogram
            
            switch Acto_select2
            
                case 'Yes'
                    try
                        Data_seg_actogram=X_fly_Raw.data(Starting_time_actogram_min:Ending_time_actogram_min,end-32+N_ch_start:end-32+N_ch_stop);
                    catch
            %            disp('Error! start day OR number of days of actogram is out of range')
            %     set(handles.M_box,'String',[oldmsgs;{'Error! start day OR number of days of actogram is out of range'}] );drawnow
            %     return
            Data_seg_actogram=X_fly_Raw.data(Starting_time_actogram_min:end,end-32+N_ch_start:end-32+N_ch_stop);
                    end
            
                    path = handles.path;   % mention your path
            myfolder = [path,Project_name,filesep,'Actogram'] ;   % new folder name
            
            if ~exist(myfolder, 'dir')
            
            mkdir(myfolder);
            end
            
            try
            for i_nch=1:N_ch
            
                 sleepmat_actogram(app, Data_seg_actogram, i_nch,N_days_actogram, geno_type,t_Lon,Monitor_number,row_header(i_nch),myfolder)
            end
            catch
              disp('Warning! Check actogram input')
            end
            
            end
            
            
            
            
            
            
            disp('*****************************************************')
            
            end
            
            
            %% ************************************************************************************* %%
            
                               %% Environment Monitor Chart ploting Begins%%
            
            %% ************************************************************************************* %%
            
            
            try
            
            myfolder = 'Environment_condition_plot' ;   % new folder name
            folder = mkdir([handles.path,Project_name,filesep,myfolder]);% create the new folder
            path_Environment  = [handles.path,Project_name,filesep,myfolder] ; % set this as path
            
            
            for i_EV=1:length(Environment_Cell)
            
                Run_number= Environment_Cell{i_EV,1};
                Environment_file= Environment_Cell{i_EV,2};
            
            Environment_monitor_name=['Monitor' num2str(Environment_file) '.txt'];
            EV_path=strcat(handles.path, Environment_monitor_name);
            T_Environment=readtable(EV_path);
            T_Environment_cleaned=table2array(T_Environment(:, end-31:end));
            
            
            FigName   = ['Environment_chart_',Run_number];
            figure('NumberTitle', 'off', 'Name',FigName);
            
            
            Time_Environment=linspace(0,length(T_Environment_cleaned)/(60*24),length(T_Environment_cleaned));
            
            title(['Environmental conditions Run-',Run_number])
            yyaxis left
            plot(Time_Environment,T_Environment_cleaned(:,4));
            ylabel('Average light over min (lux)')
            
            yyaxis right
            plot(Time_Environment,T_Environment_cleaned(:,9)/10);
            ylabel('Average temperature over min (degC)')
            xlabel('Days')
            
            
            
                   temp_EV=[path_Environment,filesep,FigName,'.png'];
                   saveas(gca,temp_EV);
                   close(FigName);
            
            
            end
            
            end
            
            %% **********************   Environment Monitoring ENDS **************************************************************************
            
            try
            [~,~,X_per_ind]=unique(R_ind_periodo(:,1));
            C_per_ind = accumarray(X_per_ind,1:size(R_ind_periodo,1),[],@(r){R_ind_periodo(r,:)});
            
            
            [nr_per,nc_per]=size(C_per_ind );
            
            Head_periodo={Triage_condn_per,'P','S', 'P-S', 'Period for P-S>0', 'Period for P-S>=10' 'Period for P-S>=30'};
            End_sp_per=NaN*ones(1,11);
            Estar_per=cell(1,11);
            Estar_per(:)={'******'};
            
            for i_per_ind=1:nr_per
            
            err_per_ind=C_per_ind{i_per_ind};
            [~, order_per_ind] = sort(cell2mat(err_per_ind(:, 1)));
            sortederr_per_ind =  sortrows(err_per_ind, [1 2 3 4]);% sort in the order of column 1,2,3,4
            N_active=nnz(~(cell2mat(sortederr_per_ind(:, 5))));
            Title_N=[Title5 ' ' num2str(N_active)];
            
            per_ind_mean2=nanmean(cell2mat(sortederr_per_ind(:, 6:end)),1);
                per_ind_std2=nanstd(cell2mat(sortederr_per_ind(:, 6:end)),0,1);
            
                N_for_SE=sum(~isnan(cell2mat(sortederr_per_ind(:, 6:end))),1); % find number of Values(not NaN in each column)
                per_ind_SE2= per_ind_std2./sqrt(N_for_SE);
            
            
            
            Column_headers_period={'Genotype', 'Run number','Monitor number','Channel'};
            
            Periodo_all=[Periodo_all;Column_headers_period Head_periodo;sortederr_per_ind;sortederr_per_ind(1,1) head_mean Title_N num2cell(per_ind_mean2);sortederr_per_ind(1,1) head_std num2cell(per_ind_std2);sortederr_per_ind(1,1) head_SE num2cell(per_ind_SE2);num2cell(End_sp_per);Estar_per;num2cell(End_sp_per)];
            Periodo_all_SC=[Periodo_all_SC;sortederr_per_ind];
            
            % Screening data for periodogram
            screen_data_period=[screen_data_period;sortederr_per_ind(1,1) num2cell(per_ind_mean2(3:4))];
            
            end
            
               screen_data_period_T=cell2table(screen_data_period,...
                 "VariableNames", ["GT" "P_S" "Period_P_S"]);
            
            screen_data_period_T_P_S = sortrows(screen_data_period_T, 'P_S');
            screen_data_period_T_Period = sortrows(screen_data_period_T, 'Period_P_S');
            
            screen_data_period_headers_New={"Genotype" "P-S" "Period(P-S>0)"};
            screen_data_period=[screen_data_period_headers_New;screen_data_period];
            
            end
            
            
            try
            % ################### P value calculation of periodogram screening begins**************************************************************************
            
            
            Periodo_all_SC=Periodo_all_SC(:,1:9); % If you need P-s>10 and P-S>20 add all
            
             Periodo_all_SC=cell2table(Periodo_all_SC,...
                  "VariableNames", ["Genotype" "Run" "Monitor" "channel" "Triage" "P" "S" "P_S" "Period_P_S"]);
            
            % Get data for the reference genotype
            ref_data_per = Periodo_all_SC(contains(Periodo_all_SC.Genotype, ref_genotype), :);
            % Extract unique genotypes
            genotypes_per = unique(Periodo_all_SC.Genotype);
            
            
            % Initialize result storage
            results_per = table('Size', [numel(genotypes_per), 3], 'VariableTypes', {'string', 'double', 'double'}, ...
                'VariableNames', {'Genotype', 'P_S_Pval', 'Period_Pval'});
            
            
            % Loop over each genotype to compare with the reference genotype
            for i = 1:numel(genotypes_per)
                % Get data for the current genotype
                current_data_per = Periodo_all_SC(strcmp(Periodo_all_SC.Genotype, genotypes_per{i}), :);
            
                    isPresent = any(ismember(ref_data_per.Genotype, genotypes_per{i}));
            
            if isPresent
              P_S_Pval = 1;
              Period_Pval = 1;
            
            else
                % Perform t-tests
                [~,P_S_Pval] = ttest2(ref_data_per.P_S, current_data_per.P_S);
                [~, Period_Pval] = ttest2(ref_data_per.Period_P_S, current_data_per.Period_P_S);
            
            end
                % Store results
                results_per.Genotype(i) = genotypes_per{i};
                results_per.P_S_Pval(i) = P_S_Pval;
                results_per.Period_Pval(i) = Period_Pval;
            end
            
            
                Pval_data_all_per=results_per;
            
            
            % ################### P value calculation of periodogram screening  Ends*************************************************************************
            
            end
            
            
            
            
            
            
            
            
            
            
            %% *****************************************************************************************
             % Combining the same genotype together
            
            if strcmp(SL_select2, 'Yes')
            
            [~,~,X]=unique(R_all(:,1));
            C = accumarray(X,1:size(R_all,1),[],@(r){R_all(r,:)});
            
            
            [nr,nc]=size(C);
            
            Ln_gp=[];
            length(C{1});
             for i_split=1:nr
             [nr_gp, nc_gp]=size(C{i_split});
             Ln_gp(i_split)=nr_gp;
             end
            
            Ln_gp_max=max(Ln_gp);
            
            Individual_head_1=NaN*ones(Ln_gp_max+2,1);
            Individual_head_mean={'Mean'};
            Individual_head_SD={'SD'};
            Individual_head_SE={'SE'};
            
            Individual_head =[num2cell(Individual_head_1);Individual_head_mean;Individual_head_SD;Individual_head_SE];
            
            
            
            for igp=1:nr
            
            err=C{igp};
            
            
            [~, order] = sort(cell2mat(err(:, 1)));
            sortederr =  sortrows(err, [1 2 3 4]);% sort in the order of column 1 and 2
            
            
            %     Fly_mean2=nanmean(cell2mat(sortederr(:, 6:end-3)),1);
            %     Fly_std2=nanstd(cell2mat(sortederr(:, 6:end-3)),0,1);
            
                  Fly_mean2=nanmean(cell2mat(sortederr(:, 6:end)),1);
                Fly_std2=nanstd(cell2mat(sortederr(:, 6:end)),0,1);
            
                N_active_array(igp)=nnz(~(cell2mat(sortederr(:, 5))));
                Title_N=[Title5 ' ' num2str(N_active_array(igp))];
            
                Fly_SE2= Fly_std2/sqrt(N_active_array(igp));
            
            
            GT_New{igp}=(err(1, 1));
            New_GT=GT_New{igp};
            
            
            
             try
            
              sleep_lost_gp1=[New_GT;num2cell(cell2mat(sortederr(:,12)));num2cell(nan((Ln_gp_max-Ln_gp(igp))+1,1));....
              num2cell(Fly_mean2(7));num2cell(Fly_std2(7));num2cell(Fly_SE2(7))];
            
            sleep_lost_gp=[sleep_lost_gp sleep_lost_gp1];
            sleep_lost_gp_final=[Individual_head sleep_lost_gp];
            
            
            Per_sleep_lost_gp1=[New_GT;num2cell(cell2mat(sortederr(:,13)));num2cell(nan((Ln_gp_max-Ln_gp(igp))+1,1));....
                num2cell(Fly_mean2(8));num2cell(Fly_std2(8));num2cell(Fly_SE2(8))];
            
            Per_sleep_lost_gp=[Per_sleep_lost_gp Per_sleep_lost_gp1];
            Per_sleep_lost_gp_final=[Individual_head Per_sleep_lost_gp];
            
            sleep_gain_gp1=[New_GT;num2cell(cell2mat(sortederr(:,14)));num2cell(nan((Ln_gp_max-Ln_gp(igp))+1,1));....
                num2cell(Fly_mean2(9));num2cell(Fly_std2(9));num2cell(Fly_SE2(9))];
            
            sleep_gain_gp=[sleep_gain_gp sleep_gain_gp1];
            sleep_gain_gp_final=[Individual_head sleep_gain_gp];
            
            Per_sleep_gain_gp1=[New_GT;num2cell(cell2mat(sortederr(:,15)));num2cell(nan((Ln_gp_max-Ln_gp(igp))+1,1));....
                num2cell(Fly_mean2(10));num2cell(Fly_std2(10));num2cell(Fly_SE2(10))];
            
            Per_sleep_gain_gp=[Per_sleep_gain_gp Per_sleep_gain_gp1];
            Per_sleep_gain_gp_final=[Individual_head Per_sleep_gain_gp];
            
            
            Latency_SD_gp1=[New_GT;num2cell(cell2mat(sortederr(:,16)));num2cell(nan((Ln_gp_max-Ln_gp(igp))+1,1));....
                num2cell(Fly_mean2(11));num2cell(Fly_std2(11));num2cell(Fly_SE2(11))];
            
            Latency_SD_gp=[Latency_SD_gp Latency_SD_gp1];
            Latency_SD_gp_final=[Individual_head Latency_SD_gp];
            
            
            
            sleep_gp1=[New_GT;num2cell(cell2mat(sortederr(:,6)));num2cell(nan((Ln_gp_max-Ln_gp(igp))+1,1));....
                num2cell(Fly_mean2(1));num2cell(Fly_std2(1)); num2cell(Fly_SE2(1))];
            
            sleep_gp=[sleep_gp sleep_gp1];
            sleep_gp_final=[Individual_head sleep_gp];
            
            
            Bout_N_gp1=[New_GT;num2cell(cell2mat(sortederr(:,7)));num2cell(nan((Ln_gp_max-Ln_gp(igp))+1,1));....
                num2cell(Fly_mean2(2));num2cell(Fly_std2(2));num2cell(Fly_SE2(2))];
            
            Bout_N_gp=[Bout_N_gp Bout_N_gp1];
            Bout_N_gp_final=[Individual_head Bout_N_gp];
            
            
            Bout_L_gp1=[New_GT;num2cell(cell2mat(sortederr(:,8)));num2cell(nan((Ln_gp_max-Ln_gp(igp))+1,1));....
                num2cell(Fly_mean2(3));num2cell(Fly_std2(3));num2cell(Fly_SE2(3))];
            
            Bout_L_gp=[Bout_L_gp Bout_L_gp1];
            Bout_L_gp_final=[Individual_head Bout_L_gp];
            
            
            Total_act_gp1=[New_GT;num2cell(cell2mat(sortederr(:,9)));num2cell(nan((Ln_gp_max-Ln_gp(igp))+1,1));....
                num2cell(Fly_mean2(4));num2cell(Fly_std2(4));num2cell(Fly_SE2(4))];
            
            Total_act_gp=[Total_act_gp Total_act_gp1];
            Total_act_gp_final=[Individual_head Total_act_gp];
            
            act_count_gp1=[New_GT;num2cell(cell2mat(sortederr(:,10)));num2cell(nan((Ln_gp_max-Ln_gp(igp))+1,1));....
                num2cell(Fly_mean2(5));num2cell(Fly_std2(5));num2cell(Fly_SE2(5))];
            
            act_count_gp=[act_count_gp act_count_gp1];
            act_count_gp_final=[Individual_head act_count_gp];
            
            Latency_gp1=[New_GT;num2cell(cell2mat(sortederr(:,11)));num2cell(nan((Ln_gp_max-Ln_gp(igp))+1,1));....
                num2cell(Fly_mean2(6));num2cell(Fly_std2(6));num2cell(Fly_SE2(6))];
            
            Latency_gp=[Latency_gp Latency_gp1];
            Latency_gp_final=[Individual_head Latency_gp];
            
            
            
            
            catch
            
              sleep_gp1=[New_GT;num2cell(cell2mat(sortederr(:,6)));num2cell(nan((Ln_gp_max-Ln_gp(igp))+1,1));....
                num2cell(Fly_mean2(1));num2cell(Fly_std2(1));num2cell(Fly_SE2(1))];
            
            sleep_gp=[sleep_gp sleep_gp1];
            sleep_gp_final=[Individual_head sleep_gp];
            
            
            Bout_N_gp1=[New_GT;num2cell(cell2mat(sortederr(:,7)));num2cell(nan((Ln_gp_max-Ln_gp(igp))+1,1));....
                num2cell(Fly_mean2(2));num2cell(Fly_std2(2));num2cell(Fly_SE2(2))];
            
            Bout_N_gp=[Bout_N_gp Bout_N_gp1];
            Bout_N_gp_final=[Individual_head Bout_N_gp];
            
            
            Bout_L_gp1=[New_GT;num2cell(cell2mat(sortederr(:,8)));num2cell(nan((Ln_gp_max-Ln_gp(igp))+1,1));....
                num2cell(Fly_mean2(3));num2cell(Fly_std2(3));num2cell(Fly_SE2(3))];
            
            Bout_L_gp=[Bout_L_gp Bout_L_gp1];
            Bout_L_gp_final=[Individual_head Bout_L_gp];
            
            
            Total_act_gp1=[New_GT;num2cell(cell2mat(sortederr(:,9)));num2cell(nan((Ln_gp_max-Ln_gp(igp))+1,1));....
                num2cell(Fly_mean2(4));num2cell(Fly_std2(4));num2cell(Fly_SE2(4))];
            
            Total_act_gp=[Total_act_gp Total_act_gp1];
            Total_act_gp_final=[Individual_head Total_act_gp];
            
            
            act_count_gp1=[New_GT;num2cell(cell2mat(sortederr(:,10)));num2cell(nan((Ln_gp_max-Ln_gp(igp))+1,1));....
                num2cell(Fly_mean2(5));num2cell(Fly_std2(5));num2cell(Fly_SE2(5))];
            
            act_count_gp=[act_count_gp act_count_gp1];
            act_count_gp_final=[Individual_head act_count_gp];
            
            Latency_gp1=[New_GT;num2cell(cell2mat(sortederr(:,11)));num2cell(nan((Ln_gp_max-Ln_gp(igp))+1,1));....
                num2cell(Fly_mean2(6));num2cell(Fly_std2(6));num2cell(Fly_SE2(6))];
            
            Latency_gp=[Latency_gp Latency_gp1];
            Latency_gp_final=[Individual_head Latency_gp];
            
            
            
            
            
             end
            
            %end
            
            
            
            
            
            if strcmp(SD_select2, 'Yes')
            
            %Triage_all=[Triage_all; L_NZ_header;sortederr(:,1:4) sortederr(:,22:24);num2cell(End_sp_T);Estar_T;num2cell(End_sp_T)];xlswrite
            
            End_sp1=NaN*ones(1,(length(Column_headers)+length(Column_headers1)));
            Estar1=cell(1,(length(Column_headers)+length(Column_headers1)));
            Estar1(:)={'******'};
                Result_all=[Result_all;Column_headers Column_headers1;sortederr(:,1:16);sortederr(1,1)  head_mean Title_N num2cell(Fly_mean2);sortederr(1,1)  head_std num2cell(Fly_std2);sortederr(1,1) head_SE num2cell(Fly_SE2);num2cell(End_sp1);Estar1;num2cell(End_sp1)];
                Result_all_SC=[Result_all_SC;sortederr(:,1:16)];
            
                summary_all=[summary_all;Column_headers Column_headers1 ;sortederr(1,1)  head_mean Title_N num2cell(Fly_mean2);sortederr(1,1)  head_std num2cell(Fly_std2); sortederr(1,1) head_SE num2cell(Fly_SE2);num2cell(End_sp1);Estar1;num2cell(End_sp1)];
                %screen_data=[screen_data;Column_headers(1) Column_headers(6:end);sortederr(1,1) num2cell(Fly_mean2)];
                screen_data=[screen_data;sortederr(1,1) num2cell(Fly_mean2)];
            else
            
            
            %Triage_all=[Triage_all;L_NZ_header;sortederr(:,1:4) sortederr(:,17:19);num2cell(End_sp_T);Estar_T;num2cell(End_sp_T)];
            End_sp2=NaN*ones(1,(length(Column_headers)));
            Estar2=cell(1,(length(Column_headers)));
            Estar2(:)={'******'};
               Result_all=[Result_all;Column_headers;sortederr(:,1:11);sortederr(1,1) head_mean Title_N num2cell(Fly_mean2);sortederr(1,1) head_std num2cell(Fly_std2);sortederr(1,1) head_SE num2cell(Fly_SE2);num2cell(End_sp2);Estar2;num2cell(End_sp2)];
               Result_all_SC=[Result_all_SC;sortederr(:,1:11)];
            
               summary_all=[summary_all;Column_headers;sortederr(1,1) head_mean Title_N num2cell(Fly_mean2);sortederr(1,1) head_std num2cell(Fly_std2);sortederr(1,1) head_SE num2cell(Fly_SE2);num2cell(End_sp2);Estar2;num2cell(End_sp2)];
               screen_data=[screen_data;sortederr(1,1) num2cell(Fly_mean2)];
            end
            %screen_data=[screen_data;sortederr(1,1) num2cell(Fly_mean2(1:6))];
            end
            %screen_data=[Column_headers(1) Column_headers(6:end);screen_data];
            %screen_data_num=cell2mat(screen_data);
            
            try
            
            
               screen_data_T=cell2table(screen_data,...
                 "VariableNames", ["GT" "sleep" "BN" "BL" "TA" "AW" "Latency" "SL" "SL_P" "SG" "SG_P" "SD_Latency"]);
            
                screen_data_SL = sortrows(screen_data_T, 'SL');
                screen_data_SLP = sortrows(screen_data_T, 'SL_P');
                screen_data_SG = sortrows(screen_data_T, 'SG');
                screen_data_SGP = sortrows(screen_data_T, 'SG_P');
                screen_data_SDL = sortrows(screen_data_T, 'SD_Latency');
            
                screen_data_sleep = sortrows(screen_data_T, 'sleep');
                screen_data_BN = sortrows(screen_data_T, 'BN');
                screen_data_BL = sortrows(screen_data_T, 'BL');
                screen_data_TA = sortrows(screen_data_T, 'TA');
                screen_data_AW = sortrows(screen_data_T, 'AW');
                screen_data_Latency = sortrows(screen_data_T, 'Latency');
            
                % screen_data_headers={"GenoType_Total_sleep" "Total_sleep" "GenoType_BoutNumber" "BoutNumber" "GenoType_BoutLength" "BoutLength" "GenoType_Total_activity" "Total_activity"...
                %     "GenoType_Activity/waking min" "Activity/waking min" "GenoType_Latency" "Latency"...
                %     "GenoType_Sleep_loss" "Sleep_loss" "GenoType_Percentage_Sleep_loss" "Percentage_Sleep_loss" "GenoType_Sleep_gain" "Sleep_gain" "GenoType_Percentage_Sleep_gain" "Percentage_Sleep_gain" "GenoType_Latency_after_SD" "Latency_after_SD"};
                %
                % screen_data_all=[screen_data_headers;
                % screen_data_sleep.GT num2cell(screen_data_sleep.sleep)...
                % screen_data_BN.GT num2cell(screen_data_BN.BN)...
                % screen_data_BL.GT num2cell(screen_data_BL.BL)...
                % screen_data_TA.GT num2cell(screen_data_TA .TA)...
                % screen_data_AW.GT num2cell(screen_data_AW.AW)...
                % screen_data_Latency.GT num2cell(screen_data_Latency.Latency)...
                % screen_data_SL.GT num2cell(screen_data_SL.SL)...
                % screen_data_SLP.GT num2cell(screen_data_SLP.SL_P)...
                % screen_data_SG.GT num2cell(screen_data_SG.SG)...
                % screen_data_SGP.GT num2cell(screen_data_SGP.SG_P)...
                % screen_data_SDL.GT num2cell(screen_data_SDL.SD_Latency)];
            
            
                   screen_data_headers_New={"Genotype" "Total_sleep"  "BoutNumber"  "BoutLength"  "Total_activity"...
                     "Activity/waking min"  "Latency" "Sleep_change_during_activation"  "Percentage_Sleep_change_during_activation"  "Sleep_change_after activation"  "Percentage_Sleep_change_after activation"  "Latency_after_activation"};
                screen_data=[screen_data_headers_New;screen_data];
            
            
            %**************************************************************************************************************************************************************
                                                             %%% P value calculation for SD experiments starts %%%%%%
            % *************************************************************************************************************************************************************
            
            try
            
            Result_all_SC=cell2table(Result_all_SC,...
                 "VariableNames", ["Genotype" "Run" "Monitor" "channel" "Triage" "Sleep_Day_min_" "BoutNumber_Day" "AverageBoutLength_Day" "TotalActivity_Day" ...
            "Activity_wakingMin"  "Latency_min_" "SleepLost_min_"  "PercentageSleepLost" ...
            "SleepGain_min_" "PercentageSleepGain"  "Latency_SD"]);
            
            % Get data for the reference genotype
            ref_data = Result_all_SC(contains(Result_all_SC.Genotype, ref_genotype), :);
            % Find the p values for the screen data with respect to refernce genotype(control)
            
            % Extract unique genotypes
            genotypes = unique(Result_all_SC.Genotype);
            
            % Initialize result storage
            results = table('Size', [numel(genotypes), 10], 'VariableTypes', {'string', 'double', 'double', 'double','double', 'double', 'double','double', 'double','double'}, ...
                'VariableNames', {'Genotype', 'Sleep_Pval', 'BoutNumber_Pval', 'BoutLength_Pval', 'TotalActivity_Pval', 'Activity_wakingMin_Pval', ...
                'Latency_Pval', 'Sleep_change_during_activataion_Pval', 'Sleep_change_after_activation_Pval', 'Latency_after_activation_Pval'});
            
            % Loop over each genotype to compare with the reference genotype
            for i = 1:numel(genotypes)
                % Get data for the current genotype
                current_data = Result_all_SC(strcmp(Result_all_SC.Genotype, genotypes{i}), :);
            
                    isPresent = any(ismember(ref_data.Genotype, genotypes{i}));
            
            if isPresent
              Sleep_Pval = 1;
              BoutNumber_Pval = 1;
              BoutLength_Pval = 1;
              TotalActivity_Pval = 1;
              Activity_wakingMin_Pval = 1;
              Latency_Pval = 1;
              Sleep_change_during_activataion_Pval = 1;
              percetage_Sleep_change_during_activataion_Pval = 1;
              Sleep_change_after_activation_Pval = 1;
              percentage_Sleep_change_after_activation_Pval = 1;
              Latency_after_activation_Pval=1;
            else
                % Perform t-tests
                [~,Sleep_Pval] = ttest2(ref_data.Sleep_Day_min_, current_data.Sleep_Day_min_);
                [~, BoutNumber_Pval] = ttest2(ref_data.BoutNumber_Day, current_data.BoutNumber_Day);
                [~, BoutLength_Pval] = ttest2(ref_data.AverageBoutLength_Day, current_data.AverageBoutLength_Day);
                [~, TotalActivity_Pval] = ttest2(ref_data.TotalActivity_Day, current_data.TotalActivity_Day);
                [~, Activity_wakingMin_Pval] = ttest2(ref_data.Activity_wakingMin, current_data.Activity_wakingMin);
                [~, Latency_Pval] = ttest2(ref_data.Latency_min_, current_data.Latency_min_);
                [~, Sleep_change_during_activataion_Pval] = ttest2(ref_data.SleepLost_min_, current_data.SleepLost_min_);
                [~, percetage_Sleep_change_during_activataion_Pval] = ttest2(ref_data.PercentageSleepLost, current_data.PercentageSleepLost);
                %[~, percetage_Sleep_change_during_activataion_Pval] = 2;
                [~, Sleep_change_after_activation_Pval] = ttest2(ref_data.SleepGain_min_, current_data.SleepGain_min_);
                [~, percentage_Sleep_change_after_activation_Pval] = ttest2(ref_data.PercentageSleepGain, current_data.PercentageSleepGain);
                [~, Latency_after_activation_Pval] = ttest2(ref_data.Latency_SD, current_data.Latency_SD);
            end
                % Store results
                results.Genotype(i) = genotypes{i};
                results.Sleep_Pval(i) = Sleep_Pval;
                results.BoutNumber_Pval(i) = BoutNumber_Pval;
                results.BoutLength_Pval(i) = BoutLength_Pval;
                results.TotalActivity_Pval(i) = TotalActivity_Pval;
                results.Activity_wakingMin_Pval(i) = Activity_wakingMin_Pval;
                results.Latency_Pval(i) = Latency_Pval;
                results.Sleep_change_during_activataion_Pval(i) = Sleep_change_during_activataion_Pval;
                %results.percetage_Sleep_change_during_activataion_Pval(i) = percetage_Sleep_change_during_activataion_Pval;
                results.Sleep_change_after_activation_Pval(i) = Sleep_change_after_activation_Pval;
                %results.percentage_Sleep_change_after_activation_Pval(i) = percentage_Sleep_change_after_activation_Pval;
                results.Latency_after_activation_Pval(i) = Latency_after_activation_Pval;
            
            end
            
            % Pval_data_headers={"Genotype" "Sleep_Pval"  "BoutNumber_Pval"  "BoutLength_Pval" "TotalActivity_Pval"  "Activity_wakingMin_Pval"  "Latency_Pval"...
            %          "SleepLost_Pval"  "percetage_SleepLoss_Pval"  "SleepGain_Pval"  "percentage_SleepGain_Pval" "Latency_SD_Pval"};
            
                Pval_data_all=results;
            
            end
            
            %**************************************************************************************************************************************************************
                                                                             % END of P value calculation for SD %
            %**************************************************************************************************************************************************************
            
            
            
            catch
                screen_data_T=cell2table(screen_data,...
                 "VariableNames", ["GT" "sleep" "BN" "BL" "TA" "AW" "Latency"]);
            
            screen_data_sleep = sortrows(screen_data_T, 'sleep');
            screen_data_BN = sortrows(screen_data_T, 'BN');
            screen_data_BL = sortrows(screen_data_T, 'BL');
            screen_data_TA = sortrows(screen_data_T, 'TA');
            screen_data_AW = sortrows(screen_data_T, 'AW');
            screen_data_Latency = sortrows(screen_data_T, 'Latency');
            
            %  screen_data_headers={"GenoType_Total_sleep" "Total_sleep" "GenoType_BoutNumber" "BoutNumber" "GenoType_BoutLength" "BoutLength" "GenoType_Total_activity" "Total_activity"...
            %         "GenoType_Activity/waking min" "Activity/waking min" "GenoType_Latency" "Latency"};
            %
            % screen_data_all=[screen_data_headers;
            %     screen_data_sleep.GT num2cell(screen_data_sleep.sleep)...
            %     screen_data_BN.GT num2cell(screen_data_BN.BN)...
            %     screen_data_BL.GT num2cell(screen_data_BL.BL)...
            %     screen_data_TA.GT num2cell(screen_data_TA .TA)...
            %     screen_data_AW.GT num2cell(screen_data_AW.AW)...
            %     screen_data_Latency.GT num2cell(screen_data_Latency.Latency)];
            
            
                   screen_data_headers_New={"Genotype" "Total_sleep"  "BoutNumber"  "BoutLength"  "Total_activity"...
                     "Activity/waking min"  "Latency"};
            
                   screen_data=[screen_data_headers_New;screen_data];
            
                %Result_all_SC=[Column_headers;Result_all_SC];
            
            
            %**************************************************************************************************************************************************************
                                                             %%% P value calculation for non-SD experiments starts %%%%%%
            % *************************************************************************************************************************************************************
            
            try
            
            Result_all_SC=cell2table(Result_all_SC,...
                 "VariableNames", ["Genotype" "Run" "Monitor" "channel" "Triage" "Sleep_Day_min_" "BoutNumber_Day" "AverageBoutLength_Day" "TotalActivity_Day" ...
            "Activity_wakingMin"  "Latency_min_"]);
            
            % Get data for the reference genotype
            ref_data = Result_all_SC(contains(Result_all_SC.Genotype, ref_genotype), :);
            
            
            % Find the p values for the screen data with respect to refernce genotype(control)
            
            % Extract unique genotypes
            genotypes = unique(Result_all_SC.Genotype);
            
            % Initialize result storage
            results = table('Size', [numel(genotypes), 7], 'VariableTypes', {'string', 'double', 'double', 'double','double', 'double', 'double'}, ...
                'VariableNames', {'Genotype', 'Sleep_Pval', 'BoutNumber_Pval', 'BoutLength_Pval', 'TotalActivity_Pval', 'Activity_wakingMin_Pval', ...
                'Latency_Pval'});
            
            
            % Loop over each genotype to compare with the reference genotype
            for i = 1:numel(genotypes)
                % Get data for the current genotype
                current_data = Result_all_SC(strcmp(Result_all_SC.Genotype, genotypes{i}), :);
            
                   isPresent = any(ismember(ref_data.Genotype, genotypes{i}));
            
            if isPresent
              Sleep_Pval = 1;
              BoutNumber_Pval = 1;
              BoutLength_Pval = 1;
              TotalActivity_Pval = 1;
              Activity_wakingMin_Pval = 1;
              Latency_Pval = 1;
            
            else
            
                % Perform t-tests
                [~, Sleep_Pval] = ttest2(ref_data.Sleep_Day_min_, current_data.Sleep_Day_min_);
                [~, BoutNumber_Pval] = ttest2(ref_data.BoutNumber_Day, current_data.BoutNumber_Day);
                [~, BoutLength_Pval] = ttest2(ref_data.AverageBoutLength_Day, current_data.AverageBoutLength_Day);
                [~, TotalActivity_Pval] = ttest2(ref_data.TotalActivity_Day, current_data.TotalActivity_Day);
                [~, Activity_wakingMin_Pval] = ttest2(ref_data.Activity_wakingMin, current_data.Activity_wakingMin);
                [~, Latency_Pval] = ttest2(ref_data.Latency_min_, current_data.Latency_min_);
            
            end
                % Store results
                results.Genotype(i) = genotypes{i};
                results.Sleep_Pval(i) = Sleep_Pval;
                results.BoutNumber_Pval(i) = BoutNumber_Pval;
                results.BoutLength_Pval(i) = BoutLength_Pval;
                results.TotalActivity_Pval(i) = TotalActivity_Pval;
                results.Activity_wakingMin_Pval(i) = Activity_wakingMin_Pval;
                results.Latency_Pval(i) = Latency_Pval;
            
            end
            
            %Pval_data_headers={"Genotype" "Sleep_Pval" "BoutNumber_Pval" "BoutLength_Pval"  "TotalActivity_Pval" "Activity_wakingMin_Pval" "Latency_Pval"};
            
                Pval_data_all=results;
            
            %**************************************************************************************************************************************************************
                                                             %%% P value calculation for non SD experiments END %%%%%%
            % *************************************************************************************************************************************************************
            
            end
            
            end
            
            
            
            if strcmp(Lifespan_select2, 'Yes')
            
            %% Lifespan result sorting and combining same genotype
            
            
            
            [~,~,X_LNZ]=unique(R_all_LNZ(:,1));
            C_LNZ = accumarray(X_LNZ,1:size(R_all_LNZ,1),[],@(r){R_all_LNZ(r,:)});
            
            
            [nr,nc]=size(C_LNZ);
            
            for igp=1:nr
            
            err_LNZ=C_LNZ{igp};
            
            
            [~, order] = sort(cell2mat(err_LNZ(:, 1)));
            sortederr_LNZ =  sortrows(err_LNZ, [1 2 3 4]);% sort in the order of column 1 and 2
            
            L_NZ_header={'Genotype', 'Run number','Monitor number','Channel', 'Last activity date', 'Last activity time', 'Life span (days)'};
            End_sp_T=NaN*ones(1,7);
            Estar_T=cell(1,7);
            Estar_T(:)={'******'};
            Triage_all=[Triage_all; L_NZ_header;sortederr_LNZ;num2cell(End_sp_T);Estar_T;num2cell(End_sp_T)];
            
            
            LF_GT{igp}=(err_LNZ(1, 1));
            New_GT_LF=LF_GT{igp};
            
            Life_gp1=[New_GT_LF;num2cell(cell2mat(sortederr_LNZ(:,7)));num2cell(nan((Ln_gp_max-Ln_gp(igp))+1,1))];
            
            Life_gp=[Life_gp Life_gp1];
            end
            
            end
            
            
            
            
            
            %% Combaining L & D results togther
            
            if strcmp(SL_select2, 'Yes')
            
            [~,~,X_L_D]=unique(R_all_L_D(:,1));
            C_L_D = accumarray(X_L_D,1:size(R_all_L_D,1),[],@(r){R_all_L_D(r,:)});
            [nr_L_D,nc_L_D]=size(C_L_D);
            
            
            
            Column_headers_L_D={'Genotype', 'Run number','Monitor number','Channel',Triage_condn, 'Sleep/L (min)', 'Sleep/D (min)', 'Bout number/L','Bout number/D', 'Average bout length/L', 'Average bout length/D', 'Total activity/L', 'Total activity/D', 'Activity/waking min/L', 'Activity/waking min /D'};
            End_sp2_L_D=NaN*ones(1,(length(Column_headers_L_D)));
            Estar2_L_D=cell(1,(length(Column_headers_L_D)));
            Estar2_L_D(:)={'******'};
            
            
            for igp_L_D=1:nr_L_D
            
            err_L_D=C_L_D{igp_L_D};
            
            
            [~, order] = sort(cell2mat(err_L_D(:, 1)));
            sortederr_L_D =  sortrows(err_L_D, [1 2 3 4]);% sort in the order of column 1 and 2
            
            
            N_active=nnz(~(cell2mat(sortederr_L_D(:, 5))));
            Title_N=[Title5 ' ' num2str(N_active)];
            
                L_D_mean2=nanmean(cell2mat(sortederr_L_D(:, 6:end)),1);
                L_D_std2=nanstd(cell2mat(sortederr_L_D(:, 6:end)),0,1);
                L_D_SE2= L_D_std2/sqrt(N_active);
            
            
            
               Result_all_L_D=[Result_all_L_D;Column_headers_L_D;sortederr_L_D;sortederr_L_D(1,1) head_mean Title_N num2cell(L_D_mean2);sortederr_L_D(1,1) head_std num2cell(L_D_std2);sortederr_L_D(1,1) head_SE num2cell(L_D_SE2);num2cell(End_sp2_L_D);Estar2_L_D;num2cell(End_sp2_L_D)];
            
            
            
            end
            end
            
            
            %try
            
                %% Combaining anticipation results togther
            
            try
            [~,~,X_anti]=unique(R_all_anti(:,1));
            C_anti = accumarray(X_anti,1:size(R_all_anti,1),[],@(r){R_all_anti(r,:)});
            [nr_anti,nc_anti]=size(C_anti);
            
            
             Ln_anti_gp=[];
            length(C_anti{1});
             for i_split=1:nr_anti
             [nr_anti_gp, nc_anti_gp]=size(C_anti{i_split});
             Ln_anti_gp(i_split)=nr_anti_gp;
             end
            
            Ln_anti_gp_max=max(Ln_anti_gp);
            
            Individual_head_11=NaN*ones(Ln_anti_gp_max+2,1);
            Individual_head_mean={'Mean'};
            Individual_head_SD={'SD'};
            Individual_head_SE={'SE'};
            Individual_head_anti =[num2cell(Individual_head_11);Individual_head_mean;Individual_head_SD;Individual_head_SE];
            
            
            
            
            Column_headers_anti={'Genotype', 'Run number','Monitor number','Channel',Triage_condn_ed,  'Morning Anticipation', 'Evening Anticipation', 'Morning Anticipation phase',  'Evening Anticipation phase'};
            End_sp2_anti=NaN*ones(1,(length(Column_headers_anti)));
            Estar2_anti=cell(1,(length(Column_headers_anti)));
            Estar2_anti(:)={'******'};
            
            
            for igp_anti=1:nr_anti
            
            err_anti=C_anti{igp_anti};
            
            
            [~, order] = sort(cell2mat(err_anti(:, 1)));
            sortederr_anti =  sortrows(err_anti, [1 2 3 4]);% sort in the order of column 1 and 2
            
            
            N_active=nnz(~(cell2mat(sortederr_anti(:, 5))));
            Title_N=[Title5 ' ' num2str(N_active)];
            
               anti_mean2=nanmean(cell2mat(sortederr_anti(:, 6:end)),1);
                anti_std2=nanstd(cell2mat(sortederr_anti(:, 6:end)),0,1);
                anti_SE2= anti_std2/sqrt(N_active);
            
            
            
            
               Result_all_anti=[Result_all_anti;Column_headers_anti;sortederr_anti;sortederr_anti(1,1) head_mean Title_N num2cell(anti_mean2);sortederr_anti(1,1) head_std num2cell(anti_std2);sortederr_anti(1,1) head_SE num2cell(anti_SE2);num2cell(End_sp2_anti);Estar2_anti;num2cell(End_sp2_anti)];
               Result_all_anti_SC=[Result_all_anti_SC;sortederr_anti];
               %% screening data for anticipation
               screen_data_anticipation=[screen_data_anticipation;sortederr_anti(1,1) num2cell(anti_mean2)];
            
            
            %% *************************************************************
            GT_New_anti{igp_anti}=(err_anti(1, 1));
            New_GT_anti=GT_New_anti{igp_anti};
            
            
            Manti_gp1=[New_GT_anti;num2cell(cell2mat(sortederr_anti(:,6)));num2cell(nan((Ln_anti_gp_max-Ln_anti_gp(igp_anti))+1,1));....
                num2cell(anti_mean2(1));num2cell(anti_std2(1));num2cell(anti_SE2(1))];
            
            Manti_gp=[Manti_gp Manti_gp1];
            Manti_gp_final=[Individual_head_anti Manti_gp];
            
            Eanti_gp1=[New_GT_anti;num2cell(cell2mat(sortederr_anti(:,7)));num2cell(nan((Ln_anti_gp_max-Ln_anti_gp(igp_anti))+1,1));....
                num2cell(anti_mean2(2));num2cell(anti_std2(2));num2cell(anti_SE2(2))];
            
            Eanti_gp=[Eanti_gp Eanti_gp1];
            Eanti_gp_final=[Individual_head_anti Eanti_gp];
            
            Manti_ph_gp1=[New_GT_anti;num2cell(cell2mat(sortederr_anti(:,8)));num2cell(nan((Ln_anti_gp_max-Ln_anti_gp(igp_anti))+1,1));....
                num2cell(anti_mean2(3));num2cell(anti_std2(3));num2cell(anti_SE2(3))];
            
            Manti_ph_gp=[Manti_ph_gp Manti_ph_gp1];
            Manti_ph_gp_final=[Individual_head_anti Manti_ph_gp];
            
            Eanti_ph_gp1=[New_GT_anti;num2cell(cell2mat(sortederr_anti(:,9)));num2cell(nan((Ln_anti_gp_max-Ln_anti_gp(igp_anti))+1,1));....
                num2cell(anti_mean2(4));num2cell(anti_std2(4));num2cell(anti_SE2(4))];
            
            Eanti_ph_gp=[Eanti_ph_gp Eanti_ph_gp1];
            Eanti_ph_gp_final=[Individual_head_anti Eanti_ph_gp];
            
            end
            
            % screening data for anticipation
            
            screen_data_anticipation_T=cell2table(screen_data_anticipation,...
                 "VariableNames", ["GT" "MA" "EA" "MAP" "EAP"]);
            
            
                screen_data_anticipation_MA = sortrows(screen_data_anticipation_T, 'MA');
                screen_data_anticipation_EA = sortrows(screen_data_anticipation_T, 'EA');
                screen_data_anticipation_MAP = sortrows(screen_data_anticipation_T, 'MAP');
                screen_data_anticipation_EAP = sortrows(screen_data_anticipation_T, 'EAP');
            
            
                screen_data_anticipation_headers_New={"Genotype" "Morning_Anticipation" "Evening_Anticipation" "Morning_anticipation_phase" "Evening_anticipation_phase"};
                screen_data_anticipation=[screen_data_anticipation_headers_New;screen_data_anticipation];
            
                % Combine anticipation screening data with other screening data
                screen_data=[screen_data, screen_data_anticipation(:,2:5)];
            
            
            end
            
            
            try
                  % Combine periodogram screening data with other screening data
            screen_data=[screen_data, screen_data_period(:,2:3)];
            end
            
            
            
            
            % ################## P values of anticipation screening results start    % %*************************************************************
            
            try
            Result_all_anti_SC=cell2table(Result_all_anti_SC,...
                  "VariableNames", ["Genotype" "Run" "Monitor" "channel" "Triage" "MA" "EA" "MAP" "EAP"]);
            
            % Get data for the reference genotype
            ref_data_anti = Result_all_anti_SC(contains(Result_all_anti_SC.Genotype, ref_genotype), :);
            % Extract unique genotypes
            genotypes_anti = unique(Result_all_anti_SC.Genotype);
            
            
            % Initialize result storage
            results_anti = table('Size', [numel(genotypes_anti), 5], 'VariableTypes', {'string', 'double', 'double','double', 'double'}, ...
                'VariableNames', {'Genotype', 'Morning_anticipation_Pval', 'Evening_anticipation_Pval', 'Morning_anticipation_Phase_Pval', 'Evening_anticipation_Phase_Pval'});
            
            
            % Loop over each genotype to compare with the reference genotype
            for i = 1:numel(genotypes_anti)
                % Get data for the current genotype
                current_data_anti = Result_all_anti_SC(strcmp(Result_all_anti_SC.Genotype, genotypes_anti{i}), :);
            
                    isPresent = any(ismember(ref_data_anti.Genotype, genotypes_anti{i}));
            
            if isPresent
              Morning_anticipation_Pval = 1;
              Evening_anticipation_Pval = 1;
              Morning_anticipation_Phase_Pval = 1;
              Evening_anticipation_Phase_Pval = 1;
            
            else
                % Perform t-tests
                [~,Morning_anticipation_Pval] = ttest2(ref_data_anti.MA, current_data_anti.MA);
                [~, Evening_anticipation_Pval] = ttest2(ref_data_anti.EA, current_data_anti.EA);
                [~,Morning_anticipation_Phase_Pval] = ttest2(ref_data_anti.MAP, current_data_anti.MAP);
                [~, Evening_anticipation_Phase_Pval] = ttest2(ref_data_anti.EAP, current_data_anti.EAP);
            
            end
                % Store results
                results_anti.Genotype(i) = genotypes_anti{i};
                results_anti.Morning_anticipation_Pval(i) = Morning_anticipation_Pval;
                results_anti.Evening_anticipation_Pval(i) = Evening_anticipation_Pval;
                results_anti.Morning_anticipation_Phase_Pval(i) = Morning_anticipation_Phase_Pval;
                results_anti.Evening_anticipation_Phase_Pval(i) = Evening_anticipation_Phase_Pval;
            end
            
            
                Pval_data_all_anti=results_anti;
            
                % ################## P values of anticipation screening results Ends % % **************************************************
            
                % Combine anticipation P val with other P values
                Pval_data_all=[Pval_data_all, Pval_data_all_anti(:, 2:5)];
            
            end
            
            try
                % Combine periodogram Pval with other P values
            Pval_data_all=[Pval_data_all, Pval_data_all_per(:, 2:3)];
            
            end
            
            
            %% ############################  Combine screening results and P-value together % *******************************************
            
            try
            % convert P value result to cell and add variable name as one column
            Pval_data_all_cell=table2cell(Pval_data_all);
            Pval_data_all_cell=[Pval_data_all.Properties.VariableNames;Pval_data_all_cell];
            
            % Get the size of the arrays
            [numRows, numCols] = size(Pval_data_all_cell);
            
            % Initialize the resulting cell array with twice the number of columns
            SC_Pval = cell(numRows, 2 * numCols);
            
            % Interleave the columns of A and B
            for col = 1:numCols
                SC_Pval(:, 2*col-1) = screen_data(:, col);  % Odd columns from screening data
                SC_Pval(:, 2*col) = Pval_data_all_cell(:, col);    % Even columns from Pvalue
            end
            
            SC_Pval(:,1)=[]; % Remove reduntant genotype column
            
            end
            
            %% ############################  Combine screening results and P-value together Finished % *******************************************
            
            
            
            %%   **************************************************************************************************************************
            
            %%   #################   saving screening results and creating screen data plot begins  #####################################
            
            
            path = handles.path;
            myfolder = 'Screening_results' ;   % new folder name
            folder = mkdir([path,Project_name,filesep,myfolder]);
            path  = [path,Project_name,filesep,myfolder] ;
            
            %save_location_analysis_par=fullfile(path, strcat(Project_name,'screening_data.xls.xls'));
            
            save_location_screening=[path,filesep, strcat(Project_name,'_screening_data.xls')] ;
            
            %writecell(screen_data_all, save_location_screening);
            
            try
                writecell(SC_Pval, save_location_screening);
            catch
            writecell(screen_data, save_location_screening);
            end
            
            % try
            % save_location_Pval=[path,filesep, strcat(Project_name,'_Hits_pVal.xls')] ;
            % writetable(Pval_data_all, save_location_Pval);
            % end
            
            
            if strcmp(SD_select2, 'Yes')
            
            FigName   = 'sleep_Loss';
            figure('NumberTitle', 'off', 'Name',FigName);
            plot(screen_data_SL.SL, '.', 'MarkerSize',20)
            ylabel('sleep Loss (min)');
            % set(gca,'XTick',1:numel(screen_data_SL.GT));
            % set(gca,'XTickLabel',screen_data_SL.GT);
            set(gca,'XTick',[]);
            
            h=gca;
                     temp1=[path,filesep,FigName,'.png'];
                   saveas(gca,temp1);
                   close(FigName);
            
            
            
                FigName   = 'sleep_Rebound';
            figure('NumberTitle', 'off', 'Name',FigName);
            plot(screen_data_SG.SG, '.', 'MarkerSize',20)
            ylabel('sleep Rebound (min)');
            % set(gca,'XTick',1:numel(screen_data_SG.GT));
            % set(gca,'XTickLabel',screen_data_SG.GT);
            h=gca;
            %h. XTickLabelRotation=90;
            set(gca,'XTick',[]);
            
            
                     temp1=[path,filesep,FigName,'.png'];
                           try
                    saveas(gca,temp1);
                    catch
                        disp("Error!! Please check genotype name. Don't use backslash '\' or '>'  for genotype naming if you need figures")
                                        set(handles.M_box,'String',[oldmsgs;{"Error!! Please check genotype name.  Don't use backslash '\' or '>'  for genotype naming if you need figures"}] );drawnow
            
                                        return
                    end
                   close(FigName);
            
            end
            
            
            
            FigName   = 'sleep';
            figure('NumberTitle', 'off', 'Name',FigName);
            plot(screen_data_sleep.sleep, '.', 'MarkerSize',20)
            ylabel('sleep/Day (min)');
            % set(gca,'XTick',1:numel(screen_data_sleep.GT));
            % set(gca,'XTickLabel',screen_data_sleep.GT);
            h=gca;
            % h. XTickLabelRotation=90;
            set(gca,'XTick',[]);
            
            
                     temp1=[path,filesep,FigName,'.png'];
                   saveas(gca,temp1);
                   close(FigName);
            
            
            
            
            
            FigName   = 'Bout_number';
            figure('NumberTitle', 'off', 'Name',FigName);
            
            plot(screen_data_BN.BN, '.', 'MarkerSize',20)
            ylabel('Bout number/Day');
            % set(gca,'XTick',1:numel(screen_data_BN.GT));
            % set(gca,'XTickLabel',screen_data_BN.GT);
            h=gca;
            % h. XTickLabelRotation=90;
            set(gca,'XTick',[]);
            
            
                    temp1=[path,filesep,FigName,'.png'];
                            try
                    saveas(gca,temp1);
                    catch
                        disp("Error!! Please check genotype name. Don't use backslash '\' or '>' for genotype naming if you need figures")
                                        set(handles.M_box,'String',[oldmsgs;{"Error!! Please check genotype name.  Don't use backslash '\' or '>' for genotype naming if you need figures"}] );drawnow
            
                                        return
                    end
              close(FigName);
            
            
            
             FigName   = 'Bout_Length';
            figure('NumberTitle', 'off', 'Name',FigName);
            
            plot(screen_data_BL.BL, '.', 'MarkerSize',20)
            ylabel('Bout Length / Day');
            % set(gca,'XTick',1:numel(screen_data_BL.GT));
            % set(gca,'XTickLabel',screen_data_BL.GT);
            h=gca;
            % h. XTickLabelRotation=90;
            set(gca,'XTick',[]);
            
            
                    temp1=[path,filesep,FigName,'.png'];
                    saveas(gca,temp1);
              close(FigName);
            
            
            
             FigName   = 'Total_activity';
            figure('NumberTitle', 'off', 'Name',FigName);
            
            plot(screen_data_TA.TA, '.', 'MarkerSize',20);
            ylabel('Total activity / Day');
            h=gca;
            set(gca,'XTick',[]);
            
                    temp1=[path,filesep,FigName,'.png'];
                    saveas(gca,temp1);
            
              close(FigName);
            
            
            FigName   = 'Activity_Waking_min';
            figure('NumberTitle', 'off', 'Name',FigName);
            
            plot(screen_data_TA.TA, '.', 'MarkerSize',20)
            ylabel('activity / Waking min');
            % set(gca,'XTick',1:numel(screen_data_AW.GT));
            % set(gca,'XTickLabel',screen_data_AW.GT);
            h=gca;
            % h. XTickLabelRotation=90;
            set(gca,'XTick',[]);
            
                    temp1=[path,filesep,FigName,'.png'];
                    saveas(gca,temp1);
              close(FigName);
            
            
            
            FigName   = 'Latency';
            figure('NumberTitle', 'off', 'Name',FigName);
            
            plot(screen_data_Latency.Latency, '.', 'MarkerSize',20)
            ylabel('Latency (min)');
            % set(gca,'XTick',1:numel(screen_data_Latency.GT));
            % set(gca,'XTickLabel',screen_data_Latency.GT);
            h=gca;
            % h. XTickLabelRotation=90;
            set(gca,'XTick',[]);
            
                    temp1=[path,filesep,FigName,'.png'];
                    saveas(gca,temp1);
              close(FigName);
            
            
            end
            
            
            % ***************************** anticipation screening plot Begins ***************************************************************************
            
            try
            screen_data_anticipation_MA.MA;
            FigName   = 'Morning_anticipation';
            figure('NumberTitle', 'off', 'Name',FigName);
            plot(screen_data_anticipation_MA.MA, '.', 'MarkerSize',20)
            ylabel('Morning Anticipation');
            h=gca;
            set(gca,'XTick',[]);
                    temp1=[path,filesep,FigName,'.png'];
                    saveas(gca,temp1);
              close(FigName);
            
            
            
              FigName   = 'Evening_anticipation';
            figure('NumberTitle', 'off', 'Name',FigName);
            plot(screen_data_anticipation_EA.EA, '.', 'MarkerSize',20)
            ylabel('Evening Anticipation');
            h=gca;
            set(gca,'XTick',[]);
                    temp1=[path,filesep,FigName,'.png'];
                    saveas(gca,temp1);
              close(FigName);
            
            
            FigName   = 'Morning_anticipation phase';
            figure('NumberTitle', 'off', 'Name',FigName);
            plot(screen_data_anticipation_MAP.MAP, '.', 'MarkerSize',20)
            ylabel('Morning Anticipation phase');
            h=gca;
            set(gca,'XTick',[]);
                    temp1=[path,filesep,FigName,'.png'];
                    saveas(gca,temp1);
              close(FigName);
            
            
              FigName   = 'Evening_anticipation_phase';
            figure('NumberTitle', 'off', 'Name',FigName);
            plot(screen_data_anticipation_EAP.EAP, '.', 'MarkerSize',20)
            ylabel('Evening Anticipation phase');
            h=gca;
            set(gca,'XTick',[]);
                    temp1=[path,filesep,FigName,'.png'];
                    saveas(gca,temp1);
              close(FigName);
            
            end
            
              %  Anticipation screening plot ENDS ***************************************************************************************************
            
            
            
            % ***************************** periodogram screening plot Begins ***************************************************************************
            
            try
            screen_data_period_T_P_S.P_S;
            FigName   = 'P-S';
            figure('NumberTitle', 'off', 'Name',FigName);
            plot(screen_data_period_T_P_S.P_S, '.', 'MarkerSize',20)
            ylabel('P-S');
            h=gca;
            set(gca,'XTick',[]);
                    temp1=[path,filesep,FigName,'.png'];
                    saveas(gca,temp1);
              close(FigName);
            
            
              FigName   = 'Period_P-S';
            figure('NumberTitle', 'off', 'Name',FigName);
            plot(screen_data_period_T_Period.Period_P_S, '.', 'MarkerSize',20)
            ylabel('Period (P-S>0)');
            h=gca;
            set(gca,'XTick',[]);
                    temp1=[path,filesep,FigName,'.png'];
                    saveas(gca,temp1);
              close(FigName);
            
            end
            
            % ***************************** periodogram screening plot Ends ***************************************************************************
            
            %% #################     Screen data plot ENDs   ####################################################################################################################
            
            
            
            
            
            
            
            %% ********************************************************************
            
            try
            [~,~,X_ed_TW]=unique(R_all_ed_TW(:,1));
            C_ed_TW = accumarray(X_ed_TW,1:size(R_all_ed_TW,1),[],@(r){R_all_ed_TW(r,:)});
            [nr_ed_TW,nc_ed_TW]=size(C_ed_TW);
            
            
            
            Column_headers_ed_TW={'Genotype', 'Run number','Monitor number','Channel',Triage_condn_ed};
            End_sp2_ed_TW=NaN*ones(1,(length(Column_headers_ed_TW)));
            Estar2_ed_TW=cell(1,(length(Column_headers_ed_TW)));
            Estar2_ed_TW(:)={'******'};
            
            
            for igp_ed_TW=1:nr_ed_TW
            
            err_ed_TW=C_ed_TW{igp_ed_TW};
            
            
            [~, order] = sort(cell2mat(err_ed_TW(:, 1)));
            sortederr_ed_TW =  sortrows(err_ed_TW, [1 2 3 4]);% sort in the order of column 1 and 2
            
            
            N_active=nnz(~(cell2mat(sortederr_ed_TW(:, 5))));
            Title_N=[Title5 ' ' num2str(N_active)];
            
            %    anti_mean2=nanmean(cell2mat(sortederr_ed_TW(:, 6:end)),1);
            %     anti_std2=nanstd(cell2mat(sortederr_ed_TW(:, 6:end)),0,1);
            %     anti_SE2= anti_std2/sqrt(N_active);
            
            
            
               Result_all_ed_TW=[Result_all_ed_TW;Column_headers_ed_TW;sortederr_ed_TW;sortederr_ed_TW(1,1) head_mean Title_N ;sortederr_ed_TW(1,1) head_std ;sortederr_ed_TW(1,1) head_SE ;num2cell(End_sp2_ed_TW);Estar2_ed_TW;num2cell(End_sp2_ed_TW)];
            
            
            end
            end
            
            
            %% ********************************************************************
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            %% ******************************************************************
            
            try
            
            Title_sleep{1}=['Sleep per ' num2str(N_inval_hr) 'hr'];
            Title_boutN{1}=['Bout number per ' num2str(N_inval_hr) 'hr'];
            Title_boutL{1}=['Average bout length per ' num2str(N_inval_hr) 'hr'];
            Title_activity{1}=['Total activity per ' num2str(N_inval_hr) 'hr'];
            
            [~,~,X_Sleep]=unique(Sleep_interval1(:,1));
            [~,~,X_BoutN]=unique(BoutN_interval1(:,1));
            [~,~,X_BoutL]=unique(BoutL_interval1(:,1));
            [~,~,X_Activity]=unique(Activity_interval1(:,1));
            
            
            C_Sleep = accumarray(X_Sleep,1:size(Sleep_interval1,1),[],@(r){Sleep_interval1(r,:)});
            C_BoutN = accumarray(X_BoutN,1:size(BoutN_interval1,1),[],@(r){BoutN_interval1(r,:)});
            C_BoutL = accumarray(X_BoutL,1:size(BoutL_interval1,1),[],@(r){BoutL_interval1(r,:)});
            C_Activity = accumarray(X_Activity,1:size(Activity_interval1,1),[],@(r){Activity_interval1(r,:)});
            
            
            
            
            
            
            
            [nr_sleep,nc_sleep]=size(C_Sleep);
            
            for i_sleep=1:nr_sleep
            
            
             err_Sleep=C_Sleep{i_sleep};
             err_BoutN=C_BoutN{i_sleep};
             err_BoutL=C_BoutL{i_sleep};
             err_Activity=C_Activity{i_sleep};
            
            
            
            [~, order_Sleep] = sort(cell2mat(err_Sleep(:, 1)));
            sortederr_Sleep =  sortrows(err_Sleep, [1 2 3 4]);% sort in the order of column 1,2,3,4
            
            [~, order_BoutN] = sort(cell2mat(err_BoutN(:, 1)));
            sortederr_BoutN =  sortrows(err_BoutN, [1 2 3 4]);% sort in the order of column 1 and 2
            
            [~, order_BoutL] = sort(cell2mat(err_BoutL(:, 1)));
            sortederr_BoutL =  sortrows(err_BoutL, [1 2 3 4]);% sort in the order of column 1 and 2
            
            [~, order_Activity] = sort(cell2mat(err_Activity(:, 1)));
            sortederr_Activity =  sortrows(err_Activity, [1 2 3 4]);% sort in the order of column 1 and 2
            
            
            %      N_active=nnz(~(cell2mat(sortederr_Sleep(:, 5))));
            %      Title_N=[Title5 ' ' num2str(N_active)];
            
            
            
            
            
                Fly_mean_sleep=nanmean(cell2mat(sortederr_Sleep(:, 6:end)),1);
                Fly_std_sleep=nanstd(cell2mat(sortederr_Sleep(:, 6:end)),0,1);
            
                Fly_mean_BoutN=nanmean(cell2mat(sortederr_BoutN(:, 6:end)),1);
                Fly_std_BoutN=nanstd(cell2mat(sortederr_BoutN(:, 6:end)),0,1);
            
                Fly_mean_BoutL=nanmean(cell2mat(sortederr_BoutL(:, 6:end)),1);
                Fly_std_BoutL=nanstd(cell2mat(sortederr_BoutL(:, 6:end)),0,1);
            
                Fly_mean_Activity=nanmean(cell2mat(sortederr_Activity(:, 6:end)),1);
                Fly_std_Activity=nanstd(cell2mat(sortederr_Activity(:, 6:end)),0,1);
            
            
            
            
            [slp_nr, slp_nc]=size(Sleep_interval1);
            End_sp3=NaN*ones(1,slp_nc);
            Estar3=cell(1,slp_nc);
            Estar3(:)={'******'};
            
            
            
            Sleep_interval=[Sleep_interval;Column_headers_inval num2cell(t_array_inval(1:end-1));sortederr_Sleep ;sortederr_Sleep(1,1) head_mean Title_N num2cell(Fly_mean_sleep);sortederr_Sleep(1,1) head_std num2cell(Fly_std_sleep);num2cell(End_sp3);Estar3;num2cell(End_sp3)];
            BoutN_interval=[BoutN_interval;Column_headers_inval num2cell(t_array_inval(1:end-1));sortederr_BoutN;sortederr_BoutN(1,1) head_mean Title_N num2cell(Fly_mean_BoutN);sortederr_BoutN(1,1) head_std num2cell(Fly_std_BoutN);num2cell(End_sp3);Estar3;num2cell(End_sp3)];
            BoutL_interval=[BoutL_interval;Column_headers_inval num2cell(t_array_inval(1:end-1));sortederr_BoutL;sortederr_BoutL(1,1) head_mean Title_N num2cell(Fly_mean_BoutL);sortederr_BoutL(1,1) head_std num2cell(Fly_std_BoutL);num2cell(End_sp3);Estar3;num2cell(End_sp3)];
            Activity_interval=[Activity_interval;Column_headers_inval num2cell(t_array_inval(1:end-1));sortederr_Activity;sortederr_Activity(1,1) head_mean Title_N num2cell(Fly_mean_Activity);sortederr_Activity(1,1) head_std num2cell(Fly_std_Activity);num2cell(End_sp3);Estar3;num2cell(End_sp3)];
            end
            
            end
            
            
            
            
            
            % Combaining eduction result
            
            try
            [~,~,X_R_ed_Ind]=unique(R_ed_Ind(:,1));
            C_R_ed_Ind = accumarray(X_R_ed_Ind,1:size(R_ed_Ind,1),[],@(r){R_ed_Ind(r,:)});
            end
            
            
            
            try
            
            [~,~,X_eduction]=unique(R_ed(:,1));
            [~,~,X_ed_Ndays]=unique(R_ed_Ndays(:,1));
            
            
            
            C_eduction = accumarray(X_eduction,1:size(R_ed,1),[],@(r){R_ed(r,:)});
            C_Ndays = accumarray(X_ed_Ndays,1:size(R_ed_Ndays,1),[],@(r){R_ed_Ndays(r,:)});
            
            
            
            
            
            
            
            
            
            %% If the individual fly is required or not
            
            
            [nr_eduction,nc_eduction]=size(C_eduction);
            
            
            
            for igp_eduction=1:nr_eduction
            
            err_eduction=C_eduction{igp_eduction};
            err_Ndays=C_Ndays{igp_eduction};
            % err_Day_by_Day=C_Day_by_Day{igp_eduction};
            
            
            
            
            [~, order] = sort(cell2mat(err_eduction(:, 1)));
            sortederr_eduction =  sortrows(err_eduction, [1 2 3]);% sort in the order of column 1 and 2
            
            
            
            [~, order_Nd] = sort(cell2mat(err_Ndays(:, 1)));
            sortederr_Ndays =  sortrows(err_Ndays, [1 2 3]);% sort in the order of column 1 and 2
            
            % Normalizing individual fly eduction result(N days)
            try
            err_Ind=C_R_ed_Ind{igp_eduction};
            [~, order_Ind] = sort(cell2mat(err_Ind(:, 1)));
            sortederr_Ind =  sortrows(err_Ind, [1 2 3 4]);% sort in the order of column 1 and 2
            Ind_Ed_Normalized=(cell2mat(sortederr_Ind(:, 6:end)))./nanmean((cell2mat(sortederr_Ind(:, 6:end))),2);
            
            %Ind_Ed_Normalized_with_title=[Run_1 Monitor_1 row_header Triaged' Ind_Ed_Normalized];
            [nrr,ncc]=size(sortederr_Ind);
            End_sp11=NaN*ones(1,ncc);
             Estar22=cell(1,ncc);
             Estar22(:)={'******'};
            
            % End_sp11=NaN*ones(1,(length(sortederr_Ind)));
            %  Estar22=cell(1,(length(sortederr_Ind)));
            %  Estar22(:)={'******'};
            
            
            t_ed_Ind=1:ncc-5;
            Ind_Ed_Normalized_all=[Ind_Ed_Normalized_all;Column_headers_inval_ed num2cell(t_ed_Ind);sortederr_Ind(:, 1:5) num2cell(Ind_Ed_Normalized);num2cell(End_sp11);Estar22;num2cell(End_sp11)];
            
            end
            
            
            
            N_active_array_ed(igp_eduction)=sum((cell2mat(sortederr_eduction(:, 4))));
            
            %% Here is the problem
            %sortederr_eduction_Avg=nanmean(cell2mat(sortederr_eduction (:, 4:end)),1);
            
            %sortederr_eduction_Avg=sum(cell2mat(sortederr_eduction (:, 4:end)),1)/N_active(igp_eduction);
            sortederr_eduction_Avg=sum(cell2mat(sortederr_eduction (:, 5:end)),1)/N_active_array_ed(igp_eduction);
            
            Normalized_eduction_bar=sortederr_eduction_Avg/nanmean(sortederr_eduction_Avg,2);
            
            
            sortederr_Ndays_Avg=nanmean(cell2mat(sortederr_Ndays (:, 4:end)),1);
            Normalized_Ndays=sortederr_Ndays_Avg/nanmean(sortederr_Ndays_Avg,2);
            
            
            
            Eduction_bar=[Eduction_bar;Normalized_eduction_bar];
            Eduction_Ndays=[Eduction_Ndays;Normalized_Ndays];
            
            
            Eduction_Day_by_Day=zeros(48,N_days_Eduction);
            
            
            for ig_av=1:N_days_Eduction;
                Eduction_Day_by_Day(:,ig_av)=Normalized_Ndays((ig_av-1)*48+1:(ig_av)*48);% 30 bin for 24 hr
            end
            
            
            
            std_error=nanstd(Eduction_Day_by_Day,0,2)/sqrt(N_days_Eduction);
            std_error_all=[std_error_all std_error];
            
            
            
            GT_New{igp_eduction}=(err_eduction(1, 1));
            New_GT=GT_New{igp_eduction};
            
            
            %ed_result_avg_std=[Title_avg Title_std;num2cell(nanmean(Day_by_Day_Edution_bin,2)) num2cell(std_error)];%Average eduction
            Title_avg=[New_GT{1}, '   ' 'Average'];
            Title_std=[New_GT{1},  '   ' 'SE'];
            
            ed_result_avg_std=[{Title_avg} {Title_std};num2cell(nanmean(Eduction_Day_by_Day,2)) num2cell(std_error)];%Average eduction
            ed_result_avg_std_all=[ed_result_avg_std_all ed_result_avg_std];
            
            
            
            end
            
            end
            
            
            
            try
                Triage_all_anti_ed=Result_all_anti(:,1:5);
            catch
                try
                Triage_all_anti_ed=Result_all_ed_TW;
                end
            end
            
            
            
            try
            Triage_results=[Result_all(:,1:5), Triage_all_anti_ed(:,5), Periodo_all(:,5)];
            catch
                try
                    Triage_results=[Result_all(:,1:5), Triage_all_anti_ed(:,5)];
                catch
                    try
                      Triage_results=[Result_all(:,1:5), Periodo_all(:,5)];
                    catch
                        try
                        Triage_results=[Triage_all_anti_ed(:,1:5), Periodo_all(:,5)];
                        catch
                            try
                                Triage_results=[Result_all(:,1:5)];
                            catch
                                try
                                    Triage_results=[Triage_all_anti_ed(:,1:5)];
                                catch
                                    try
                                    Triage_results=[Periodo_all(:,1:5)];
                                    catch
                                        disp("Please select atleast one analysis")
                                        set(handles.M_box,'String',[oldmsgs;{'Please select atleast one analysis'}] );drawnow
            
                                        return
                    end
                end
                        end
                    end
                end
                end
            end
            
            
            
            disp('Please wait...')
            
            set(handles.M_box,'String',[oldmsgs;{'Please wait...'}] );drawnow
            
            
            
            % set(handles.M_box, 'string', 'Please wait...')
            %f = msgbox('Please wait...');
            
            
             switch periodo_fig_select2
            
                 case 'Yes'
            %% Combaining periodogram data together
            
            %try
            [~,~,X_periodo]=unique( periodo_data_all(1,:));
            
            C_periodo = accumarray(X_periodo,1:size(periodo_data_all,2),[],@(r){periodo_data_all(:,r)});
            
            [nr_p,nc_p]=size(C_periodo);
            
            A_periodo=[];
            for i=1:nr_p
            
                avg=mean(cell2mat(C_periodo{i}(2:end,:)),2);
               avg1=[(C_periodo{i}(1,1));num2cell(avg)];
                A_periodo=[A_periodo avg1];
            
            end
            
            
            
            
            
             N_bin=round(length(A_periodo)/30);
            
             ydata=[];
             for i_bin=1:N_bin
                 ydata(i_bin,:)=sum(cell2mat(A_periodo((i_bin-1)*30+2:(i_bin*30)+1,:)));% for 30 min bin 0.5hr resulution
             end
            
            
            
            
            %Edution_bin_Nflies(i_bin,:)=sum(Data_Eduction_Nflies((i_bin-1)*30+1:(i_bin*30),:));
            
            
            [nr_pp,nc_pp]=size(A_periodo);
            for i_pp=1:nc_pp
            
            % activity_data=cell2mat(A_periodo(2:end,i_pp));
            % Time_data=1:1:length(activity_data);
            
            activity_data=ydata(1:end,i_pp);
            Time_data=1:30:length(A_periodo)-1;
            
            
            nop=5;
            figure;
            
            % periodogram calculations
            
             [period_1, power1, sigline1]=chi2_periodogram(app, Time_data,activity_data,nop,100,'perms',[3], A_periodo(1,i_pp), alpha_val, per_min, per_max);
            
            Title_power=[ A_periodo{1,i_pp}, '   ' 'Qp'];
            Title_sig=[ A_periodo{1,i_pp},  '   ' 'Sig.line'];
            
            periodogram_ploting_data_1=[{Title_power} {Title_sig}; num2cell(power1') num2cell(sigline1')];
            periodogram_ploting_data_2=[periodogram_ploting_data_2 periodogram_ploting_data_1];
            periodogram_ploting_data_3=[{'period'};num2cell(period_1')];
            periodogram_ploting_data=[periodogram_ploting_data_3 periodogram_ploting_data_2];
            
            end
            
            
             %% saving periodogram
            
            path1 = handles.path;   % mention your path
            myfolder1 = 'periodogram figures' ;   % new folder name
            % folder = mkdir([path1,filesep,myfolder1]) ;
            % path1  = [path1,filesep,myfolder1] ;
            
            
            folder1 = mkdir([path1,Project_name,filesep,myfolder1]);
            path1  = [path1,Project_name,filesep,myfolder1] ;
            
            
             for k = 1:nr_p
            
                   GT_name= A_periodo(1,k);
            
                    figure(k);
                    temp1=[path1,filesep,GT_name{1},'.png'];
                    try
                    saveas(gca,temp1);
                   close(figure(k))
                    catch
                        disp("Error!! Please check genotype name. Don't use backslash '\' or '>'  for genotype naming if you need periodogram figures")
                                        set(handles.M_box,'String',[oldmsgs;{"Error!! Please check genotype name.  Don't use backslash '\' or '>' for genotype naming if you need periodogram figures"}] );drawnow
            
                                        return
                    end
            
             end
            
            
             end
            
            
            Analysis_info1=['L_ON =' num2str(t_Lon)];
            Analysis_info2=['L_Off =' num2str(t_Loff)];
            Analysis_info3=['Start day of sleep =' num2str(start_day_sleep)];
            Analysis_info4=['Number of days of sleep analysis =' num2str(N_days_analysis)];
            
            Analysis_info5=['Start day of SD =' num2str(SD_start_day)];
            Analysis_info6=['Start time of SD =' num2str(SD_start_hr)];
            Analysis_info7=['SD duration =' num2str(SD_duration)];
            Analysis_info8=['SL interval =' num2str(SD_interval)];
            Analysis_info81=['Rebound interval =' num2str(SG_interval)];
            
            
            Analysis_info9=['Start day of eduction/anticipation =' num2str(Eduction_startday)];
            Analysis_info10=['Number of days of eduction/anticipation =' num2str(N_days_Eduction)];
            
            Analysis_info11=['Sleep inteval =' num2str(N_inval_hr)];
            Analysis_info12=['Periodogram significant level =' num2str(Significance_level)];
            Analysis_info13=['Minimum period =' num2str(per_min)];
            Analysis_info14=['Maximum period =' num2str(per_max)];
            Analysis_info15=['Start day of periodogram =' num2str(st_day_per)];
            Analysis_info16=['Number of days of periodogram =' num2str(N_day_per)];
            Analysis_info17=['Start day of actogram =' num2str(start_day_actogram)];
            Analysis_info18=['Number of days of actogram =' num2str(N_days_actogram)];
            
            
            
            TXT_ana=TXT1(:,1:6);
            Ana_sp=NaN*ones(1,5);
            Analysis_title=['Analysis parameters'];
            Analysis_g=['Genotype specification information'];
            %  Analysis_info=[{Analysis_title};{Analysis_info1};{Analysis_info2};{Analysis_info3};{Analysis_info4};{Analysis_info5};{Analysis_info6};...
            %      {Analysis_info7};{Analysis_info8};{Analysis_info9};{Analysis_info10};{Analysis_info11};{Analysis_info12};{Analysis_info13};{Analysis_info14}];
            %
            %
            Analysis_info=[{Analysis_g} num2cell(Ana_sp);TXT_ana;{Analysis_title} num2cell(Ana_sp);...
                {Analysis_info1} num2cell(Ana_sp);{Analysis_info2} num2cell(Ana_sp);{Analysis_info3} num2cell(Ana_sp);...
                {Analysis_info4} num2cell(Ana_sp);{Analysis_info5} num2cell(Ana_sp);{Analysis_info6} num2cell(Ana_sp);...
                 {Analysis_info7} num2cell(Ana_sp);{Analysis_info8} num2cell(Ana_sp);{Analysis_info81} num2cell(Ana_sp);...
            
                 {Analysis_info9} num2cell(Ana_sp);{Analysis_info10} num2cell(Ana_sp);{Analysis_info11} num2cell(Ana_sp);...
                 {Analysis_info12} num2cell(Ana_sp);{Analysis_info13} num2cell(Ana_sp);{Analysis_info14} num2cell(Ana_sp);...
                 {Analysis_info15} num2cell(Ana_sp);{Analysis_info16} num2cell(Ana_sp); {Analysis_info17} num2cell(Ana_sp);...
                 {Analysis_info18} num2cell(Ana_sp)];
            
            
            
            %% save results
            
            
            path1 = handles.path;   % mention your path
            myfolder1 = Project_name ;   % new folder name
            folder1 = mkdir([path1,filesep,myfolder1]) ;
            Result_folder  = [path1,filesep,myfolder1] ;
            
            save_location_analysis_par=fullfile(Result_folder, strcat(Project_name,'_Analysis_parameters.xls'));
            save_location_24=fullfile(Result_folder, strcat(Project_name,'_Sleep_analysis_results_24.xls'));
            save_location_summary_24=fullfile(Result_folder, strcat(Project_name,'_Sleep_analysis_summary.xls'));
            save_location_eduction=fullfile(Result_folder, strcat(Project_name,'_Eduction_result.xls'));
            save_location_Triage_info=fullfile(Result_folder, strcat(Project_name,'_Triage_info.xls'));
            
            save_location_L_D=fullfile(Result_folder, strcat(Project_name,'_L&D_Results.xls'));
            save_location_lifespan=fullfile(Result_folder, strcat(Project_name,'_Lifespan_date_time.xls'));
            save_location_lifespangp=fullfile(Result_folder, strcat(Project_name,'_lifespan.xls'));
            save_location_Eduction_Ind=fullfile(Result_folder, strcat(Project_name,'_Eduction_individual_flies.xls'));
            save_location_periodogram=fullfile(Result_folder, strcat(Project_name,'_Periodogram_results_individual_flies.xls'));
            save_location_per_ploting=fullfile(Result_folder, strcat(Project_name,'_Periodogram_data_for_custom_plot.xls'));
            
            
            
            save_location_Sleep_Int=fullfile(Result_folder, strcat(Project_name,'_Sleep_interval.xls'));
            save_location_BoutN_Int=fullfile(Result_folder, strcat(Project_name,'_BoutN_interval.xls'));
            save_location_BoutL_Int=fullfile(Result_folder, strcat(Project_name,'_BoutL_interval.xls'));
            save_location_Activity_Int=fullfile(Result_folder, strcat(Project_name,'_Activity_interval.xls'));
            save_location_anticipation=fullfile(Result_folder, strcat(Project_name,'_Anticipation_results.xls'));
            save_location_sleep=fullfile(Result_folder, strcat(Project_name,'_Sleep_24h.xls'));
            save_location_BoutN=fullfile(Result_folder, strcat(Project_name,'_Bout_Number_24h.xls'));
            save_location_BoutL=fullfile(Result_folder, strcat(Project_name,'_Average_Bout_Length_24h.xls'));
            save_location_Totalact=fullfile(Result_folder, strcat(Project_name,'_Total_activity_24h.xls'));
            save_location_Activity_count=fullfile(Result_folder, strcat(Project_name,'_Activity_per_waking_min_24h.xls'));
            save_location_Latency=fullfile(Result_folder, strcat(Project_name,'_Latency_24h.xls'));
            save_location_Manti=fullfile(Result_folder, strcat(Project_name,'_Morning_anticipation_24h.xls'));
            save_location_Eanti=fullfile(Result_folder, strcat(Project_name,'_Evening_anticipation_24h.xls'));
            save_location_Manti_ph=fullfile(Result_folder, strcat(Project_name,'_Morning_anticipation_phase_24h.xls'));
            save_location_Eanti_ph=fullfile(Result_folder, strcat(Project_name,'_Evening_anticipation_phase_24h.xls'));
            save_location_sleep_lost=fullfile(Result_folder, strcat(Project_name,'_Sleep_lost_24h.xls'));
            save_location_per_sleep_lost=fullfile(Result_folder, strcat(Project_name,'_Percentage_Sleep_lost_24h.xls'));
            save_location_sleep_gain=fullfile(Result_folder, strcat(Project_name,'_Sleep_gain_24h.xls'));
            save_location_per_sleep_gain=fullfile(Result_folder, strcat(Project_name,'_Percentage_Sleep_gain_24h.xls'));
            save_location_per_sleep_gain=fullfile(Result_folder, strcat(Project_name,'_Percentage_Sleep_gain_24h.xls'));
            save_location_Latency_SD=fullfile(Result_folder, strcat(Project_name,'_Latency_SD_24h.xls'));
            
            % save_location_SC=fullfile(Result_folder, strcat(Project_name,'_Results_for_screening.xls'));
            
            
            
            
            %% Initialisation of POI Libs
            % Add Java POI Libs to matlab javapath
            
            % javaaddpath('poi_library/poi-3.8-20120326.jar');
            % javaaddpath('poi_library/poi-ooxml-3.8-20120326.jar');
            % javaaddpath('poi_library/poi-ooxml-schemas-3.8-20120326.jar');
            % javaaddpath('poi_library/xmlbeans-2.3.0.jar');
            % javaaddpath('poi_library/dom4j-1.6.1.jar');
            % javaaddpath('poi_library/stax-api-1.0.1.jar');
            
            
            %xlswrite(save_location_analysis_par,Analysis_info);
            writecell(Analysis_info,save_location_analysis_par);
            writecell(Triage_results,save_location_Triage_info);
            
            try
                writecell(Result_all,save_location_24);
                %writetable(Result_all_SC,save_location_SC);
                writecell(summary_all,save_location_summary_24);
                writecell(Result_all_L_D,save_location_L_D);
            end
            
            
            try
            writecell(Triage_all,save_location_lifespan);
            writecell(Life_gp,save_location_lifespangp);
            end
            
            try
            writecell(ed_result_avg_std_all,save_location_eduction);
            end
            
            try
                writecell(Result_all_anti,save_location_anticipation);
            end
            
            try
                writecell(Ind_Ed_Normalized_all,save_location_Eduction_Ind);
            end
            
            
            try
                writecell(Periodo_all,save_location_periodogram);
            end
            
            try
            writecell(periodogram_ploting_data,save_location_per_ploting);
            end
            %*********************************************************
            
            %per interval results
            try
            writecell(Sleep_interval,save_location_Sleep_Int);
            writecell(BoutN_interval,save_location_BoutN_Int);
            writecell(BoutL_interval,save_location_BoutL_Int);
            writecell(Activity_interval,save_location_Activity_Int);
            end
            
            
            %*************************************************************
            
            
            
            if strcmp(Outfile_select2, 'Yes')
            
            try
            writecell(sleep_gp_final,save_location_sleep);
            writecell(Bout_N_gp_final,save_location_BoutN);
            writecell(Bout_L_gp_final,save_location_BoutL);
            writecell(Total_act_gp_final,save_location_Totalact);
            writecell(act_count_gp_final,save_location_Activity_count);
            writecell(Latency_gp_final,save_location_Latency);
            
            end
            
            
            try
            writecell(sleep_lost_gp_final,save_location_sleep_lost);
            writecell(Per_sleep_lost_gp_final,save_location_per_sleep_lost);
            writecell(sleep_gain_gp_final,save_location_sleep_gain);
            writecell(Per_sleep_gain_gp_final,save_location_per_sleep_gain);
            writecell(Latency_SD_gp_final,save_location_Latency_SD);
            end
            
            %% Write Result/24hr
            
            try
            writecell(Manti_gp_final,save_location_Manti);
            writecell(Eanti_gp_final,save_location_Eanti);
            writecell(Manti_ph_gp_final,save_location_Manti_ph);
            writecell(Eanti_ph_gp_final,save_location_Eanti_ph);
            end
            
            end
            
            %% ************************************************************ ************ *********************************************************************************************************
                                                                                         %% Eduction Plot start%%
            %% ***********************************************************************************************************************************************************************************
            
            if strcmp(Eduction_select2, 'Yes')
            
            ed_result=Eduction_bar';
            
            ed_result_Ndays=Eduction_Ndays';
            
             tp1=1:1:24-Photo_period;
            tp2=(24-Photo_period)+1:1:(24-Photo_period)+(2*Photo_period);
            tp3=(24-Photo_period)+(2*Photo_period)+1:1:48;
            
            
            bin_1=ed_result(1:24-Photo_period,:);
            
            
            bin_2=ed_result((24-Photo_period)+1:(24-Photo_period)+(2*Photo_period),:);
            %Normalized_err2=err_result((24-Photo_period)+1:(24-Photo_period)+(2*Photo_period),:);
            
            bin_3=ed_result((24-Photo_period)+(2*Photo_period)+1:48,:);
            %Normalized_err3=err_result((24-Photo_period)+(2*Photo_period)+1:48,:);
            
            [npr,npc]=size(ed_result);
            
            for i_plot=1:npc
                Normalized_err1=std_error_all(1:24-Photo_period,i_plot);% upper limit for error bar
                L1=zeros(size(Normalized_err1));% Lower limit
                Normalized_err2=std_error_all((24-Photo_period)+1:(24-Photo_period)+(2*Photo_period),i_plot);
                L2=zeros(size(Normalized_err2));% Lower limit
                Normalized_err3=std_error_all((24-Photo_period)+(2*Photo_period)+1:48,i_plot);
                L3=zeros(size(Normalized_err3));% Lower limit
            
            figure;
            h1=bar(tp1,bin_1(:,i_plot),1);
            set(h1,'FaceColor','k')
              hold on
            er1 = errorbar(tp1,bin_1(:,i_plot),L1,Normalized_err1,'k');
            set(er1,'linestyle','none');
             hold on
            h2=bar(tp2,bin_2(:,i_plot),1);
             set(h2,'FaceColor','w')
             hold on
            er2 = errorbar(tp2,bin_2(:,i_plot),L2,Normalized_err2,'k');
            set(er2,'linestyle','none');
             hold on
            h3=bar(tp3,bin_3(:,i_plot),1);
             set(h3,'FaceColor','k')
             hold on
            er3 = errorbar(tp3,bin_3(:,i_plot),L3, Normalized_err3,'k');
            set(er3,'linestyle','none');
             set(gca, 'xtick', 2:2:48);
            
            New_GT=GT_New{i_plot};
            
               Title_plot=[New_GT{1}, '  '  'N= ' num2str(N_active_array_ed(i_plot))];
              title(Title_plot)
            
               ylabel('Normalized activity / 30 min')
              xlabel('30 min bin for 24 hr')
            
            
              figure;
              plot(ed_result_Ndays(:,i_plot))
              %title(geno_type1{i_plot}, num2str(N_active(i_plot)))
              title(Title_plot)
              ylabel('Normalized activity / 30 min')
              xlabel('Time( hr) ')
            
            end
            
            
            path = handles.path;   % mention your path
            myfolder = 'Eduction figures' ;   % new folder name
            % folder = mkdir([path,filesep,myfolder]) ;
            %path  = [path,filesep,myfolder] ;
            folder = mkdir([path,Project_name,filesep,myfolder]);
            path  = [path,Project_name,filesep,myfolder] ;
            
            
            n_gt=2*length(GT_New) ;
            plot_name={'-plot1', '-plot2'};
             for k = 1:n_gt
                 i_round=round(k/2);
                 ip=mod(k,2)+1;
                   GT_name=GT_New{i_round};
                   p_name=plot_name{ip};
                    figure(k);
                    temp=[path,filesep,GT_name{1},p_name,'.png'];
                            try
                    saveas(gca,temp);
                    close( figure(k) );
                    catch
                        disp("Error!! Please check genotype name. Don't use backslash '\' or '>'  for genotype naming if you need figures")
                                        set(handles.M_box,'String',[oldmsgs;{"Error!! Please check genotype name.  Don't use backslash '\' or '>' for genotype naming if you need figures"}] );drawnow
            
                                        return
                    end
             end
            
            
            
            end
            
            %% ************************************************************ ************ *********************************************************************************************************
                                                                                         %% Eduction Plot End %%
            %% ***********************************************************************************************************************************************************************************
            
            
            
            
            disp('Analysis COMPLETED! Results exported to xls file')
            set(handles.M_box,'String',[oldmsgs;{'Analysis COMPLETED!'}] );drawnow
            
            catch
                disp('Unknown Error! Please check the input(eg: Missing toolbox, genotype specification file, Monitor file, days...)')
            set(handles.M_box,'String',[oldmsgs;{'Unknown Error! Please check the input(eg: Missing toolbox, genotype specification file, days...)'}] );drawnow
            return;
            end
        end

        % Selection changed function: Triage_date
        function Triage_date_SelectionChangedFcn(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            
            
            T_date_select1=get(handles.Triage_date, 'SelectedObject');
            T_date_select2=get(T_date_select1, 'String');
        end

        % Selection changed function: Triage_date_ed
        function Triage_date_ed_SelectionChangedFcn(app, event)
            % --- Executes when selected object is changed in Triage_date_ed.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            
            T_date_ed_select1=get(handles.Triage_date_ed, 'SelectedObject');
            T_date_ed_select2=get(T_date_ed_select1, 'String');
        end

        % Selection changed function: Triage_date_per
        function Triage_date_per_SelectionChangedFcn(app, event)
            % --- Executes when selected object is changed in Triage_date_per.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            
            T_date_per_select1=get(handles.Triage_date_per, 'SelectedObject');
            T_date_per_select2=get(T_date_per_select1, 'String');
        end

        % Selection changed function: analysis_per_interval
        function analysis_per_interval_SelectionChangeFcn(app, event)
            % --- Executes when selected object is changed in analysis_per_interval.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            
            Interval_select1=get(handles.analysis_per_interval, 'SelectedObject');
            Interval_select2=get(Interval_select1, 'String');
            
            
            switch Interval_select2
                case 'Yes'
                   set(handles.per_Interval, 'Enable', 'on');
                case 'No'
                   set(handles.per_Interval, 'Enable', 'off');
             end
        end

        % Selection changed function: cum_fig
        function cum_fig_SelectionChangedFcn(app, event)
            % --- Executes when selected object is changed in cum_fig.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            
            cum_fig_select1=get(handles.cum_fig, 'SelectedObject');
            cum_fig_select2=get(cum_fig_select1, 'String');
        end

        % Selection changed function: lifespan
        function lifespan_SelectionChangedFcn(app, event)
            % --- Executes when selected object is changed in lifespan.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            
            
            
            
            Lifespan_select1=get(handles.lifespan, 'SelectedObject');
            Lifespan_select2=get(Lifespan_select1, 'String');
            
            
            switch Lifespan_select2
            
                case 'Yes'
            set(handles.LF_th, 'Enable', 'on');
            set(handles.Max_LF, 'Enable', 'on');
            
               case 'No'
             set(handles.LF_th, 'Enable', 'off');
             set(handles.Max_LF, 'Enable', 'off');
            
            end
        end

        % Button pushed function: open_file
        function open_file_Callback(app, event)
            % --- Executes on button press in open_file.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            
            warning off
            
            [file_name path_name]=uigetfile({'*.xlsx'},'File Selector');
            Fullpathname=strcat(path_name, file_name);
            handles.Fullpath = Fullpathname;
            set(handles.File_Name,'String',Fullpathname);
            
            handles.path=path_name;
            
            
            
            
            
            
            set(handles.t_LON, 'Enable', 'on');
            set(handles.t_Loff1, 'Enable', 'on');
            set(handles.Start, 'Enable', 'on');
            
            
            set(handles.N_days_ana, 'Enable', 'on');
            set(handles.st_day_sleep, 'Enable', 'on');
            
            
            set(handles.Int_No, 'Enable', 'on');
            set(handles.Int_Yes, 'Enable', 'on');
            
            set(handles.SD_No, 'Enable', 'on');
            set(handles.SD_Yes, 'Enable', 'on');
            
            
            
            
            
                   set(handles.Threshold2, 'Enable', 'on');
                   set(handles.Threshold3, 'Enable', 'on');
                   set(handles.Threshold4, 'Enable', 'on');
            
            
            set(handles.T_date_No, 'Enable', 'on');
            set(handles.T_date_Yes, 'Enable', 'on');
            
            
            
            set(handles.Ed_No, 'Enable', 'on');
            set(handles.Ed_Yes, 'Enable', 'on');
            
            set(handles.periodo_No, 'Enable', 'on');
            set(handles.periodo_Yes, 'Enable', 'on');
            
            
            
            set(handles.project_name, 'Enable', 'on');
            
            
            
            set(handles.SL_No, 'Enable', 'on');
            set(handles.SL_Yes, 'Enable', 'on');
            
            set(handles.T_date_ed_No, 'Enable', 'on');
             set(handles.T_date_ed_Yes, 'Enable', 'on');
            
            
                   set(handles.st_ed_anti, 'Enable', 'on');
                   set(handles.N_ed_anti, 'Enable', 'on');
            %        set(handles.Ind_Ed_No, 'Enable', 'on');
            %        set(handles.Ind_Ed_Yes, 'Enable', 'on');
                    set(handles.Anti_Yes, 'Enable', 'on');
                   set(handles.Anti_No, 'Enable', 'on');
                   set(handles.Eduction_Yes, 'Enable', 'on');
                   set(handles.Eduction_No, 'Enable', 'on');
            
            
            
                     set(handles.act_Yes, 'Enable', 'on');
                     set(handles.act_No, 'Enable', 'on');
            
            
              set(handles.LF_Yes, 'Enable', 'on');
             set(handles.LF_No, 'Enable', 'on');
            
            
             set(handles.Out_Yes, 'Enable', 'on');
             set(handles.Out_No, 'Enable', 'on');
            
            
            set(handles.LF_th, 'Enable', 'off');
            set(handles.Max_LF, 'Enable', 'off');
            
                   set(handles.SL_ch1, 'Enable', 'on');
                   set(handles.SL_ch2, 'Enable', 'on');
                   set(handles.SL_ch3, 'Enable', 'on');
            
            
                    set(handles.ed_ch1, 'Enable', 'on');
                   set(handles.ed_ch2, 'Enable', 'on');
                   set(handles.ed_ch3, 'Enable', 'on');
            
                      set(handles.Threshold3_ed, 'Enable', 'on');
                     set(handles.Threshold2_ed, 'Enable', 'on');
                     set(handles.Threshold4_ed, 'Enable', 'on');
            
            
            guidata(hObject, handles);
        end

        % Selection changed function: periodogram
        function periodogram_SelectionChangedFcn(app, event)
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            
            
            periodo_select1=get(handles.periodogram, 'SelectedObject');
            periodo_select2=get(periodo_select1, 'String');
            
            
            
            
            switch periodo_select2
                case 'Yes'
                   set(handles.sig_val, 'Enable', 'on');
                   set(handles.start_per, 'Enable', 'on');
                   set(handles.N_per, 'Enable', 'on');
                   set(handles.min_period, 'Enable', 'on');
                   set(handles.max_period, 'Enable', 'on');
            
                   set(handles.per_ch1, 'Enable', 'on');
                   set(handles.per_ch2, 'Enable', 'on');
                   set(handles.per_ch3, 'Enable', 'on');
            
            
            
            
            
                     set(handles.T_date_No_per, 'Enable', 'on');
                     set(handles.T_date_Yes_per, 'Enable', 'on');
                     set(handles.Threshold3_per, 'Enable', 'on');
                      set(handles.Threshold4_per, 'Enable', 'on');
                       set(handles.Threshold2_per, 'Enable', 'on');
                      set(handles.PF_Yes, 'Enable', 'on');
                      set(handles.PF_No, 'Enable', 'on');
            
                case 'No'
                   set(handles.sig_val, 'Enable', 'off');
                   set(handles.start_per, 'Enable', 'off');
                   set(handles.N_per, 'Enable', 'off');
                   set(handles.min_period, 'Enable', 'off');
                   set(handles.max_period, 'Enable', 'off');
            
            
                     %set(handles.Triage_text3, 'Enable', 'off');
            
            
                     set(handles.T_date_No_per, 'Enable', 'off');
                     set(handles.T_date_Yes_per, 'Enable', 'off');
                     set(handles.Threshold3_per, 'Enable', 'off');
                      set(handles.Threshold4_per, 'Enable', 'off');
                       set(handles.Threshold2_per, 'Enable', 'off');
            
                       set(handles.PF_Yes, 'Enable', 'off');
                      set(handles.PF_No, 'Enable', 'off');
            
                   set(handles.per_ch1, 'Enable', 'off');
                   set(handles.per_ch2, 'Enable', 'off');
                   set(handles.per_ch3, 'Enable', 'off');
            
             end
        end

        % Selection changed function: sleep_ana
        function sleep_ana_SelectionChangedFcn(app, event)
            % --- Executes when selected object is changed in sleep_ana.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            
            
            SL_select1=get(handles.sleep_ana, 'SelectedObject');
            SL_select2=get(SL_select1, 'String');
            
            
            switch SL_select2
            
                case 'Yes'
            
             set(handles.N_days_ana, 'Enable', 'on');
            set(handles.st_day_sleep, 'Enable', 'on');
            
            
            set(handles.Int_No, 'Enable', 'on');
            set(handles.Int_Yes, 'Enable', 'on');
            
            set(handles.SD_No, 'Enable', 'on');
            set(handles.SD_Yes, 'Enable', 'on');
            
            
            
                   set(handles.Threshold2, 'Enable', 'on');
                   set(handles.Threshold3, 'Enable', 'on');
                   set(handles.Threshold4, 'Enable', 'on');
            
            set(handles.T_date_No, 'Enable', 'on');
            set(handles.T_date_Yes, 'Enable', 'on');
            
                   set(handles.SL_ch1, 'Enable', 'on');
                   set(handles.SL_ch2, 'Enable', 'on');
                   set(handles.SL_ch3, 'Enable', 'on');
            
            case 'No'
            
            set(handles.N_days_ana, 'Enable', 'off');
            set(handles.st_day_sleep, 'Enable', 'off');
            
            
            set(handles.Int_No, 'Enable', 'off');
            set(handles.Int_Yes, 'Enable', 'off');
            
            set(handles.SD_No, 'Enable', 'off');
            set(handles.SD_Yes, 'Enable', 'off');
            
            
                   set(handles.Threshold2, 'Enable', 'off');
                   set(handles.Threshold3, 'Enable', 'off');
                   set(handles.Threshold4, 'Enable', 'off');
            
                   set(handles.SL_ch1, 'Enable', 'off');
                   set(handles.SL_ch2, 'Enable', 'off');
                   set(handles.SL_ch3, 'Enable', 'off');
            
            set(handles.T_date_No, 'Enable', 'off');
            set(handles.T_date_Yes, 'Enable', 'off');
            
            
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create figure1 and hide until all components are created
            app.figure1 = uifigure('Visible', 'off');
            colormap(app.figure1, 'parula');
            app.figure1.Position = [555 275 1037 587];
            app.figure1.Name = 'sleepmat';
            app.figure1.Resize = 'off';
            app.figure1.HandleVisibility = 'callback';
            app.figure1.Tag = 'figure1';

            % Create Start
            app.Start = uibutton(app.figure1, 'push');
            app.Start.ButtonPushedFcn = createCallbackFcn(app, @Start_Callback, true);
            app.Start.Tag = 'Start';
            app.Start.BackgroundColor = [0 1 0];
            app.Start.FontSize = 10.6666666666667;
            app.Start.Enable = 'off';
            app.Start.Position = [56.3333333333333 128.333333333333 112 25.3333333333333];
            app.Start.Text = 'Start';

            % Create uipanel14
            app.uipanel14 = uipanel(app.figure1);
            app.uipanel14.Title = 'General input';
            app.uipanel14.BackgroundColor = [0 1 1];
            app.uipanel14.Tag = 'uipanel14';
            app.uipanel14.FontSize = 21.3333333333333;
            app.uipanel14.Position = [4 163 230 423];

            % Create open_file
            app.open_file = uibutton(app.uipanel14, 'push');
            app.open_file.ButtonPushedFcn = createCallbackFcn(app, @open_file_Callback, true);
            app.open_file.Tag = 'open_file';
            app.open_file.BackgroundColor = [0.0588235294117647 1 1];
            app.open_file.FontSize = 13.3333333333333;
            app.open_file.Position = [6 346 68 31.3333333333333];
            app.open_file.Text = 'Open file';

            % Create File_Name
            app.File_Name = uieditfield(app.uipanel14, 'text');
            app.File_Name.Tag = 'File_Name';
            app.File_Name.HorizontalAlignment = 'center';
            app.File_Name.FontSize = 10.6666666666667;
            app.File_Name.BackgroundColor = [0.0588235294117647 1 1];
            app.File_Name.Enable = 'off';
            app.File_Name.Position = [74 347 146.666666666667 30.6666666666667];

            % Create text4
            app.text4 = uilabel(app.uipanel14);
            app.text4.Tag = 'text4';
            app.text4.BackgroundColor = [0.0588235294117647 1 1];
            app.text4.HorizontalAlignment = 'center';
            app.text4.VerticalAlignment = 'top';
            app.text4.WordWrap = 'on';
            app.text4.FontSize = 10.6666666666667;
            app.text4.Position = [5 246 187.333333333333 16.6666666666666];
            app.text4.Text = 'Light on time (Clock time, 24 hr format)';

            % Create t_LON
            app.t_LON = uieditfield(app.uipanel14, 'text');
            app.t_LON.Tag = 't_LON';
            app.t_LON.HorizontalAlignment = 'center';
            app.t_LON.FontSize = 10.6666666666667;
            app.t_LON.Enable = 'off';
            app.t_LON.Position = [196 241 24 30];
            app.t_LON.Value = '8';

            % Create text5
            app.text5 = uilabel(app.uipanel14);
            app.text5.Tag = 'text5';
            app.text5.BackgroundColor = [0.0588235294117647 1 1];
            app.text5.HorizontalAlignment = 'center';
            app.text5.VerticalAlignment = 'top';
            app.text5.WordWrap = 'on';
            app.text5.FontSize = 10.6666666666667;
            app.text5.Position = [5 212 187.333333333333 19.3333333333333];
            app.text5.Text = 'Light off time (Clock time, 24 hr format)';

            % Create t_Loff1
            app.t_Loff1 = uieditfield(app.uipanel14, 'text');
            app.t_Loff1.Tag = 't_Loff1';
            app.t_Loff1.HorizontalAlignment = 'center';
            app.t_Loff1.FontSize = 10.6666666666667;
            app.t_Loff1.Enable = 'off';
            app.t_Loff1.Position = [197 207 24.6666666666667 30];
            app.t_Loff1.Value = '20';

            % Create P_name
            app.P_name = uilabel(app.uipanel14);
            app.P_name.Tag = 'P_name';
            app.P_name.BackgroundColor = [0.0588235294117647 1 1];
            app.P_name.HorizontalAlignment = 'center';
            app.P_name.VerticalAlignment = 'top';
            app.P_name.WordWrap = 'on';
            app.P_name.FontSize = 13.3333333333333;
            app.P_name.Position = [5 277 84 28.6666666666667];
            app.P_name.Text = 'Project Name';

            % Create project_name
            app.project_name = uieditfield(app.uipanel14, 'text');
            app.project_name.Tag = 'project_name';
            app.project_name.HorizontalAlignment = 'center';
            app.project_name.FontSize = 10.6666666666667;
            app.project_name.Enable = 'off';
            app.project_name.Position = [90 280 126.666666666667 30];
            app.project_name.Value = 'project-1';

            % Create text35
            app.text35 = uilabel(app.uipanel14);
            app.text35.Tag = 'text35';
            app.text35.BackgroundColor = [0.0588235294117647 1 1];
            app.text35.HorizontalAlignment = 'center';
            app.text35.VerticalAlignment = 'top';
            app.text35.WordWrap = 'on';
            app.text35.FontSize = 10.6666666666667;
            app.text35.Position = [20 317 187.333333333333 16.6666666666667];
            app.text35.Text = 'Upload a geneotype specification file';

            % Create lifespan_th
            app.lifespan_th = uilabel(app.uipanel14);
            app.lifespan_th.Tag = 'lifespan_th';
            app.lifespan_th.BackgroundColor = [0.0588235294117647 1 1];
            app.lifespan_th.HorizontalAlignment = 'center';
            app.lifespan_th.VerticalAlignment = 'top';
            app.lifespan_th.WordWrap = 'on';
            app.lifespan_th.FontSize = 10.6666666666667;
            app.lifespan_th.Position = [21 116 98.6666666666667 16.6666666666667];
            app.lifespan_th.Text = 'Activity threshold';

            % Create text46
            app.text46 = uilabel(app.uipanel14);
            app.text46.Tag = 'text46';
            app.text46.BackgroundColor = [0.0588235294117647 1 1];
            app.text46.HorizontalAlignment = 'center';
            app.text46.VerticalAlignment = 'top';
            app.text46.WordWrap = 'on';
            app.text46.FontSize = 10.6666666666667;
            app.text46.Position = [133 112 84.6666666666667 31.3333333333333];
            app.text46.Text = 'Maximum days of lifespan';

            % Create lifespan
            app.lifespan = uibuttongroup(app.uipanel14);
            app.lifespan.SelectionChangedFcn = createCallbackFcn(app, @lifespan_SelectionChangedFcn, true);
            app.lifespan.Title = 'Lifespan';
            app.lifespan.BackgroundColor = [0 1 1];
            app.lifespan.Tag = 'lifespan';
            app.lifespan.FontSize = 13.3333333333333;
            app.lifespan.Position = [22 152 153 47];

            % Create LF_Yes
            app.LF_Yes = uiradiobutton(app.lifespan);
            app.LF_Yes.Tag = 'LF_Yes';
            app.LF_Yes.Enable = 'off';
            app.LF_Yes.Text = 'Yes';
            app.LF_Yes.FontSize = 10.6666666666667;
            app.LF_Yes.Position = [12 0 54.6666666666667 24.6666666666667];

            % Create LF_No
            app.LF_No = uiradiobutton(app.lifespan);
            app.LF_No.Tag = 'LF_No';
            app.LF_No.Enable = 'off';
            app.LF_No.Text = 'No';
            app.LF_No.FontSize = 10.6666666666667;
            app.LF_No.Position = [88 5 54.6666666666667 24.6666666666667];
            app.LF_No.Value = true;

            % Create LF_th
            app.LF_th = uieditfield(app.uipanel14, 'text');
            app.LF_th.Tag = 'LF_th';
            app.LF_th.HorizontalAlignment = 'center';
            app.LF_th.FontSize = 10.6666666666667;
            app.LF_th.Enable = 'off';
            app.LF_th.Position = [49 83 28 27.3333333333333];
            app.LF_th.Value = '2';

            % Create Max_LF
            app.Max_LF = uieditfield(app.uipanel14, 'text');
            app.Max_LF.Tag = 'Max_LF';
            app.Max_LF.HorizontalAlignment = 'center';
            app.Max_LF.FontSize = 10.6666666666667;
            app.Max_LF.Enable = 'off';
            app.Max_LF.Position = [164 82 26.6666666666667 28];
            app.Max_LF.Value = '11';

            % Create Out_files
            app.Out_files = uibuttongroup(app.uipanel14);
            app.Out_files.SelectionChangedFcn = createCallbackFcn(app, @Out_files_SelectionChangedFcn, true);
            app.Out_files.Title = 'Individual output files';
            app.Out_files.BackgroundColor = [0 1 1];
            app.Out_files.Tag = 'Out_files';
            app.Out_files.FontSize = 13.3333333333333;
            app.Out_files.Position = [36 15 157 55];

            % Create Out_Yes
            app.Out_Yes = uiradiobutton(app.Out_files);
            app.Out_Yes.Tag = 'Out_Yes';
            app.Out_Yes.Enable = 'off';
            app.Out_Yes.Text = 'Yes';
            app.Out_Yes.FontSize = 10.6666666666667;
            app.Out_Yes.Position = [11 4 54.4 24.7692307692308];

            % Create Out_No
            app.Out_No = uiradiobutton(app.Out_files);
            app.Out_No.Tag = 'Out_No';
            app.Out_No.Enable = 'off';
            app.Out_No.Text = 'No';
            app.Out_No.FontSize = 10.6666666666667;
            app.Out_No.Position = [87 4 54.6666666666667 24.6666666666667];
            app.Out_No.Value = true;

            % Create uipanel16
            app.uipanel16 = uipanel(app.figure1);
            app.uipanel16.Title = 'Circadian analysis';
            app.uipanel16.BackgroundColor = [0.4706 0.6706 0.1882];
            app.uipanel16.Tag = 'uipanel16';
            app.uipanel16.FontSize = 21.3333333333333;
            app.uipanel16.Position = [547 0 487 587];

            % Create Ed_plot
            app.Ed_plot = uibuttongroup(app.uipanel16);
            app.Ed_plot.SelectionChangedFcn = createCallbackFcn(app, @Ed_plot_SelectionChangeFcn, true);
            app.Ed_plot.Title = 'Circadian analysis';
            app.Ed_plot.BackgroundColor = [0.466666666666667 0.674509803921569 0.188235294117647];
            app.Ed_plot.Tag = 'Ed_plot';
            app.Ed_plot.FontSize = 13.3333333333333;
            app.Ed_plot.Position = [149 509 163 47];

            % Create Ed_Yes
            app.Ed_Yes = uiradiobutton(app.Ed_plot);
            app.Ed_Yes.Tag = 'Ed_Yes';
            app.Ed_Yes.Enable = 'off';
            app.Ed_Yes.Text = 'Yes';
            app.Ed_Yes.FontSize = 10.6666666666667;
            app.Ed_Yes.Position = [22 4 54.6666666666667 24.6666666666667];
            app.Ed_Yes.Value = true;

            % Create Ed_No
            app.Ed_No = uiradiobutton(app.Ed_plot);
            app.Ed_No.Interruptible = 'off';
            app.Ed_No.Tag = 'Ed_No';
            app.Ed_No.Enable = 'off';
            app.Ed_No.Text = 'No';
            app.Ed_No.FontSize = 10.6666666666667;
            app.Ed_No.Position = [98 5 53.3333333333333 24.6666666666667];

            % Create Anticipation
            app.Anticipation = uibuttongroup(app.uipanel16);
            app.Anticipation.SelectionChangedFcn = createCallbackFcn(app, @Anticipation_SelectionChangedFcn, true);
            app.Anticipation.Title = 'Anticipation analysis';
            app.Anticipation.BackgroundColor = [0.466666666666667 0.674509803921569 0.188235294117647];
            app.Anticipation.Tag = 'Anticipation';
            app.Anticipation.FontSize = 10.6666666666667;
            app.Anticipation.Position = [10 455 117 43];

            % Create Anti_Yes
            app.Anti_Yes = uiradiobutton(app.Anticipation);
            app.Anti_Yes.Tag = 'Anti_Yes';
            app.Anti_Yes.Enable = 'off';
            app.Anti_Yes.Text = 'Yes';
            app.Anti_Yes.FontSize = 10.6666666666667;
            app.Anti_Yes.Position = [6 3 46 24.6666666666667];
            app.Anti_Yes.Value = true;

            % Create Anti_No
            app.Anti_No = uiradiobutton(app.Anticipation);
            app.Anti_No.Interruptible = 'off';
            app.Anti_No.Tag = 'Anti_No';
            app.Anti_No.Enable = 'off';
            app.Anti_No.Text = 'No';
            app.Anti_No.FontSize = 10.6666666666667;
            app.Anti_No.Position = [56 4 47.3333333333333 24.6666666666667];

            % Create Eduction
            app.Eduction = uibuttongroup(app.uipanel16);
            app.Eduction.SelectionChangedFcn = createCallbackFcn(app, @Eduction_SelectionChangedFcn, true);
            app.Eduction.Title = 'Eduction analysis';
            app.Eduction.BackgroundColor = [0.466666666666667 0.674509803921569 0.188235294117647];
            app.Eduction.Tag = 'Eduction';
            app.Eduction.FontSize = 10.6666666666667;
            app.Eduction.Position = [142 453 103 47];

            % Create Eduction_Yes
            app.Eduction_Yes = uiradiobutton(app.Eduction);
            app.Eduction_Yes.Tag = 'Eduction_Yes';
            app.Eduction_Yes.Enable = 'off';
            app.Eduction_Yes.Text = 'Yes';
            app.Eduction_Yes.FontSize = 10.6666666666667;
            app.Eduction_Yes.Position = [12 7 46.6666666666667 24.6666666666667];

            % Create Eduction_No
            app.Eduction_No = uiradiobutton(app.Eduction);
            app.Eduction_No.Tag = 'Eduction_No';
            app.Eduction_No.Enable = 'off';
            app.Eduction_No.Text = 'No';
            app.Eduction_No.FontSize = 10.6666666666667;
            app.Eduction_No.Position = [52 7 38.6666666666667 24.6666666666667];
            app.Eduction_No.Value = true;

            % Create st_ed_anti
            app.st_ed_anti = uieditfield(app.uipanel16, 'text');
            app.st_ed_anti.Tag = 'st_ed_anti';
            app.st_ed_anti.HorizontalAlignment = 'center';
            app.st_ed_anti.FontSize = 10.6666666666667;
            app.st_ed_anti.Enable = 'off';
            app.st_ed_anti.Position = [37 360 30.6666666666667 30];
            app.st_ed_anti.Value = '1';

            % Create N_ed_anti
            app.N_ed_anti = uieditfield(app.uipanel16, 'text');
            app.N_ed_anti.Tag = 'N_ed_anti';
            app.N_ed_anti.HorizontalAlignment = 'center';
            app.N_ed_anti.FontSize = 10.6666666666667;
            app.N_ed_anti.Enable = 'off';
            app.N_ed_anti.Position = [159 356 30.6666666666667 30];
            app.N_ed_anti.Value = '4';

            % Create text30
            app.text30 = uilabel(app.uipanel16);
            app.text30.Tag = 'text30';
            app.text30.BackgroundColor = [0.466666666666667 0.674509803921569 0.188235294117647];
            app.text30.HorizontalAlignment = 'center';
            app.text30.VerticalAlignment = 'top';
            app.text30.WordWrap = 'on';
            app.text30.FontSize = 10.6666666666667;
            app.text30.Position = [20 248 72 23.3333333333334];
            app.text30.Text = 'Start day of periodogram';

            % Create start_per
            app.start_per = uieditfield(app.uipanel16, 'text');
            app.start_per.Tag = 'start_per';
            app.start_per.HorizontalAlignment = 'center';
            app.start_per.FontSize = 10.6666666666667;
            app.start_per.Enable = 'off';
            app.start_per.Position = [33 207 30.6666666666667 30];
            app.start_per.Value = '5';

            % Create Triage_date_ed
            app.Triage_date_ed = uibuttongroup(app.uipanel16);
            app.Triage_date_ed.SelectionChangedFcn = createCallbackFcn(app, @Triage_date_ed_SelectionChangedFcn, true);
            app.Triage_date_ed.Title = 'Specific date for eduction/anticipation triage';
            app.Triage_date_ed.BackgroundColor = [0.466666666666667 0.674509803921569 0.188235294117647];
            app.Triage_date_ed.Tag = 'Triage_date_ed';
            app.Triage_date_ed.FontSize = 10.6666666666667;
            app.Triage_date_ed.Position = [254 458 227 47];

            % Create T_date_ed_Yes
            app.T_date_ed_Yes = uiradiobutton(app.Triage_date_ed);
            app.T_date_ed_Yes.Tag = 'T_date_ed_Yes';
            app.T_date_ed_Yes.Enable = 'off';
            app.T_date_ed_Yes.Text = 'Yes';
            app.T_date_ed_Yes.FontSize = 10.6666666666667;
            app.T_date_ed_Yes.Position = [25 5 54.6666666666667 24.6666666666667];

            % Create T_date_ed_No
            app.T_date_ed_No = uiradiobutton(app.Triage_date_ed);
            app.T_date_ed_No.Tag = 'T_date_ed_No';
            app.T_date_ed_No.Enable = 'off';
            app.T_date_ed_No.Text = 'No';
            app.T_date_ed_No.FontSize = 10.6666666666667;
            app.T_date_ed_No.Position = [118 4 53.3333333333333 24.6666666666667];
            app.T_date_ed_No.Value = true;

            % Create periodogram
            app.periodogram = uibuttongroup(app.uipanel16);
            app.periodogram.SelectionChangedFcn = createCallbackFcn(app, @periodogram_SelectionChangedFcn, true);
            app.periodogram.Title = 'Periodogram analysis';
            app.periodogram.BackgroundColor = [0.466666666666667 0.674509803921569 0.188235294117647];
            app.periodogram.Tag = 'periodogram';
            app.periodogram.FontSize = 10.6666666666667;
            app.periodogram.Position = [6 281 153 47];

            % Create periodo_Yes
            app.periodo_Yes = uiradiobutton(app.periodogram);
            app.periodo_Yes.Tag = 'periodo_Yes';
            app.periodo_Yes.Enable = 'off';
            app.periodo_Yes.Text = 'Yes';
            app.periodo_Yes.FontSize = 10.6666666666667;
            app.periodo_Yes.Position = [17 5 54.6666666666667 24];

            % Create periodo_No
            app.periodo_No = uiradiobutton(app.periodogram);
            app.periodo_No.Tag = 'periodo_No';
            app.periodo_No.Enable = 'off';
            app.periodo_No.Text = 'No';
            app.periodo_No.FontSize = 10.6666666666667;
            app.periodo_No.Position = [90 6 53.3333333333333 24.6666666666667];
            app.periodo_No.Value = true;

            % Create sig_level
            app.sig_level = uilabel(app.uipanel16);
            app.sig_level.Tag = 'sig_level';
            app.sig_level.BackgroundColor = [0.466666666666667 0.674509803921569 0.188235294117647];
            app.sig_level.HorizontalAlignment = 'center';
            app.sig_level.VerticalAlignment = 'top';
            app.sig_level.WordWrap = 'on';
            app.sig_level.FontSize = 10.6666666666667;
            app.sig_level.Position = [183 248 74 24];
            app.sig_level.Text = 'Significance level';

            % Create sig_val
            app.sig_val = uieditfield(app.uipanel16, 'text');
            app.sig_val.Tag = 'sig_val';
            app.sig_val.HorizontalAlignment = 'center';
            app.sig_val.FontSize = 10.6666666666667;
            app.sig_val.Enable = 'off';
            app.sig_val.Position = [197 205 44 32];
            app.sig_val.Value = '0.01';

            % Create min
            app.min = uilabel(app.uipanel16);
            app.min.Tag = 'min';
            app.min.BackgroundColor = [0.466666666666667 0.674509803921569 0.188235294117647];
            app.min.HorizontalAlignment = 'center';
            app.min.VerticalAlignment = 'top';
            app.min.WordWrap = 'on';
            app.min.FontSize = 10.6666666666667;
            app.min.Position = [42 168 50.6666666666667 24];
            app.min.Text = 'Minimum period';

            % Create min_period
            app.min_period = uieditfield(app.uipanel16, 'text');
            app.min_period.Tag = 'min_period';
            app.min_period.HorizontalAlignment = 'center';
            app.min_period.FontSize = 10.6666666666667;
            app.min_period.Enable = 'off';
            app.min_period.Position = [53 128 30.6666666666667 30];
            app.min_period.Value = '14';

            % Create text34
            app.text34 = uilabel(app.uipanel16);
            app.text34.Tag = 'text34';
            app.text34.BackgroundColor = [0.466666666666667 0.674509803921569 0.188235294117647];
            app.text34.HorizontalAlignment = 'center';
            app.text34.VerticalAlignment = 'top';
            app.text34.WordWrap = 'on';
            app.text34.FontSize = 10.6666666666667;
            app.text34.Position = [153 168 50.6666666666667 24];
            app.text34.Text = 'Maximum period';

            % Create max_period
            app.max_period = uieditfield(app.uipanel16, 'text');
            app.max_period.Tag = 'max_period';
            app.max_period.HorizontalAlignment = 'center';
            app.max_period.FontSize = 10.6666666666667;
            app.max_period.Enable = 'off';
            app.max_period.Position = [162 126 30.6666666666667 30];
            app.max_period.Value = '34';

            % Create Triage_date_per
            app.Triage_date_per = uibuttongroup(app.uipanel16);
            app.Triage_date_per.SelectionChangedFcn = createCallbackFcn(app, @Triage_date_per_SelectionChangedFcn, true);
            app.Triage_date_per.Title = 'Specific date for periodogram triage';
            app.Triage_date_per.BackgroundColor = [0.466666666666667 0.674509803921569 0.188235294117647];
            app.Triage_date_per.Tag = 'Triage_date_per';
            app.Triage_date_per.FontSize = 10.6666666666667;
            app.Triage_date_per.Position = [279 213 192 39];

            % Create T_date_Yes_per
            app.T_date_Yes_per = uiradiobutton(app.Triage_date_per);
            app.T_date_Yes_per.Tag = 'T_date_Yes_per';
            app.T_date_Yes_per.Enable = 'off';
            app.T_date_Yes_per.Text = 'Yes';
            app.T_date_Yes_per.FontSize = 10.6666666666667;
            app.T_date_Yes_per.Position = [24 6 54.6666666666667 19.3333333333333];

            % Create T_date_No_per
            app.T_date_No_per = uiradiobutton(app.Triage_date_per);
            app.T_date_No_per.Tag = 'T_date_No_per';
            app.T_date_No_per.Enable = 'off';
            app.T_date_No_per.Text = 'No';
            app.T_date_No_per.FontSize = 10.6666666666667;
            app.T_date_No_per.Position = [128 5 53.3333333333333 19.3333333333333];
            app.T_date_No_per.Value = true;

            % Create text31
            app.text31 = uilabel(app.uipanel16);
            app.text31.Tag = 'text31';
            app.text31.BackgroundColor = [0.466666666666667 0.674509803921569 0.188235294117647];
            app.text31.HorizontalAlignment = 'center';
            app.text31.VerticalAlignment = 'top';
            app.text31.WordWrap = 'on';
            app.text31.FontSize = 10.6666666666667;
            app.text31.Position = [95 245 88 28];
            app.text31.Text = 'Number of days of periodogram';

            % Create N_per
            app.N_per = uieditfield(app.uipanel16, 'text');
            app.N_per.Tag = 'N_per';
            app.N_per.HorizontalAlignment = 'center';
            app.N_per.FontSize = 10.6666666666667;
            app.N_per.Enable = 'off';
            app.N_per.Position = [122 207 30.6666666666667 30];
            app.N_per.Value = '7';

            % Create Acto
            app.Acto = uibuttongroup(app.uipanel16);
            app.Acto.SelectionChangedFcn = createCallbackFcn(app, @Acto_SelectionChangedFcn, true);
            app.Acto.Title = 'Actogram';
            app.Acto.BackgroundColor = [0.466666666666667 0.674509803921569 0.188235294117647];
            app.Acto.Tag = 'Acto';
            app.Acto.FontSize = 10.6666666666667;
            app.Acto.Position = [73 75 129 47];

            % Create act_Yes
            app.act_Yes = uiradiobutton(app.Acto);
            app.act_Yes.Tag = 'act_Yes';
            app.act_Yes.Enable = 'off';
            app.act_Yes.Text = 'Yes';
            app.act_Yes.FontSize = 10.6666666666667;
            app.act_Yes.Position = [17 5 54.6666666666667 24];

            % Create act_No
            app.act_No = uiradiobutton(app.Acto);
            app.act_No.Tag = 'act_No';
            app.act_No.Enable = 'off';
            app.act_No.Text = 'No';
            app.act_No.FontSize = 10.6666666666667;
            app.act_No.Position = [68 6 53.3333333333333 24.6666666666667];
            app.act_No.Value = true;

            % Create text37
            app.text37 = uilabel(app.uipanel16);
            app.text37.Tag = 'text37';
            app.text37.BackgroundColor = [0.466666666666667 0.674509803921569 0.188235294117647];
            app.text37.HorizontalAlignment = 'center';
            app.text37.VerticalAlignment = 'top';
            app.text37.WordWrap = 'on';
            app.text37.FontSize = 10.6666666666667;
            app.text37.Position = [31 44 122 23.3333333333333];
            app.text37.Text = 'Start day of actogram';

            % Create st_act
            app.st_act = uieditfield(app.uipanel16, 'text');
            app.st_act.Tag = 'st_act';
            app.st_act.HorizontalAlignment = 'center';
            app.st_act.FontSize = 10.6666666666667;
            app.st_act.Enable = 'off';
            app.st_act.Position = [72 14 30.6666666666667 30];
            app.st_act.Value = '1';

            % Create N_act
            app.N_act = uieditfield(app.uipanel16, 'text');
            app.N_act.Tag = 'N_act';
            app.N_act.HorizontalAlignment = 'center';
            app.N_act.FontSize = 10.6666666666667;
            app.N_act.Enable = 'off';
            app.N_act.Position = [221 9 30.6666666666667 30];
            app.N_act.Value = '11';

            % Create text38
            app.text38 = uilabel(app.uipanel16);
            app.text38.Tag = 'text38';
            app.text38.BackgroundColor = [0.466666666666667 0.674509803921569 0.188235294117647];
            app.text38.HorizontalAlignment = 'center';
            app.text38.VerticalAlignment = 'top';
            app.text38.WordWrap = 'on';
            app.text38.FontSize = 10.6666666666667;
            app.text38.Position = [195 38 88 28];
            app.text38.Text = 'Number of days of actogram';

            % Create Period_figure
            app.Period_figure = uibuttongroup(app.uipanel16);
            app.Period_figure.SelectionChangedFcn = createCallbackFcn(app, @Period_figure_SelectionChangedFcn, true);
            app.Period_figure.Title = ' Periodogram figure';
            app.Period_figure.BackgroundColor = [0.466666666666667 0.674509803921569 0.188235294117647];
            app.Period_figure.Tag = 'Period_figure';
            app.Period_figure.FontSize = 10.6666666666667;
            app.Period_figure.Position = [184 281 158 47];

            % Create PF_Yes
            app.PF_Yes = uiradiobutton(app.Period_figure);
            app.PF_Yes.Tag = 'PF_Yes';
            app.PF_Yes.Enable = 'off';
            app.PF_Yes.Text = 'Yes';
            app.PF_Yes.FontSize = 10.6666666666667;
            app.PF_Yes.Position = [17 5 54.6666666666667 24];

            % Create PF_No
            app.PF_No = uiradiobutton(app.Period_figure);
            app.PF_No.Tag = 'PF_No';
            app.PF_No.Enable = 'off';
            app.PF_No.Text = 'No';
            app.PF_No.FontSize = 10.6666666666667;
            app.PF_No.Position = [94 5 53.3333333333333 24.6666666666667];
            app.PF_No.Value = true;

            % Create Triaged_cndn_ed
            app.Triaged_cndn_ed = uibuttongroup(app.uipanel16);
            app.Triaged_cndn_ed.Title = 'Triage condition for sleep eduction analysis';
            app.Triaged_cndn_ed.BackgroundColor = [0.466666666666667 0.674509803921569 0.188235294117647];
            app.Triaged_cndn_ed.Tag = 'Triaged_cndn_ed';
            app.Triaged_cndn_ed.FontSize = 10.6666666666667;
            app.Triaged_cndn_ed.Position = [239 335 241 119];

            % Create ed_ch1
            app.ed_ch1 = uicheckbox(app.Triaged_cndn_ed);
            app.ed_ch1.Tag = 'ed_ch1';
            app.ed_ch1.Enable = 'off';
            app.ed_ch1.Text = ' Last day of activity <';
            app.ed_ch1.FontSize = 10.6666666666667;
            app.ed_ch1.Position = [20 81 124.666666666667 19.3333333333333];
            app.ed_ch1.Value = true;

            % Create ed_ch2
            app.ed_ch2 = uicheckbox(app.Triaged_cndn_ed);
            app.ed_ch2.Tag = 'ed_ch2';
            app.ed_ch2.Enable = 'off';
            app.ed_ch2.Text = ' Last day activity(Lon-6 to Lon+6) <';
            app.ed_ch2.FontSize = 10.6666666666667;
            app.ed_ch2.Position = [9 44 192 19.3333333333333];

            % Create ed_ch3
            app.ed_ch3 = uicheckbox(app.Triaged_cndn_ed);
            app.ed_ch3.Tag = 'ed_ch3';
            app.ed_ch3.Enable = 'off';
            app.ed_ch3.Text = 'Activity per 30 min <';
            app.ed_ch3.FontSize = 10.6666666666667;
            app.ed_ch3.Position = [11 13 131.333333333333 20.6666666666667];

            % Create Threshold3_ed
            app.Threshold3_ed = uieditfield(app.Triaged_cndn_ed, 'text');
            app.Threshold3_ed.Tag = 'Threshold3_ed';
            app.Threshold3_ed.HorizontalAlignment = 'center';
            app.Threshold3_ed.FontSize = 10.6666666666667;
            app.Threshold3_ed.Enable = 'off';
            app.Threshold3_ed.Position = [140 74 28.6666666666667 30];
            app.Threshold3_ed.Value = '5';

            % Create Threshold2_ed
            app.Threshold2_ed = uieditfield(app.Triaged_cndn_ed, 'text');
            app.Threshold2_ed.Tag = 'Threshold2_ed';
            app.Threshold2_ed.HorizontalAlignment = 'center';
            app.Threshold2_ed.FontSize = 10.6666666666667;
            app.Threshold2_ed.Enable = 'off';
            app.Threshold2_ed.Position = [195 40 28.6666666666667 30];
            app.Threshold2_ed.Value = '2';

            % Create Threshold4_ed
            app.Threshold4_ed = uieditfield(app.Triaged_cndn_ed, 'text');
            app.Threshold4_ed.Tag = 'Threshold4_ed';
            app.Threshold4_ed.HorizontalAlignment = 'center';
            app.Threshold4_ed.FontSize = 10.6666666666667;
            app.Threshold4_ed.Enable = 'off';
            app.Threshold4_ed.Position = [130 10 26.6666666666667 30];
            app.Threshold4_ed.Value = '1';

            % Create text26
            app.text26 = uilabel(app.uipanel16);
            app.text26.Tag = 'text26';
            app.text26.BackgroundColor = [0.466666666666667 0.674509803921569 0.188235294117647];
            app.text26.HorizontalAlignment = 'center';
            app.text26.VerticalAlignment = 'top';
            app.text26.WordWrap = 'on';
            app.text26.FontSize = 10.6666666666667;
            app.text26.Position = [7 397 112 26];
            app.text26.Text = 'Start day of eduction/anticipation';

            % Create text27
            app.text27 = uilabel(app.uipanel16);
            app.text27.Tag = 'text27';
            app.text27.BackgroundColor = [0.466666666666667 0.674509803921569 0.188235294117647];
            app.text27.HorizontalAlignment = 'center';
            app.text27.VerticalAlignment = 'top';
            app.text27.WordWrap = 'on';
            app.text27.FontSize = 10.6666666666667;
            app.text27.Position = [125 397 107 28];
            app.text27.Text = 'Number of days of eduction / anticipation';

            % Create uibuttongroup27
            app.uibuttongroup27 = uibuttongroup(app.uipanel16);
            app.uibuttongroup27.Title = 'Triage condition for periodogram';
            app.uibuttongroup27.BackgroundColor = [0.466666666666667 0.674509803921569 0.188235294117647];
            app.uibuttongroup27.Tag = 'uibuttongroup27';
            app.uibuttongroup27.FontSize = 10.6666666666667;
            app.uibuttongroup27.Position = [254 87 233 115];

            % Create per_ch1
            app.per_ch1 = uicheckbox(app.uibuttongroup27);
            app.per_ch1.Tag = 'per_ch1';
            app.per_ch1.Enable = 'off';
            app.per_ch1.Text = ' Last day of activity <';
            app.per_ch1.FontSize = 10.6666666666667;
            app.per_ch1.Position = [7 74 124.666666666667 19.3333333333333];
            app.per_ch1.Value = true;

            % Create per_ch2
            app.per_ch2 = uicheckbox(app.uibuttongroup27);
            app.per_ch2.Tag = 'per_ch2';
            app.per_ch2.Enable = 'off';
            app.per_ch2.Text = ' Last day activity(Lon-6 to Lon+6) <';
            app.per_ch2.FontSize = 10.6666666666667;
            app.per_ch2.Position = [6 39 192 19.3333333333333];

            % Create per_ch3
            app.per_ch3 = uicheckbox(app.uibuttongroup27);
            app.per_ch3.Tag = 'per_ch3';
            app.per_ch3.Enable = 'off';
            app.per_ch3.Text = 'Activity per 30 min <';
            app.per_ch3.FontSize = 10.6666666666667;
            app.per_ch3.Position = [6 10 132 19.3333333333333];

            % Create Threshold3_per
            app.Threshold3_per = uieditfield(app.uibuttongroup27, 'text');
            app.Threshold3_per.Tag = 'Threshold3_per';
            app.Threshold3_per.HorizontalAlignment = 'center';
            app.Threshold3_per.FontSize = 10.6666666666667;
            app.Threshold3_per.Enable = 'off';
            app.Threshold3_per.Position = [128 69 30.6666666666667 30];
            app.Threshold3_per.Value = '5';

            % Create Threshold2_per
            app.Threshold2_per = uieditfield(app.uibuttongroup27, 'text');
            app.Threshold2_per.Tag = 'Threshold2_per';
            app.Threshold2_per.HorizontalAlignment = 'center';
            app.Threshold2_per.FontSize = 10.6666666666667;
            app.Threshold2_per.Enable = 'off';
            app.Threshold2_per.Position = [186 32 24.6666666666667 30];
            app.Threshold2_per.Value = '2';

            % Create Threshold4_per
            app.Threshold4_per = uieditfield(app.uibuttongroup27, 'text');
            app.Threshold4_per.Tag = 'Threshold4_per';
            app.Threshold4_per.HorizontalAlignment = 'center';
            app.Threshold4_per.FontSize = 10.6666666666667;
            app.Threshold4_per.Enable = 'off';
            app.Threshold4_per.Position = [122 7 34 30];
            app.Threshold4_per.Value = '1';

            % Create uipanel15
            app.uipanel15 = uipanel(app.figure1);
            app.uipanel15.Title = 'Sleep analysis';
            app.uipanel15.BackgroundColor = [0.929411764705882 0.694117647058824 0.125490196078431];
            app.uipanel15.Tag = 'uipanel15';
            app.uipanel15.FontSize = 21.3333333333333;
            app.uipanel15.Position = [234 1 313 586];

            % Create sleep_ana
            app.sleep_ana = uibuttongroup(app.uipanel15);
            app.sleep_ana.SelectionChangedFcn = createCallbackFcn(app, @sleep_ana_SelectionChangedFcn, true);
            app.sleep_ana.Title = ' Sleep analysis';
            app.sleep_ana.BackgroundColor = [0.929411764705882 0.694117647058824 0.125490196078431];
            app.sleep_ana.Tag = 'sleep_ana';
            app.sleep_ana.FontSize = 13.3333333333333;
            app.sleep_ana.Position = [15 503 127 47];

            % Create SL_Yes
            app.SL_Yes = uiradiobutton(app.sleep_ana);
            app.SL_Yes.Tag = 'SL_Yes';
            app.SL_Yes.Enable = 'off';
            app.SL_Yes.Text = 'Yes';
            app.SL_Yes.FontSize = 10.6666666666667;
            app.SL_Yes.Position = [13 6 54.4 24.7692307692308];
            app.SL_Yes.Value = true;

            % Create SL_No
            app.SL_No = uiradiobutton(app.sleep_ana);
            app.SL_No.Tag = 'SL_No';
            app.SL_No.Enable = 'off';
            app.SL_No.Text = 'No';
            app.SL_No.FontSize = 10.6666666666667;
            app.SL_No.Position = [64 4 54.6666666666667 24.6666666666667];

            % Create text32
            app.text32 = uilabel(app.uipanel15);
            app.text32.Tag = 'text32';
            app.text32.BackgroundColor = [0.929411764705882 0.694117647058824 0.125490196078431];
            app.text32.HorizontalAlignment = 'center';
            app.text32.VerticalAlignment = 'top';
            app.text32.WordWrap = 'on';
            app.text32.FontSize = 10.6666666666667;
            app.text32.Position = [8 465 90.6666666666667 36.6666666666667];
            app.text32.Text = 'Start day of sleep analysis';

            % Create st_day_sleep
            app.st_day_sleep = uieditfield(app.uipanel15, 'text');
            app.st_day_sleep.Tag = 'st_day_sleep';
            app.st_day_sleep.HorizontalAlignment = 'center';
            app.st_day_sleep.FontSize = 10.6666666666667;
            app.st_day_sleep.Enable = 'off';
            app.st_day_sleep.Position = [36 440 30.6666666666667 30.0000000000001];
            app.st_day_sleep.Value = '1';

            % Create text2
            app.text2 = uilabel(app.uipanel15);
            app.text2.Tag = 'text2';
            app.text2.BackgroundColor = [0.929411764705882 0.694117647058824 0.125490196078431];
            app.text2.HorizontalAlignment = 'center';
            app.text2.VerticalAlignment = 'top';
            app.text2.WordWrap = 'on';
            app.text2.FontSize = 10.6666666666667;
            app.text2.Position = [126 466 90.6666666666667 36.6666666666667];
            app.text2.Text = 'Number of days of sleep analysis';

            % Create N_days_ana
            app.N_days_ana = uieditfield(app.uipanel15, 'text');
            app.N_days_ana.Tag = 'N_days_ana';
            app.N_days_ana.HorizontalAlignment = 'center';
            app.N_days_ana.FontSize = 10.6666666666667;
            app.N_days_ana.Enable = 'off';
            app.N_days_ana.Position = [145 442 30.6666666666667 30.0000000000001];
            app.N_days_ana.Value = '4';

            % Create analysis_per_interval
            app.analysis_per_interval = uibuttongroup(app.uipanel15);
            app.analysis_per_interval.SelectionChangedFcn = createCallbackFcn(app, @analysis_per_interval_SelectionChangeFcn, true);
            app.analysis_per_interval.Title = 'Sleep Analysis/Interval';
            app.analysis_per_interval.BackgroundColor = [0.929411764705882 0.694117647058824 0.125490196078431];
            app.analysis_per_interval.Tag = 'analysis_per_interval';
            app.analysis_per_interval.FontSize = 10.6666666666667;
            app.analysis_per_interval.Position = [159 505 131 47];

            % Create Int_Yes
            app.Int_Yes = uiradiobutton(app.analysis_per_interval);
            app.Int_Yes.Tag = 'Int_Yes';
            app.Int_Yes.Enable = 'off';
            app.Int_Yes.Text = 'Yes';
            app.Int_Yes.FontSize = 10.6666666666667;
            app.Int_Yes.Position = [10 5 54.6666666666667 24.6666666666667];

            % Create Int_No
            app.Int_No = uiradiobutton(app.analysis_per_interval);
            app.Int_No.Tag = 'Int_No';
            app.Int_No.Enable = 'off';
            app.Int_No.Text = 'No';
            app.Int_No.FontSize = 10.6666666666667;
            app.Int_No.Position = [57 5 44.6666666666667 24.6666666666667];
            app.Int_No.Value = true;

            % Create text21
            app.text21 = uilabel(app.uipanel15);
            app.text21.Tag = 'text21';
            app.text21.BackgroundColor = [0.929411764705882 0.694117647058824 0.125490196078431];
            app.text21.HorizontalAlignment = 'center';
            app.text21.VerticalAlignment = 'top';
            app.text21.WordWrap = 'on';
            app.text21.FontSize = 10.6666666666667;
            app.text21.Position = [234 482 73.3333333333334 14];
            app.text21.Text = 'Per interval (hr)';

            % Create per_Interval
            app.per_Interval = uieditfield(app.uipanel15, 'text');
            app.per_Interval.Tag = 'per_Interval';
            app.per_Interval.HorizontalAlignment = 'center';
            app.per_Interval.FontSize = 10.6666666666667;
            app.per_Interval.Enable = 'off';
            app.per_Interval.Position = [250 442 30.6666666666667 30];
            app.per_Interval.Value = '1';

            % Create Sleep_dep
            app.Sleep_dep = uibuttongroup(app.uipanel15);
            app.Sleep_dep.SelectionChangedFcn = createCallbackFcn(app, @Sleep_dep_SelectionChangeFcn, true);
            app.Sleep_dep.Title = 'Sleep Deprivation';
            app.Sleep_dep.BackgroundColor = [0.929411764705882 0.694117647058824 0.125490196078431];
            app.Sleep_dep.Tag = 'Sleep_dep';
            app.Sleep_dep.FontSize = 10.6666666666667;
            app.Sleep_dep.Position = [16 386 107 47];

            % Create SD_Yes
            app.SD_Yes = uiradiobutton(app.Sleep_dep);
            app.SD_Yes.Tag = 'SD_Yes';
            app.SD_Yes.Enable = 'off';
            app.SD_Yes.Text = 'Yes';
            app.SD_Yes.FontSize = 10.6666666666667;
            app.SD_Yes.Position = [6 6 54.6666666666667 24.6666666666667];

            % Create SD_No
            app.SD_No = uiradiobutton(app.Sleep_dep);
            app.SD_No.Tag = 'SD_No';
            app.SD_No.Enable = 'off';
            app.SD_No.Text = 'No';
            app.SD_No.FontSize = 10.6666666666667;
            app.SD_No.Position = [56 6 39.3333333333333 24.6666666666667];
            app.SD_No.Value = true;

            % Create text9
            app.text9 = uilabel(app.uipanel15);
            app.text9.Tag = 'text9';
            app.text9.BackgroundColor = [0.929411764705882 0.694117647058824 0.125490196078431];
            app.text9.HorizontalAlignment = 'center';
            app.text9.VerticalAlignment = 'top';
            app.text9.WordWrap = 'on';
            app.text9.FontSize = 10.6666666666667;
            app.text9.Position = [6 286 60.6666666666667 15.3333333333333];
            app.text9.Text = 'SD day';

            % Create SD_day
            app.SD_day = uieditfield(app.uipanel15, 'text');
            app.SD_day.Tag = 'SD_day';
            app.SD_day.HorizontalAlignment = 'center';
            app.SD_day.FontSize = 10.6666666666667;
            app.SD_day.Enable = 'off';
            app.SD_day.Position = [20 250 30.6666666666667 30];
            app.SD_day.Value = '3';

            % Create text10
            app.text10 = uilabel(app.uipanel15);
            app.text10.Tag = 'text10';
            app.text10.BackgroundColor = [0.929411764705882 0.694117647058824 0.125490196078431];
            app.text10.HorizontalAlignment = 'center';
            app.text10.VerticalAlignment = 'top';
            app.text10.WordWrap = 'on';
            app.text10.FontSize = 10.6666666666667;
            app.text10.Position = [82 275 108 42.6666666666667];
            app.text10.Text = 'SD starting time (Clock time, 24 hr format)';

            % Create SD_hr
            app.SD_hr = uieditfield(app.uipanel15, 'text');
            app.SD_hr.Tag = 'SD_hr';
            app.SD_hr.HorizontalAlignment = 'center';
            app.SD_hr.FontSize = 10.6666666666667;
            app.SD_hr.Enable = 'off';
            app.SD_hr.Position = [116 248 30.6666666666667 30];
            app.SD_hr.Value = '20';

            % Create text11
            app.text11 = uilabel(app.uipanel15);
            app.text11.Tag = 'text11';
            app.text11.BackgroundColor = [0.929411764705882 0.694117647058824 0.125490196078431];
            app.text11.HorizontalAlignment = 'center';
            app.text11.VerticalAlignment = 'top';
            app.text11.WordWrap = 'on';
            app.text11.FontSize = 10.6666666666667;
            app.text11.Position = [190 286 76 16];
            app.text11.Text = 'SD duration';

            % Create SD_dur
            app.SD_dur = uieditfield(app.uipanel15, 'text');
            app.SD_dur.Tag = 'SD_dur';
            app.SD_dur.HorizontalAlignment = 'center';
            app.SD_dur.FontSize = 10.6666666666667;
            app.SD_dur.Enable = 'off';
            app.SD_dur.Position = [216 248 30.6666666666667 30];
            app.SD_dur.Value = '12';

            % Create SD_Th
            app.SD_Th = uibuttongroup(app.uipanel15);
            app.SD_Th.SelectionChangedFcn = createCallbackFcn(app, @SD_Th_SelectionChangeFcn, true);
            app.SD_Th.Title = ' Minimum SD% threshold';
            app.SD_Th.BackgroundColor = [0.929411764705882 0.694117647058824 0.125490196078431];
            app.SD_Th.Tag = 'SD_Th';
            app.SD_Th.FontSize = 10.6666666666667;
            app.SD_Th.Position = [16 188 162 47];

            % Create SD_Th_Yes
            app.SD_Th_Yes = uiradiobutton(app.SD_Th);
            app.SD_Th_Yes.Tag = 'SD_Th_Yes';
            app.SD_Th_Yes.Enable = 'off';
            app.SD_Th_Yes.Text = 'Yes';
            app.SD_Th_Yes.FontSize = 10.6666666666667;
            app.SD_Th_Yes.Position = [13 6 54.4 24.7692307692308];
            app.SD_Th_Yes.Value = true;

            % Create SD_Th_No
            app.SD_Th_No = uiradiobutton(app.SD_Th);
            app.SD_Th_No.Tag = 'SD_Th_No';
            app.SD_Th_No.Enable = 'off';
            app.SD_Th_No.Text = 'No';
            app.SD_Th_No.FontSize = 10.6666666666667;
            app.SD_Th_No.Position = [88 5 54.6666666666667 24.6666666666667];

            % Create text12
            app.text12 = uilabel(app.uipanel15);
            app.text12.Tag = 'text12';
            app.text12.BackgroundColor = [0.929411764705882 0.694117647058824 0.125490196078431];
            app.text12.HorizontalAlignment = 'center';
            app.text12.VerticalAlignment = 'top';
            app.text12.WordWrap = 'on';
            app.text12.FontSize = 10.6666666666667;
            app.text12.Position = [194 221 76 15.3333333333333];
            app.text12.Text = 'Minimum SD %';

            % Create SD_per
            app.SD_per = uieditfield(app.uipanel15, 'text');
            app.SD_per.Tag = 'SD_per';
            app.SD_per.HorizontalAlignment = 'center';
            app.SD_per.FontSize = 10.6666666666667;
            app.SD_per.BackgroundColor = [0.929411764705882 0.694117647058824 0.125490196078431];
            app.SD_per.Enable = 'off';
            app.SD_per.Position = [214 189 30.6666666666667 30];
            app.SD_per.Value = '90';

            % Create Triage_date
            app.Triage_date = uibuttongroup(app.uipanel15);
            app.Triage_date.SelectionChangedFcn = createCallbackFcn(app, @Triage_date_SelectionChangedFcn, true);
            app.Triage_date.Title = ' Specific date for triage';
            app.Triage_date.BackgroundColor = [0.929411764705882 0.694117647058824 0.125490196078431];
            app.Triage_date.Tag = 'Triage_date';
            app.Triage_date.FontSize = 10.6666666666667;
            app.Triage_date.Position = [20 138 164 47];

            % Create T_date_Yes
            app.T_date_Yes = uiradiobutton(app.Triage_date);
            app.T_date_Yes.Tag = 'T_date_Yes';
            app.T_date_Yes.Enable = 'off';
            app.T_date_Yes.Text = 'Yes';
            app.T_date_Yes.FontSize = 10.6666666666667;
            app.T_date_Yes.Position = [20 6 54.6666666666667 24.6666666666667];

            % Create T_date_No
            app.T_date_No = uiradiobutton(app.Triage_date);
            app.T_date_No.Tag = 'T_date_No';
            app.T_date_No.Enable = 'off';
            app.T_date_No.Text = 'No';
            app.T_date_No.FontSize = 10.6666666666667;
            app.T_date_No.Position = [92 5 53.3333333333333 24.6666666666667];
            app.T_date_No.Value = true;

            % Create Triaged_cndn
            app.Triaged_cndn = uibuttongroup(app.uipanel15);
            app.Triaged_cndn.Title = 'Triage condition for sleep analysis';
            app.Triaged_cndn.BackgroundColor = [0.929411764705882 0.694117647058824 0.125490196078431];
            app.Triaged_cndn.Tag = 'Triaged_cndn';
            app.Triaged_cndn.FontSize = 10.6666666666667;
            app.Triaged_cndn.Position = [4 11 244 125];

            % Create SL_ch1
            app.SL_ch1 = uicheckbox(app.Triaged_cndn);
            app.SL_ch1.Tag = 'SL_ch1';
            app.SL_ch1.Enable = 'off';
            app.SL_ch1.Text = ' Last day of activity <';
            app.SL_ch1.FontSize = 10.6666666666667;
            app.SL_ch1.Position = [6 86 124.666666666667 19.3333333333333];
            app.SL_ch1.Value = true;

            % Create SL_ch2
            app.SL_ch2 = uicheckbox(app.Triaged_cndn);
            app.SL_ch2.Tag = 'SL_ch2';
            app.SL_ch2.Enable = 'off';
            app.SL_ch2.Text = ' Last day activity(Lon-6 to Lon+6) <';
            app.SL_ch2.FontSize = 10.6666666666667;
            app.SL_ch2.Position = [6 50 192 19.3333333333333];

            % Create SL_ch3
            app.SL_ch3 = uicheckbox(app.Triaged_cndn);
            app.SL_ch3.Tag = 'SL_ch3';
            app.SL_ch3.Enable = 'off';
            app.SL_ch3.Text = 'Activity per 30 min <';
            app.SL_ch3.FontSize = 10.6666666666667;
            app.SL_ch3.Position = [5 13 132 19.3333333333333];

            % Create Threshold4
            app.Threshold4 = uieditfield(app.Triaged_cndn, 'text');
            app.Threshold4.Tag = 'Threshold4';
            app.Threshold4.HorizontalAlignment = 'center';
            app.Threshold4.FontSize = 10.6666666666667;
            app.Threshold4.Enable = 'off';
            app.Threshold4.Position = [128 8 29.3333333333333 30];
            app.Threshold4.Value = '1';

            % Create Threshold2
            app.Threshold2 = uieditfield(app.Triaged_cndn, 'text');
            app.Threshold2.Tag = 'Threshold2';
            app.Threshold2.HorizontalAlignment = 'center';
            app.Threshold2.FontSize = 10.6666666666667;
            app.Threshold2.Enable = 'off';
            app.Threshold2.Position = [193 44 28.6666666666667 30];
            app.Threshold2.Value = '2';

            % Create Threshold3
            app.Threshold3 = uieditfield(app.Triaged_cndn, 'text');
            app.Threshold3.Tag = 'Threshold3';
            app.Threshold3.HorizontalAlignment = 'center';
            app.Threshold3.FontSize = 10.6666666666667;
            app.Threshold3.Enable = 'off';
            app.Threshold3.Position = [126 80 28.6666666666667 30];
            app.Threshold3.Value = '5';

            % Create text47
            app.text47 = uilabel(app.uipanel15);
            app.text47.Tag = 'text47';
            app.text47.BackgroundColor = [0.929411764705882 0.694117647058824 0.125490196078431];
            app.text47.HorizontalAlignment = 'center';
            app.text47.VerticalAlignment = 'top';
            app.text47.WordWrap = 'on';
            app.text47.FontSize = 10.6666666666667;
            app.text47.Position = [40 362 68 16];
            app.text47.Text = 'SL interval';

            % Create SD_int
            app.SD_int = uieditfield(app.uipanel15, 'text');
            app.SD_int.Tag = 'SD_int';
            app.SD_int.HorizontalAlignment = 'center';
            app.SD_int.FontSize = 10.6666666666667;
            app.SD_int.Enable = 'off';
            app.SD_int.Position = [56 330 28 30];
            app.SD_int.Value = '12';

            % Create SG_int
            app.SG_int = uieditfield(app.uipanel15, 'text');
            app.SG_int.Tag = 'SG_int';
            app.SG_int.HorizontalAlignment = 'center';
            app.SG_int.FontSize = 10.6666666666667;
            app.SG_int.Enable = 'off';
            app.SG_int.Position = [199 322 28 30];
            app.SG_int.Value = '12';

            % Create text48
            app.text48 = uilabel(app.uipanel15);
            app.text48.Tag = 'text48';
            app.text48.BackgroundColor = [0.929411764705882 0.694117647058824 0.125490196078431];
            app.text48.HorizontalAlignment = 'center';
            app.text48.VerticalAlignment = 'top';
            app.text48.WordWrap = 'on';
            app.text48.FontSize = 10.6666666666667;
            app.text48.Position = [183 356 68 23.3333333333334];
            app.text48.Text = 'Rebound interval';

            % Create cum_fig
            app.cum_fig = uibuttongroup(app.uipanel15);
            app.cum_fig.SelectionChangedFcn = createCallbackFcn(app, @cum_fig_SelectionChangedFcn, true);
            app.cum_fig.Title = 'Cumulative sleep figure';
            app.cum_fig.BackgroundColor = [0.929411764705882 0.694117647058824 0.125490196078431];
            app.cum_fig.Tag = 'cum_fig';
            app.cum_fig.FontSize = 10.6666666666667;
            app.cum_fig.Position = [145 386 143 47];

            % Create cum_Yes
            app.cum_Yes = uiradiobutton(app.cum_fig);
            app.cum_Yes.Tag = 'cum_Yes';
            app.cum_Yes.Enable = 'off';
            app.cum_Yes.Text = 'Yes';
            app.cum_Yes.FontSize = 10.6666666666667;
            app.cum_Yes.Position = [13 9 54.4 24.7692307692308];

            % Create cum_No
            app.cum_No = uiradiobutton(app.cum_fig);
            app.cum_No.Tag = 'cum_No';
            app.cum_No.Enable = 'off';
            app.cum_No.Text = 'No';
            app.cum_No.FontSize = 10.6666666666667;
            app.cum_No.Position = [66 6 54.6666666666667 24.6666666666667];
            app.cum_No.Value = true;

            % Create M_box
            app.M_box = uitextarea(app.figure1);
            app.M_box.Tag = 'M_box';
            app.M_box.HorizontalAlignment = 'center';
            app.M_box.FontSize = 13.3333333333333;
            app.M_box.BackgroundColor = [0.301960784313725 0.745098039215686 0.933333333333333];
            app.M_box.Enable = 'off';
            app.M_box.Position = [13 54.3333333333333 212 69.3333333333333];
            app.M_box.Value = {'Status update'};

            % Create uipanel17
            app.uipanel17 = uipanel(app.figure1);
            app.uipanel17.Title = 'Vesrion information';
            app.uipanel17.Tag = 'uipanel17';
            app.uipanel17.FontSize = 10.6666666666667;
            app.uipanel17.Position = [15 2 204 51];

            % Create text43
            app.text43 = uilabel(app.uipanel17);
            app.text43.Tag = 'text43';
            app.text43.HorizontalAlignment = 'center';
            app.text43.VerticalAlignment = 'top';
            app.text43.WordWrap = 'on';
            app.text43.FontSize = 10.6666666666667;
            app.text43.Position = [39 16 106.666666666667 13.3333333333333];
            app.text43.Text = 'V2024.2 (Gamma)';

            % Create text44
            app.text44 = uilabel(app.uipanel17);
            app.text44.Tag = 'text44';
            app.text44.HorizontalAlignment = 'center';
            app.text44.VerticalAlignment = 'top';
            app.text44.WordWrap = 'on';
            app.text44.FontSize = 10.6666666666667;
            app.text44.Position = [23 4 134 14];
            app.text44.Text = 'Date modified : 10/09/2024';

            % Show the figure after all components are created
            app.figure1.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = sleepmat(varargin)

            runningApp = getRunningApp(app);

            % Check for running singleton app
            if isempty(runningApp)

                % Create UIFigure and components
                createComponents(app)

                % Register the app with App Designer
                registerApp(app, app.figure1)

                % Execute the startup function
                runStartupFcn(app, @(app)sleepmat_OpeningFcn(app, varargin{:}))
            else

                % Focus the running singleton app
                figure(runningApp.figure1)

                app = runningApp;
            end

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.figure1)
        end
    end
end