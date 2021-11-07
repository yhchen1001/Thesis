%% Clear workspace and all figures
clear all;clc;close all

get(0,'Factory');
set(0,'defaultfigurecolor',[1 1 1]);

%% Load Data
[FileName,Path2RadarData,filter_index]=uigetfile('05  ThreeHumans_50m_400M_all_walking_SomeHandsInPockets','D:\EEE4022S\Dataset\OneDrive_2021-08-23\07  Radar datasets');
RadarData = load([Path2RadarData filesep FileName]);

%%====== Extract Range Profiles before, after Equalisation and after Notch filtering ===============================================================================================================

RangeProfiles_BeforeEq = RadarData.RangeLines_BeforeEq;
RangeProfiles_AfterEq = RadarData.RangeLines_AfterEq;
RangeProfiles_AfterEqNotch = RadarData.RangeLines_AfterEQ_Notch;

figure; imagesc(20*log10(abs(RangeProfiles_AfterEqNotch)));

% Remove unwanted data
RangeProfiles_AfterEqNotch = RangeProfiles_AfterEqNotch(:, 400: end);        % Dr Abdul Gaffar - for first dataset 01  SingleHuman_50m_400M__07_12_2011_walking


%% Extract other radar parameters

 PRF_Hz = RadarData.Info.PRF_Hz;
 Bandwidth_Hz = RadarData.Info.Bandwidth_Hz;
 RangeStart_m = RadarData.Info.RangeStart_m;
 BlindRange_m = RadarData.Info.BlindRange_m;
 
[NumOfPulses,NumOfRangeBins]=size(RangeProfiles_AfterEqNotch);


%% =============Plot Range Profiles=====================================================================================================================================

fontsize1 = 12;
clims = [-40 0];

% Normalise data to have a peak of 0dB or 1 in linear scale
[MaxRangeLine MaxIdx] = max(max(abs(RangeProfiles_AfterEqNotch)));

% Plot range lines
figure; axes('fontsize',fontsize1);
imagesc(20*log10(abs(RangeProfiles_AfterEqNotch)./MaxRangeLine),clims);
colorbar;
xlabel('Range (bins)','fontsize',fontsize1);
ylabel('Number of pulses','fontsize',fontsize1);
% %title('Range lines: after Eq Notch (1 person)','fontsize',fontsize1);


%%=======Input============================================================================================================================================


%Maintain Number of Columns
M = size(RangeProfiles_AfterEqNotch,1);                                         
 
%Maintain Number of Rows
N = size(RangeProfiles_AfterEqNotch,2);

%start bin - this variable is used to set the plotting axis at corresponding range bin 
start_bin = 400;

%define reference cell full reference cells = window_size * 2
Window_size = 18;
    
%Guard cells
Guard_cells = 4;

%Obtain the total length of Refence Cells
Reference_Cell = (Window_size) * 2;

%define desired Pfa
Pfa = 0.0000001;

%% Radar Detection
disp('Range Detection CA-CFAR')

%Data after Power Law 
Abs_Data = abs(RangeProfiles_AfterEqNotch).^2;

%Obtains the Data after power law data set size, number of rows and cloumn
dimensions = size(Abs_Data);

% sets a threshold array with initally as 0s with the same dimension of the data set
T_CACFAR = zeros(dimensions(1),dimensions(2));

%apply the CA_CFAR_function to column by column, the output returns an array of threshold values into the threshold array 
for i = 1:dimensions(1)
	T_CACFAR(i,:) = CA_CFAR_Function(Pfa, Window_size, Guard_cells, Abs_Data(i,:));
end

%detection logic, the data after power law array is substract with the threshold array
%(if negative, return 0 meaning no target;if positive return 1 meaning target detected)
Detections_rt = double((Abs_Data-T_CACFAR)>0);
%%----------------------------------------------------------------------------------------------------------------------------
%define PRI
PRI = 1/PRF_Hz;

%define light of speed
c = 3e8;

%convert from range bins to range (m)
range = 400*(c/(2*Bandwidth_Hz)):1*(c/(2*Bandwidth_Hz)):(NumOfRangeBins+399)*(c/(2*Bandwidth_Hz));

%convert num of pulse axis to time(s) axis
t= 0:PRI:PRI*M;

%plot time-domain range profiles
figure; axes('fontsize',fontsize1);
imagesc(range,t,20*log10(abs(RangeProfiles_AfterEqNotch)./MaxRangeLine),clims);
colorbar;
xlabel('Range(m)','fontsize',fontsize1);
ylabel('Time(s)','fontsize',fontsize1);
%title('Time-Domain Detection with X mark (single person walking)','fontsize',fontsize1);


%% marking the time-domain range profiles by plotting with 'x' on the target detected
x_plot_time_range(Detections_rt,RangeProfiles_AfterEqNotch,c,Bandwidth_Hz,PRI,start_bin);


%% plotting the range doppler map using CA-CFAR technique========================================================================= 

%define window_length
window_length = 512;                                                        % Dr Abdul Gaffar - must be a power of 2

%conver window_length to CPI
CoherentProcessingInterval_s = window_length*1/PRF_Hz;                      % Dr Abdul Gaffar

hop_length = window_length;

%count = 5;      % Dr Abdul Gaffar : limiting how many Range-Doppler images are formed 

%defines how many sectors the range profile is splitted into according to window_length
count = 1 + floor((size(RangeProfiles_AfterEqNotch,1)-window_length)/hop_length);
start = 1;


%%set a detection array with dimension of number of range-Doppler maps generated and range bins
Detection_Col = zeros(count,size(RangeProfiles_AfterEqNotch,2));

%plotting the range doppler map with x mark
for k = 1:count
    
    subset = RangeProfiles_AfterEqNotch(start:window_length+start-1,:);		%take the segment of full data
    hamm_window = repmat(hamming(window_length),1,size(subset,2));		%apply hamming window to supress
    window_set = subset.*hamm_window;						
    subset_fft = fftshift(fft(window_set,[],1),1);                              % Dr Abdul Gaffar ---YH comment: apply fft to convert to frequency domain
    %subset_fft = fftshift(fft(subset),1); 
    abs_subset = abs(subset_fft).^2;						% apply power law
    s_dimension = size(abs_subset);
    T_subset = zeros(s_dimension(1),s_dimension(2));				%set threshold array
    
    %apply to CA-CFAR function to get estimated threshold
    for i = 1:s_dimension(1)
        T_subset(i,:) = CA_CFAR_Function(Pfa, Window_size, Guard_cells, abs_subset(i,:));
    end
    
    n = size(RangeProfiles_AfterEqNotch,1);
    start = start + hop_length;  
    
    fs = (-window_length/2:1:(window_length/2-1))*PRF_Hz/window_length;         % Dr Abdul Gaffar ---YH comment: set up frequency axis
    
    %plot range-Doppler map   
     figure; 
     axes('fontsize',fontsize1);
     imagesc(range,fs,20*log10(abs(subset_fft)./MaxRangeLine),clims);                     % Dr Abdul Gaffar (removed axes)
     colorbar; 
     xlabel('Range(m)','fontsize',fontsize1);
     ylabel('Frequency (Hz)','fontsize',fontsize1);
%%     title('Range_Doppler Map',fontsize1);
    
    Detections_subset = double((abs_subset-T_subset)>0);			% apply detection logic
    
    %marking the plot with x
    x_plot_range_doppler(Detections_subset, abs_subset, c, Bandwidth_Hz, PRF_Hz, window_length,start_bin);
    
    Detection_Col(k,:) = sum(Detections_subset);                %put all detection of the range-Doppler map to 1 row according to the range bin
                
end


%% applying the detection from the doppler to the time domain

figure; axes('fontsize',fontsize1);
imagesc(range,t,20*log10(abs(RangeProfiles_AfterEqNotch)./MaxRangeLine),clims);
%imagesc(20*log10(abs(RangeProfiles_AfterEqNotch)./MaxRangeLine),clims);
colorbar;
xlabel('Range(m)','fontsize',fontsize1);
ylabel('Time(s)','fontsize',fontsize1);
%title('Range lines: after Eq Notch with X mark (Frequency)','fontsize',fontsize1);
hold on;

window_length_td = 512;
CoherentProcessingInterval_s_td = window_length_td*1/PRF_Hz; 
start_td = 1;
hop_length_td = window_length_td;

range_detection_plot(Detection_Col, count, start_td, hop_length_td,start_bin,c,Bandwidth_Hz,PRF_Hz);

%% plotting functions

function x_plot_pulses_rangebins(detection_arr, dt_array)

A = size(dt_array,1);   
B = size(dt_array,2); 
    
for a = 1:A   
    for b = 1:B   
        if detection_arr(a,b) > 0
            text(b,a,'X');
        end
    end
end

end

function x_plot_time_range(detection_arr, dt_array, c, bandwidth, PRI, start_bin)

M = size(dt_array,1);                                         
N = size(dt_array,2);

for m = 1:M   
    for n = start_bin:N+start_bin-1   
       if detection_arr(m,n-start_bin+1) > 0
            text(n*(c/(2*bandwidth)),m*PRI,'X');
        end
    end
end

end

function x_plot_range_doppler(detection_arr, dt_array, c, bandwidth, PRF, window_length,start_bin)

Y = size(dt_array,1);                                         
X = size(dt_array,2);

for y = 1:Y   
    for x = 1:X   
       if detection_arr(y,x) > 0
            text((x+start_bin-1)*(c/(2*bandwidth)),(y-(window_length/2)-1)*(PRF/window_length),'X');
        end
    end
end

end

function range_detection_plot(detection_arr, count, start_td, hop_length_td,start_bin,c,bandwidth,PRF)


for d = 1:count             %go through the rows of the detection_col
    for e = 1:size(detection_arr,2)             
        if detection_arr(d,e)>0
            plot((e+start_bin-1)*(c/(2*bandwidth)),(start_td:start_td + hop_length_td - 1)*(1/PRF),'kx','MarkerSize',12); %plots the full CPI column corresponding to range bins
        end
    end
    start_td = start_td + hop_length_td;
end

end




% %% Clear workspace and all figures
% clear all;clc;close all
% 
% %% Load Data
% [FileName,Path2RadarData,filter_index]=uigetfile('01  SingleHuman_50m_400M__07_12_2011_walking','D:\EEE4022S\Dataset\OneDrive_2021-08-23\07  Radar datasets');
% RadarData = load([Path2RadarData filesep FileName]);
% 
% %% Extract Range Profiles before, after Equalisation and after Notch filtering
% 
% RangeProfiles_BeforeEq = RadarData.RangeLines_BeforeEq;
% RangeProfiles_AfterEq = RadarData.RangeLines_AfterEq;
% RangeProfiles_AfterEqNotch = RadarData.RangeLines_AfterEQ_Notch;
% 
% %% Extract other radar parameters
% 
%  PRF_Hz = RadarData.Info.PRF_Hz;
%  Bandwidth_Hz = RadarData.Info.Bandwidth_Hz;
%  RangeStart_m = RadarData.Info.RangeStart_m;
%  BlindRange_m = RadarData.Info.BlindRange_m;
%  
% [NumOfPulses,NumOfRangeBins]=size(RangeProfiles_AfterEqNotch);
% 
% %% Plot Range Profiles
% 
% fontsize1 = 12;
% clims = [-40 0];
% 
% % Normalise data to have a peak of 0dB or 1 in linear scale
% [MaxRangeLine MaxIdx] = max(max(abs(RangeProfiles_AfterEqNotch)));
% 
% % Plot range lines
% figure; axes('fontsize',fontsize1);
% imagesc(20*log10(abs(RangeProfiles_AfterEqNotch)./MaxRangeLine),clims);
% colorbar;
% xlabel('Range (bins)','fontsize',fontsize1);
% ylabel('Number of pulses','fontsize',fontsize1);
% title('Range lines: after Eq Notch','fontsize',fontsize1);
% 
%  %Maintain Number of Columns
% M = size(RangeProfiles_AfterEqNotch,1);                                         
%  
% %Maintain Number of Rows
% N = size(RangeProfiles_AfterEqNotch,2);   
% 
% %define Training Cells
% Window_size = 16;
%     
% %Guard cells
% Guard_cells = 2;
% 
% %Obtain the total length of Refence Cells
% Reference_Cell = (Window_size) * 2;
% 
% %define desired Pfa
% Pfa = 0.00001;
% 
% %T_CACFAR = zeros(M,N);
% 
% %Data after Power Law 
% %Abs_Data = abs(RangeProfiles_AfterEqNotch).^2;
% 
% %alpha_CA = Reference_Cell*(Pfa^(-1/Reference_Cell)-1);
% 
% %% Radar Detection
% 
% disp('Range Detection')
% 
% %Data after Power Law 
% Abs_Data = abs(RangeProfiles_AfterEqNotch).^2;
% 
% dimensions = size(Abs_Data);
% 
% T_CACFAR = zeros(dimensions(1),dimensions(2));
% 
% for i = 1:dimensions(1)
% 	T_CACFAR(i,:) = CA_CFAR_Function(Pfa, Window_size, Guard_cells, Abs_Data(i,:));
% end
% 
% Detections_rt = double((Abs_Data-T_CACFAR)>0);
% 
% 
% % for j = 1:M
% %     RL_Data = Abs_Data(j,:);  % Ra_Data, size 534
% %     
% %     data_length = length(RL_Data);
% %         
% %     for i = 1:data_length
% %         if(i < Window_size+Guard_cells+1 || i > data_length - (Window_size+Guard_cells+1))
% %             T_CACFAR(i) = 0;
% %             continue;
% %         end
% %         CUT = i;
% %         
% %         RefCells_lagging = RL_Data(CUT-Window_size-Guard_cells:CUT-1-Guard_cells);
% %         RefCells_leading = RL_Data(CUT+1+Guard_cells:CUT+Window_size+Guard_cells);
% %         Sum_Reference_cells = sum(RefCells_lagging) + sum(RefCells_leading);  
% %         g = Sum_Reference_cells *(1/Reference_Cell);
% % 
% %         T_CACFAR(j,i) = alpha_CA * g;
% % 
% %     end    
% % 
% % end
% % 
% % Detections_rt = double((Abs_Data-T_CACFAR)>0);
% 
% 
% 
% PRI = 1/PRF_Hz;
% 
% c = 3e8;
% range = 0:1*(c/(2*Bandwidth_Hz)):NumOfRangeBins*(c/(2*Bandwidth_Hz));
% t= 0:PRI:PRI*M;
% 
% figure; axes('fontsize',fontsize1);
% imagesc(range,t,20*log10(abs(RangeProfiles_AfterEqNotch)./MaxRangeLine),clims);
% colorbar;
% xlabel('Range(m)','fontsize',fontsize1);
% ylabel('Time(s)','fontsize',fontsize1);
% title('Time Domain with X','fontsize',fontsize1);
% 
% 
% %marking the plot with x
% for y = 1:M   
%     for x = 1:N   
%        if Detections_rt(y,x) > 0
%             text(x*(c/(2*Bandwidth_Hz)),y*PRI,'X');
%         end
%     end
% end
% hold off
% 
% window_length = 256;
% hop_length = window_length;
% count = 1 + floor((size(RangeProfiles_AfterEqNotch,1)-window_length)/hop_length);
% start = 1;
% 
% 
% for k = 1:count
%     
%     subset = RangeProfiles_AfterEqNotch(start:window_length+start-1,:);
%     subset_fft = fft(subset,[],1);
%     abs_subset = abs(subset_fft).^2;
%     s_dimension = size(abs_subset);
%     T_subset = zeros(s_dimension(1),s_dimension(2));
%     for i = 1:s_dimension(1)
%         T_subset(i,:) = CA_CFAR_Function(Pfa, Window_size, Guard_cells, abs_subset(i,:));
%     end
% 
%     n = size(RangeProfiles_AfterEqNotch,1);
%     
%     fs = PRF_Hz*(0:1:256)/n;
% 
%     figure; 
%     axes('fontsize',fontsize1);
%     imagesc(range,fs,20*log10(abs(subset_fft)./MaxRangeLine),clims);
%     colorbar;
%     xlabel('Range(m)','fontsize',fontsize1);
%     ylabel('Frequency','fontsize',fontsize1);
%     title('Range_Doppler Map',fontsize1);
%     start = start + hop_length;
%     A = size(abs_subset,1);   
%     B = size(abs_subset,2); 
%     
%     Detections_rt = double((abs_subset-T_subset)>0);
%     
%     %marking the plot with x
%     for a = 1:A   
%         for b = 1:B   
%             if Detections_rt(a,b) > 0
%                 text(b*(c/(2*Bandwidth_Hz)),a*PRF_Hz/n,'X');
%             end
%         end
%     end
%     hold off
% 
%     
% end







