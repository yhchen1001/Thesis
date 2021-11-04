% clear all;
% close all;
% 
% % Input parameters
% Pfa = 10^-3;
% Data_length = 100000;
% RefWindow = 32;
% Guard_cells = 1;
% 
% RefCells = zeros(RefWindow,Data_length);
% 
% % Generate simulated data
% %y_complex  = randn(1,Data_length) + 1i*randn(1,Data_length);
% 
% y_complex  = (normrnd(0,1,1,Data_length) + 1i*normrnd(0,1,1,Data_length))/sqrt(2);
% %power law detector
% DataAfterPowerLawDetector = abs(y_complex).^2;
% Detection = zeros(1,Data_length);
% T_total = zeros(1,Data_length);
% CUT_Total = zeros(1,Data_length);
% 
% StartCUT = 33;
% StopCUT = Data_length-RefWindow;
% 
% for i = StartCUT:StopCUT  % Start at the index of 33 and end at an index of before CUT + RefWindow > data length
%     
%     CUT = i;   % Dr Abdul Gaffar. For now, only considering CUT = 50. %YH: set the CUT to be i and it would slide to the next cell of interest
%     
%     RefCells_lagging = DataAfterPowerLawDetector( (CUT-RefWindow/2-Guard_cells):(CUT-1-Guard_cells));  % Total length = RefWindow/2
%     RefCells_leading = DataAfterPowerLawDetector( (CUT+1+Guard_cells):(CUT+RefWindow/2+Guard_cells));  % Total length = RefWindow/2
%     Sum_Reference_cells = sum(RefCells_lagging) + sum(RefCells_leading);    
%     
%     g = Sum_Reference_cells./RefWindow;
% 
%     %Scaling factor
%     alpha = RefWindow*(Pfa^(-1/(RefWindow))-1);
% 
%     % Compute threshold
%     T = g*alpha;
%     T_total(i) = T;
%     
% 
%     CUT_Value = DataAfterPowerLawDetector(CUT);
%     
%     CUT_Total(i) =  CUT_Value;
%     
% end
% 
% CUT_Vector = DataAfterPowerLawDetector(StartCUT:StopCUT);
% Threshold_vector = T_total(StartCUT:StopCUT);
% 
% FalseAlarm = (CUT_Vector > Threshold_vector);             % A FalseAlarm value of 0 means a false target was not detected
% NumOfFalseAlarm = sum(FalseAlarm);                        % A FalseAlarm value of 1 means that a false target was detected
% PFA_SimulatedData = NumOfFalseAlarm/(length(CUT_Vector));
% 
% PFA_error = abs((Pfa - PFA_SimulatedData)/Pfa)*100

% X = [1 2 3; 4 5 6; 7 8 9]; 
% figure; 
% imagesc(X); 
% colorbar; 
% hold on; 
% plot(1,1, 'kx', 'MarkerSize',12); 
% plot(2,3,'kx', 'MarkerSize',12)



% clear all;
% close all;
% 
% % Input parameters
% Pfa = 10^-3;
% % Data_length = 100000;
% RefWindow = 32;
% Guard_cells = 1;
% 
% %RefCells = zeros(RefWindow,Data_length);
% 
% % Generate 2-D simulated data
% %y_complex  = randn(1,Data_length) + 1i*randn(1,Data_length);
% 
% row = 1000;
% column = 100;
% NumOfPoints = row*column;
% 
% y_complex  = (normrnd(0,1,row,column) + 1i*normrnd(0,1,row,column))/sqrt(2);
% %power law detector
% DataAfterPowerLawDetector = abs(y_complex).^2;
% T_total = zeros(row,column);
% CUT_Total = zeros(row,column);
% 
% StartCUT = 33;
% StopCUT = column-RefWindow;
% 
% for i = StartCUT:StopCUT  % Start at the index of 33 and end at an index of before CUT + RefWindow > data length
%     
%     CUT = i;   % Dr Abdul Gaffar. For now, only considering CUT = 50. %YH: set the CUT to be i and it would slide to the next cell of interest
%     
%     RefCells_lagging = DataAfterPowerLawDetector( (CUT-RefWindow/2-Guard_cells):(CUT-1-Guard_cells));  % Total length = RefWindow/2
%     RefCells_leading = DataAfterPowerLawDetector( (CUT+1+Guard_cells):(CUT+RefWindow/2+Guard_cells));  % Total length = RefWindow/2
%     Sum_Reference_cells = sum(RefCells_lagging) + sum(RefCells_leading);    
%     
%     g = Sum_Reference_cells./RefWindow;
% 
%     %Scaling factor
%     alpha = RefWindow*(Pfa^(-1/(RefWindow))-1);
% 
%     % Compute threshold
%     T = g*alpha;
%     T_total(i) = T;
%     
% 
%     CUT_Value = DataAfterPowerLawDetector(CUT);
%     
%     CUT_Total(i) =  CUT_Value;
%     
% end
% 
% CUT_Vector = DataAfterPowerLawDetector(StartCUT:StopCUT);
% Threshold_vector = T_total(StartCUT:StopCUT);
% 
% FalseAlarm = (CUT_Vector > Threshold_vector);             % A FalseAlarm value of 0 means a false target was not detected
% NumOfFalseAlarm = sum(FalseAlarm);                        % A FalseAlarm value of 1 means that a false target was detected
% PFA_SimulatedData = NumOfFalseAlarm/(length(CUT_Vector));
% 
% PFA_error = abs((Pfa - PFA_SimulatedData)/Pfa)*100


% %% Cell Averaging CFAR Detection Algorithm
% clear all;  
% close all;
% 
% PFA = 1e-3; % Abdul Gaffar 1e-3 = 1/1000             
%  
% row = 1000;
% column = 100;
% 
% a = normrnd(0,1,row,column); %Chosen by Gaussian PDF
% b = normrnd(0,1,row,column); %Chosen by Gaussian PDF
% x = (a + (1i*b))/sqrt(2); %a - I, b - Q --> NOISE SO all detections are false alarms??
% 
% SLD = abs(x).^2; % Square Law Detector
%   
%     %[row, column] = size(SLD);
%     NumDataPoints = row*column;
%     
%     %% Arrays to store detection positions
%     row_det = []; %detection postiions in row
%     column_det = []; %detection postiions in column
%     
%     %% Parameters
%     window = 32; %window size
%     N = window*2; %Number of Reference Cells
%     guard_cells = 2; %Guard cells - 3/4
%     
%     CFAR_T =  []; % intialise array for each CUT
%     counter = 0;
%     
%     %% Determine CFAR for each CUT
%     %% Apply algorithm along each column in spectrogram
%     for c = 1:column
%         
%         power = SLD(1:row,c); %Single column
%         
%         %No false alarm outside the window and guard cell region
%         region = window + guard_cells +1;
%         for r = region:row-region
%             
%             CUT = power(r); %power of CUT
%             %CA-CFAR
%             lag_window = power(r-window-guard_cells:r-1-guard_cells);
%             lead_window = power(r+1+guard_cells:r+window+guard_cells);
%             
%             %Calculate CA-CFAR Threshold
%             %1. Calculate the interference statistic
%             g_CA = (sum(lag_window)+ sum(lead_window))./N;
%             %2. Calculate CFAR constant
%             alpha_CA = N*(PFA^(-1/(N))-1);
%             %CA-CFAR Threshold
%             threshold = g_CA*alpha_CA;
%             
%             CFAR_T = [CFAR_T;threshold];
%             
%             if (threshold < CUT)
%                 row_det = [row_det; r];
%                 column_det = [column_det;c];
%                 counter = counter + 1;
%             end 
% 
%         end
%   
%     end
%     
%     %% Plot r(CUT) and threshold - 100000*10e-3 = 100 false alarms
%     PFA_sim = counter/NumDataPoints;
%     Error = abs(((PFA_sim-PFA)/PFA)*100) % Abdul Gaffar abs()



clear all;
close all;

Pfa = 1e-4;


Iterations = 1e6;
N = 24;


Reference_cells = zeros(N,Iterations);

for i = 1:N
    I = randn(1,Iterations);
    Q = randn(1,Iterations);
    Reference_cells(i,:) = (I + 1j*Q)/sqrt(2);
end

Reference_cells_AD = abs(Reference_cells).^2;
Sum_Reference_cells = sum(Reference_cells_AD,1);
g = Sum_Reference_cells./N;
%Scaling factor
alpha = N*(Pfa^(-1/(N))-1);


%CA Threshold
T = g*alpha;

I_test = randn(Iterations,1);
Q_test = randn(Iterations,1);
Test_Signal = (I_test + 1j*Q_test)/sqrt(2);
Test_AD = abs(Test_Signal).^2;

False_Alarms = sum((Test_AD.'-T)>0);
Simulated_Pfa = False_Alarms/Iterations

Pfa_Error = abs(100*(Simulated_Pfa-Pfa)/Pfa)


