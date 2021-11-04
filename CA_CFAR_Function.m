%Function For CA-CFAR, Taking in Values of Pfa, Window Size, Guard Cells,
%and the signal,outputing the Threshold
% 
% function [T_CFAR] = CA_CFAR_Function(Pfa,WinSize,GuardCell,RL_Data)
% 
% Data_length = length(RL_Data);
% 
% RefCells = WinSize * 2;
% 
% StartCUT = RefCells + 1;
% StopCUT = Data_length - RefCells;
% 
% T_total = zeros(1,Data_length);
% CUT_Total = zeros(1,Data_length);
% 
% DataAfterPowerLawDetector = RL_Data;
% 
% for i = StartCUT:StopCUT  % Start at the index of 33 and end at an index of before CUT + RefWindow > data length
%     
%     CUT = i;   % Dr Abdul Gaffar. For now, only considering CUT = 50. %YH: set the CUT to be i and it would slide to the next cell of interest
%     
%     RefCells_lagging = DataAfterPowerLawDetector( (CUT - WinSize- GuardCell):(CUT-1-GuardCell));  % Total length = RefWindow/2
%     RefCells_leading = DataAfterPowerLawDetector( (CUT+1+GuardCell):(CUT + WinSize));  % Total length = RefWindow/2
%     Sum_Reference_cells = sum(RefCells_lagging) + sum(RefCells_leading);    
%     
%     g = Sum_Reference_cells./RefCells;
% 
%     %Scaling factor
%     alpha = RefCells*(Pfa^(-1/(RefCells))-1);
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
% FalseAlarm = (CUT_Vector > Threshold_vector); 
% NumOfFalseAlarm = sum(FalseAlarm);
% 
% PFA_SimulatedData = NumOfFalseAlarm/(length(CUT_Vector));
% PFA_error = abs((Pfa - PFA_SimulatedData)/Pfa)*100
% 
% T_CFAR = Threshold_vector;
% 
% end


function [T_CA_CFAR] = CA_CFAR_Function(Pfa, Win_Size, Guard_Cells, Abs_Signal)

RefWin = (Win_Size)*2;
Data_length = length(Abs_Signal);

T_CFAR = zeros(Data_length,1);
for i = 1:Data_length
    CUT = i;
	if(CUT < Win_Size+Guard_Cells+1 || CUT > Data_length - (Win_Size+Guard_Cells+1))
		T_CFAR(CUT) = 1;
		continue;
	end
	RefCells_lagging = Abs_Signal(CUT-Win_Size-Guard_Cells:CUT-1-Guard_Cells);
	RefCells_leading = Abs_Signal(CUT+1+Guard_Cells:CUT+Win_Size+Guard_Cells);
	
	Sum_Reference_cells = sum(RefCells_lagging)+sum(RefCells_leading);
	g = Sum_Reference_cells./RefWin;
	alpha = RefWin*(Pfa^(-1/(RefWin))-1);
	T_CFAR(CUT) = g*alpha;
	
end

T_CA_CFAR = T_CFAR;

end


