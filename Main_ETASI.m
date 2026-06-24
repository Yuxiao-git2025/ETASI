%% Loading-Data
% It is better to analyze the aftershock data (like Ridgecrest or Ding'ri)
clc;
[time,mag,Mmin]=LoadGansu();
Param0 = [0.1, 0.005, 1.6, 1.2, 1.5, 1.2, 0.01];
Result = cal_ETASI(time, mag, Mmin, Param0);
%% Plot Results
ETASI_PlotResult(Result);