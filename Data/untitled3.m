%JasonSecula
%   Muon Lab Report
%   12/5/16
name = 'Jason Secula';
disp(name)
lab = 'Muon Decay Lab';
disp(lab)
%Trial_Two is an array with all the data organized from least to greatest
subplot(3,3,1)
plot(Trial_Two)
title('Line Plot of each decay event recorded during trial two')
xlabel('Time in Seconds Mirco-seconds (X 10^-6)'), ylabel('Number of Decay event'), grid
legend('Plot Trial Two Variable')
%Trial_Two is an array with all the data organized from least to greatest
%NumTrial Describes an array from 1 to 2862 (Number of events)
subplot(3,3,2)
plot(Trial_Two, NumTrial)
title('Line Plot of each decay event recorded during trial two')
xlabel('Time in Seconds Mirco-seconds (X 10^-6)'), ylabel('Number of Decay event'), grid
legend('Plot Trial 2 vs. 1 through 2826 Variable')
%Liner scatter plot of Trial two graphed against itself
subplot(3,3,3)
scatter(Trial_Two, Trial_Two)
title('Number of Muon Decay event N(t).')
xlabel('Time in Seconds Mirco-seconds (X 10^-6)'), ylabel('Number of Decay event'), grid
legend('Trial2 vs. Trial2')
%Scatter plot of Ln(t) vs. t
subplot(3,3,4)
scatter(lnNt, Trial_Two)
title('Number of Muon Decay event Ln(t) vs t.')
xlabel('Time in Seconds Mirco-seconds (X 10^-6)'), ylabel('Number of Decay event'), grid
legend(' Ln(t) vs. Time')
%------------Trial Two Histogram----------------------------------
subplot(3,3,5)                              % 3 graphs tall, 3 wide.
h1 = histogram(Trial_Two);                  %In spot 5, a histogram of trial two is plotted.
title(' Number of Muon Decay event N(t). ') %Graph Details.
xlabel('Time in Micro-seconds (X 10^-6)'), ylabel('Number of Decays events'), grid
legend('Histogram of Trial 2')
hold on
h1.BinWidth = 0.2;              %Currently the Bin Width is set to 0.2MuSeconds thick
h1.FaceColor = [0.7 0.7 0.7];   %Displays a silver color
h1.EdgeColor = 'b';             %'b'-Blue outline.
%------------------
% coeffs = polyfit(Trial_Two,NumTrial,1);
% 
% Coefficients = table(coeffs, 'VariableNames',{'Coefficient'}) %#ok<NOPTS>
% 
% %Get fitted values
% fittedX = linspace(min(Trial_Two), max(Trial_Two), 200);
% fittedY = polyval(coeffs, fittedX);
% %Plot the fitted line
% hold on;
% plot(fittedX, fittedY, 'r-', 'LineWidth', 3);
%----------------------
subplot(3,3,9)
h1.Normalization = 'probability';
disp(h1.Normalization)
h1.BinWidth = 0.2;
ban = 'Bandwith is';
disp(h1.BinWidth)
hold on
binranges = 0.012:19.96;
[bincounts] = histc(Trial_Two,binranges);
br = 'Bins range frm 0.012muSec to 19.96musec';                   %BinRange
disp(br)                                                        %BinRange
disp(binranges)                                                 %BinRange
 
bc = 'Number of Decay per Bin:';                                %BinCount
disp(bc)                                                        %BinCount
disp(bincounts)                                                 %BinCount
%--------------------Background Histogram-------------------------
BackGrnd1 = bincounts - 143.15;             %20Bins, default by matLab
BackGrnd2 = bincounts - 71.575;             %40Bins, 0.5 width
BackGrnd3 = bincounts - 28.63;              %100Bins, 0.2 width
bg = 'Background = number of decays/number of bins';
disp(bg)    
BackGround = table(BackGrnd1,BackGrnd2,BackGrnd3,'VariableNames',{'Default','FortyBin','WunHunitdBin'})%#ok<NOPTS>
 
subplot(3,3,8)
h3 = histogram(BackGrnd1);
hold on
h4 = histogram(BackGrnd2);
hold on
h5 = histogram(BackGrnd3);
hold on
title(' Number of Muon Decay event N(t) between 0.12 and 20 microseconds ')
xlabel('Time in Micro-seconds (X 10^-6)'), ylabel('Number of Decays events'), grid
legend('Back Ground of 2863 decay events/20bins', 'Back Ground of 2863 decay events/40bins','Back Ground of 2863 decay events/100bins')
 
subplot(3,3,7)
h6 = histogram(BackGrnd);
title(' No. of Muon Decay event N(t) between 0.12 and 20 microsec ')
xlabel('Time in Micro-seconds (X 10^-6)'), ylabel('Number of Decays events'), grid
% ------------All the Histogram----------------------------------
subplot(3,3,6)
h7 = histogram(Trial_Two);
hold on
h9 = histogram(BackGrnd);
hold on
title(' All the grams. ')
xlabel('Time in Micro-seconds (X 10^-6)'), ylabel('Number of Decays events'), grid
legend('All of the Histograms')
average_Decay = mean(BackGrnd);               %Average
sa1 = 'Average Decay for values greater than 12microseconds hours. '; 
disp(sa1)                                     %Average
disp( average_Decay )                       %Average
DataCollected = table(average_Decay2,standard_deviation2,maximum_Decay2,minimum_Decay2,'VariableNames',{'AverageDecay ','StandardDeviation','MaximiumDecay','MinimumDecay'}) %#ok<NOPTS>
diary ('MuonTest1')                         %Save File
%whos                                       %FileSize
win = 'complete';
disp(lab)
disp(win)

