%JasonSecula
%   Muon Lab Report
%   12/5/16
%-------------------------------------------------------------------
%Variables
d1 = '-----------------------------------------------------------------';
win = 'complete';
%Trial_Two is the raw data: events of decay
name = 'Jason Secula';
disp(name)
lab = 'Muon Decay Lab';
disp(lab)
disp(d1)
mario = 'Okay, here we go!';
disp(mario)
disp(d1)
%-----------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------
%------------------------------Trial Two Histogram----------------------------------
%-----------------------------------------------------------------------------------
figure('Name','Trial_Two, Defaulted by MatLab')                    % Figure 3
h1 = histogram(Trial_Two);                  %In spot 5, a histogram of trial two is plotted.
%-------------------------
title(' Number of Muon Decay event N(t). ') %Graph Details.
xlabel('Time in Micro-seconds (X 10^-6)'), ylabel('Number of Decays events'), grid
legend({'Histogram of Trial 2','BinRange vs. BinWidth'})
disp(win)
disp(d1)
hold on
%-----------------------------------------------------------------------------------
%-----------------------------Normalized, Width 0.2---------------------------------
%-----------------------------------------------------------------------------------
figure('Name','Trial_Two, Banwidth of 0.2MicroSeconds')              % Figure 4
h2 = histogram(Trial_Two); 
h2.BinWidth = 0.2;              %Currently the Bin Width is set to 0.2MuSeconds thick
h2.FaceColor = [0.8 0.8 0.8];   
h2.EdgeColor = 'b';             %'b'-Blue outline.
hold on
disp(win)
disp(d1)
ban = 'Bandwith is';
disp(ban)
disp(h2.BinWidth)
disp(d1)
hold on
%----------------------
BinRng1 = 0.012:0.2:19.96;        %start, interval (0.2microseconds), end
[BinCnts,ind] = histc(Trial_Two,BinRng1);
br = 'Bins range frm 0.012muSec to 19.96musec';                   %BinRange
disp(br)                                                        %BinRange
disp(BinRng1)                                                 %BinRange
 
bc = 'Number of Decay per Bin at 0.2Width: ';                                %BinCount
disp(bc)                                                        %BinCount
disp(BinCnts)                                                 %BinCount
disp(d1)
hold off
%-----------------------------------------------------------------------------------
%-----------------------------BackGround, Width 0.2---------------------------------
%-----------------------------------------------------------------------------------
figure('Name','BackGround for 0.2 width, 100 bin')     
p1 = polyfit(PntTwoCnt,PntTwoRng,1);
disp(p1)
Coefficients1 = table(p1, 'VariableNames',{'CountVsRange1'}) %#ok<NOPTS>
disp(d1)
 
[p2,S,mu] = polyfit(PntTwoCnt,PntTwoRng,1);
disp(p2)
p2S = 'PolyFit Range vs. Bin Count, 0.2MicroSeconds version 2';
Coefficients2 = table(p2, 'VariableNames',{'CountVsRange2'}) %#ok<NOPTS>
disp(d1)
disp(S)
disp(mu)
p2mu = 'PolyFit Range vs. Bin Count, 0.2MicroSeconds version 2';
disp(p2mu)
 
Ply = table(p1,p2,'VariableNames',{' PolyFit 1 ',' Poly Fit 2'}) %#ok<NOPTS>
disp(d1)
 
%--------------------------------
%------------PLOT----------------
%--------------------------------
fittedX = linspace(min(Trial_Two), max(Trial_Two),100) %#ok<NOPTS>
fittedY = polyval(p1, fittedX);
fitx = 'Fitted X coordinates';
disp(fitx)
disp(fittedX)                                                 %BinCount
disp(d1)
fity = 'Fitted Y coordinates';
disp(fity)
disp(fittedY)                                                 %BinCount
disp(d1)
%Plot the fitted line
hold on;
 
subplot(2,2,1)
bar(h2.');
hold on;
plot(fittedX, fittedY, 'or', 'LineWidth', 3);
 
fittedX2 = linspace(min(Trial_Two), max(Trial_Two),100) %#ok<NOPTS>
fittedY2 = polyval(p2, fittedX);
fitx = 'Fitted X coordinates';
disp(fitx)
disp(fittedX2)                                                 %BinCount
disp(d1)
fity = 'Fitted Y coordinates';
disp(fity)
disp(fittedY2)                                                 %BinCount
disp(d1)
%Plot the fitted line
hold on;
 
subplot(2,2,2)
bar(h2.');
hold on;
plot(fittedX, fittedY, 'or', 'LineWidth', 3);
 
hold off
%--------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------
%-----------------------------Normalized, Width0.5 ------------------------------------
%--------------------------------------------------------------------------------------
figure('Name','Trial_Two, Banwidth of 0.5MicroSeconds')          % Figure 5
h3 = histogram(Trial_Two); 
h3.BinWidth = 0.5;              %Currently the Bin Width is set to 0.5MuSeconds thick
h3.FaceColor = [0.5 0.5 0.5];   
h3.EdgeColor = 'k';             %'k'-Blue outline.
hold on 
 
h3.Normalization = 'probability';
disp(h3.Normalization)
 
h3.BinWidth = 0.5;
ban = 'Bandwith is';
disp(ban)
disp(h3.BinWidth)
disp(d1)
hold on
%----------------------
BinRng1 = 0.012:0.5:19.96;        %start, interval (0.2microseconds), end
[bincounts] = histc(Trial_Two,BinRng1);
br = 'Bins range frm 0.012muSec to 19.96musec';                   %BinRange
disp(br)                                                        %BinRange
disp(BinRng1)                                                 %BinRange
 
bc = 'Number of Decay per Bin at 0.5Width: ';                   %BinCount
disp(bc)                                                        %BinCount
disp(bincounts)                                                 %BinCount
disp(d1)
legend('show')
hold off
 
 
 
 
 
 
 
 
 
 
 
 
% ------------------
disp( average_Decay )                       %Average
disp(d1)
%-----------------------------------------------------------------
DataCollected = table(average_Decay2,standard_deviation2,maximum_Decay2,minimum_Decay2,'VariableNames',{'AverageDecay ','StandardDeviation','MaximiumDecay','MinimumDecay'}) %#ok<NOPTS>
%-----------------------------------------------------------------
diary ('MuonHistoT2v9')                       %Save File
%whos                                       %FileSize
%-----------------------------------------------------------------
disp(d1)
disp(lab)
disp(win)
disp(d1)
%----------------------------------------------------------------
