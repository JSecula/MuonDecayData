%JasonSecula
%   Muon Lab Report
%   12/5/16
%-------------------------------------------------------------------
%Variables
d1 = '-----------------------------------------------------------------';
win = 'completed';
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
%------------------------------Trial Two Histogram----------------------------------
%-----------------------------------------------------------------------------------
figure('Name','Trial_ Two, Defaulted by MatLab')                    % Figure 1
h1 = histogram(Trial_Two);                  %In spot 5, a histogram of trial two is plotted.
%-------------------------
title(' Number of Muon Decay event N(t). ') %Graph Details.
h1.FaceColor = [0 0 1];   
h1.EdgeColor = [0 0 1];     
xlabel('Time in Micro-seconds (X 10^-6)'), ylabel('Number of Decays events'), grid
legend('show')
%-----------------------------------------------------------------------------------
%-----------------------------Normalized, Width 0.2---------------------------------
%-----------------------------------------------------------------------------------
figure('Name','Trial_ Two, Banwidth of 0.2MicroSeconds')              % Figure 2
subplot(2,2,2)
h2 = histogram(Trial_Two); 
h2.BinWidth = 0.2;              %Currently the Bin Width is set to 0.2MuSeconds thick
h2.FaceColor = [0 0.7 0.5];   
h2.EdgeColor = [0 0.7 0.5];             %
 
BinRng1 = 0.012:0.2:19.96;        %start, interval (0.2microseconds), end
[BinCnts,ind] = histc(Trial_Two,BinRng1);
br = 'Bins range frm 0.012muSec to 19.96musec';                   %BinRange
disp(br)                                                        %BinRange
disp(BinRng1)                                                 %BinRange
 
bc = 'Number of Decay per Bin at 0.2Width: ';                                %BinCount
disp(bc)                                                        %BinCount
disp(BinCnts)                                                 %BinCount
disp(d1)
hold on
 
title(' No. of Muon Decay event N(t) between 0.12 and 20 microsec ')
xlabel('Time in Micro-seconds (X 10^-6)'), ylabel('Number of Decays events'), grid
legend('show')
 
%-----------------------------------------------------------------------------------
%-----------------------------BackGround, Width 0.2---------------------------------
%-----------------------------------------------------------------------------------
subplot(2,2,1)
p1 = plot(PntTwoCnt, PntTwoRng,'Color', [0.3 0.5 0.6]);
subplot(3,2,1)
plot(PntTwoCnt, PntTwoRng, 'Color', [1 0 0]);
minus = BinCnts - (PntTwoCnt/PntTwoRng);
hold on
title(' Minus Background')
xlabel('Time in Micro-seconds (X 10^-6)'), ylabel('Number of Decays events'), grid
legend('show')
 
subplot(2,2,3)
h4 = histogram(minus); 
h4.BinWidth = 0.2;              %Currently the Bin Width is set to 0.2MuSeconds thick
h4.FaceColor = [0.5 0.5 0.5];   
h4.EdgeColor = [0.5 0.5 0.5];             %
hold on
title(' Count (rise) vs. Range (run) ')
xlabel('Time in Micro-seconds (X 10^-6)'), ylabel('Number of Decays events'), grid
legend('show')
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
hold on
%----------------------
BinRng2 = 0.012:0.5:19.96;        %start, interval (0.5microseconds), end
[BinCount] = histc(Trial_Two,BinRng2);                 %BinRange
disp(br)                                                        %BinRange
disp(BinRng2)                                                 %BinRange
 
bc = 'Number of Decay per Bin at 0.5Width: ';                   %BinCount
disp(bc)                                                        %BinCount
disp(BinCount)                                                 %BinCount
title(' No. of Muon Decay event N(t) between 0.12 and 20 microsec ')
xlabel('Time in Micro-seconds (X 10^-6)'), ylabel('Number of Decays events'), grid
legend('show')
 
%-----------------------------------------------------------------------------------
%-----------------------------BackGround, Width 0.5---------------------------------
%-----------------------------------------------------------------------------------
hold on
subplot(2,2,4)
p3 = plot(PntFiveCnt, PntFiveRng);
title(' Count (rise) vs. Range (run) ')
xlabel('Time in Micro-seconds (X 10^-6)'), ylabel('Number of Decays events'), grid
legend('show')
hold on
 
subplot(2,2,1)
plot(PntFiveCnt, PntFiveRng, 'Color', [1 0 0]);
minus = BinCount - PntFiveCnt/PntFiveRng;
title('The Minus BackGround Histogram ')
xlabel('Time in Micro-seconds (X 10^-6)'), ylabel('Number of Decays events'), grid
legend('show')
hold on
 
 
subplot(2,2,3)
h5 = histogram(minus); 
h5.BinWidth = 0.5;              %Currently the Bin Width is set to 0.5MuSeconds thick
h5.FaceColor = [0.5 0.5 0.5];   
h5.EdgeColor = [0.5 0.5 0.5];    
xlabel('Time in Micro-seconds (X 10^-6)'), ylabel('Number of Decays events'), grid
legend('show')
title(' Minus Back Ground ')
 
hold off
%-----------------------------------------------------------------------------------
%------------------------------------------MATH-------------------------------------
%-----------------------------------------------------------------------------------
 
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
