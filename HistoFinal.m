%JasonSecula
%   Muon Lab Report
%   12/5/16
%--------------------------------------------------------------------------------------------------------------------------------------------------------------
%Variables
%--------------------------------------------------------------------------------------------------------------------------------------------------------------
d1 = '--------------------------------------------------------------------------------------------------------------';
win = 'Complete';
name = 'Jason Secula';
lab = 'Muon Decay Lab';
str = date;
mario = 'Okay, here we go!';

%----------------
disp(name)      %           First Object displayed
disp(lab)       %           Name, Lab, Date, Kick off
disp(str)       %
disp(d1)        %
disp(mario)     %
disp(d1)        %
%----------------

Expected = 2.19703;
NoD = 2862;                                 %number of decays
SumDecay = sum(MuD);                        %Decay's added (1.2015e+04)
SrtD = sqrt(NoD);                           %SquareRoot of the Number of Decay for uncertainty

Mu = mean(MuD);                             %Average
yln = log(MuD);                             %this would be ln(t) not ln(N(t)) ???
sigma = std(MuD);                           %StandardD (4.986757523444430)
maximumD = max(MuD);                        %Max
minimumD = min(MuD);                        %Min
ProbDisthD = fitdist(MuD, 'Normal');        %Normal Fit MuonDecay Events
llambda = (NoD/SumDecay);                   %(llambda = 0.2382)

%--------------------------------------------------------------------------------------------------------------------------------------------------------------
muDecayEvents = table(Mu,sigma,minimumD,maximumD,'VariableNames',{'AverageDecay ','StandardDeviation','MinimumDecay','MaximiumDecay'}) %#ok<NASGU>
disp(d1)
%--------------------------------------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------Default Histogram-------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------------------------------------------------------

%Variables for Trial Two Histogram with Default Bin Width. 
DeltaNhd = NoD/20;                        %Change in N = NumberOfCounts/BinCount
SrtBhd = sqrt(BinCnt);                    %SquareRoot of the Number of Decay for uncertainty
BGhD = 646.66;                          %Background = (Sum DecayEvent from 12 to 20MuSec)/(Number of Bins 12 to 20 MuS)
BinRD = 0.04:1:19.96;                %Bin Range
betahd0 = 1;
rzX = [12.54;13.54;14.54;15.54;16.54;17.54;18.54;19.54]; %BinCenters (12 to 20 microSec)
rzY = [41;44;42;39;34;34;54;36]; %BinRange(12 to 20 microSec)
BinCounthD = [897;529;319;203;168;96;80;64;56;46;47;35;40;44;42;42;30;36;54;0];
%--------------------------------------------------------------------------------------------------------------------------------------------------------------
%Figure 1)
figure('Name','Muon Decay, Defaulted by MatLab')

hD = histogram(MuD);           %#ok<*NOPTS>
hD.FaceColor = [0.7 0.7 0.7];   
hD.EdgeColor = [0 0 0];  
BinCnt = histcounts(MuD);
        %[counts,centers] = histogram-bin counts of data set(MuD),
                    %with range BinRD
BinRD = BinRD(:);               %(:) converts Bin Range to column
BinCnt = BinCnt(:);             %(:) converts Bin Count into a column
CentrBn = (BinRD+0.5)             %Center of Each Bin

hold on
title(' No. of Muon Decay event N(t) between 0.12 and 20 microsec ')
xlabel('Time in Micro-seconds (X 10^-6)')
grid
hold on

Xp1 = CentrBn(:,1); %Interval the DecayEventData set down collum one
Yp1 = BinCounthD(:,1); %Interval the Bin Range set down collum one
format long
LLambdap1 = Xp1\Yp1;       %llambda is slope
yCalcp1 = LLambdap1*Xp1;

% p1 =plot(Xp1,Yp1,':p','LineWidth',1,'MarkerSize',3);
% hold on

%-----------------------------------------------------------------------------------
%-----------------------------------------BackGround of Default Histogram-----------
%-----------------------------------------------------------------------------------
ylabel('Number of Decays events')
p2 = bar(rzX,rzY);
p2.FaceColor = [0 0 1];   
p2.EdgeColor = [0 0 1];

legend('Histogram of MuonDecay Events',' Center of Bin vs. BinCound')
hold off

%Figure 2)
figure('Name','Non-Linear Least Square Fit')
LeastSqrtFitTot = lsqnonneg(CentrBn,BinCounthD);   %(4.0319)
LeastSqrtFitBG = lsqnonneg(rzX,rzY);             %(2.4711)

FittMLS = times(MuD,LeastSqrtFitTot);
FittMLSBG = times(MuD,LeastSqrtFitBG);

p3 = histogram(FittMLS);
p3.FaceColor = [0.5 0.9 0.2 ];   
p3.EdgeColor = [0 0 0];
hold on

p4 = histogram(FittMLSBG);
p4.FaceColor = [1 1 0.7];   
p4.EdgeColor = [0 0 0];
hold on

legend([p3 p4],'NonLinear Square Fit','BackGround NonLinear Square Fit')
hold off

%--------------------------------------------------------------------------------------------------------------------------------------------------------------
[Xpos,resnorm,residual] = lsqnonneg(CentrBn,BinCounthD); %Obtain the solution and residual information.
  NormTest = norm(residual)^2;   %^LineSquareNonNegativeFunction
        %^Test Residual Value for Validity        
%Nonnegative Least Squares with Nondefault Options
options = optimset('Display','final');

%Figure 3)
figure;
hD1 = histogram(MuD);
hold on
hDr = bar(CentrBn, residual);     
hDr.FaceColor = [0.3 0.8 1];   
hDr.EdgeColor = [0 0 0];  

title(' (Chi^2) : residual Muon Decay event N(t) between 0.12 and 20 microsec ')
xlabel('Time in Micro-seconds (X 10^-6)')
ylabel('Number of Decay Events')
grid
legend([hD1 hDr], 'Muon Decay Events', 'Chi^2')
hold on

%--------------------------------------------------------------------------------------------------------------------------------------------------------------
%Figure 4)
figure('Name','MuonDecayEvents Fitted to Line'); 

hhD2 = histogram(MuD); %#ok<*NOPTS> 
hhD2.BinWidth = 1;              %Currently the Bin Width is set to 0.5MuSeconds thick
hhD2.FaceColor = [0.8 0.8 0.8];   
hhD2.EdgeColor = 'k';             %'k'-Blue outline.
hold on

[hightD,CentrBn] = hist(MuD,BinRD);   %BinCnt,CntrD
    %    F(R) - F(L)  =  h*(R - L).
hold on
nhD = length(CentrBn);          %number of bars (Equals 20)
whD = CentrBn(2)-CentrBn(1);     %Muon Decay split into bits by a factor of 1.
                                    %^^^^^^^^Defined by BinRange^^^^(Equals 1)
thD = linspace(CentrBn(1)-(whD/2), CentrBn(20)+whD/2, nhD+1);

phD = fix(nhD/2);       %"Fix" rounds value to nearest whole value
lnP = plot(CentrBn([phD phD]),[0 hightD(phD)],'b:','LineWidth', 2, 'MarkerSize',2);
h = text(CentrBn(phD)-.2,hightD(phD)/2,'     h');
h1 = text(CentrBn(1)-.2,hightD(1)/2,'     hi');
h2 = text(CentrBn(2)-.2,hightD(2)/2,'     hi+1');
h3 = text(CentrBn(3)-.2,hightD(3)/2,'     hi+2');
%h(i) its height 
hold on

dt = diff(thD); %dt(i its width
Fvals = cumsum([0.04,hightD.*dt]); % |F| is complete cubic spline interpolant.
F = spline(thD, [0.04, Fvals, 19.96]); %Bounds From 0.04 to 19.96
DivF = fnder(F);  % computes its first derivative
fnplt(DivF, 'r', 2)   %Function Plot
hold on

title(' No. of Muon Decay event N(t) between 0.04 and 19.96 microsec ')
xlabel('Time in Micro-seconds (X 10^-6)')
ylabel('Number of Decays events')
grid on
legend('Muon Decay Events','h center','CubicSplineFit');
hold on
%--------------------------------------------------------------------------------------------------------------------------------------------------------------
MuBC = mean(BinCnt);                   %Average Bin Height
sigmaBC = std(BinCnt);                 %StandardD
maximumDBC = max(BinCnt);              %Max
minimumDBC = min(BinCnt);              %Min
%--------------------------------------------------------------------------------------------------------------------------------------------------------------
disp(d1)
BinCountData = table(MuBC,sigmaBC,minimumDBC,maximumDBC,'VariableNames',{'AverageDecay','StandardDeviation','MinimumDecay','MaximiumDecay'})
disp(d1)
%--------------------------------------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------------
%Variables for Normalized, Width0.5 Histogram with Default Bin Width. 
DeltaN5 = NoD/40;                        %Change in N = NumberOfCounts/BinCount
%SrtB5 = sqrt(BinCnt5);                    %SquareRoot of the Number of Decay for uncertainty
BinR5 = 0.04:0.5:19.96;                %Bin Range
BinR5 = BinR5(:); 
%--------------------------------------------------------------------------------------
%-----------------------------Normalized, Width0.5 ------------------------------------
%--------------------------------------------------------------------------------------
% Figure 5
figure('Name','Trial_Two, Banwidth of 0.5MicroSeconds')          
h5 = histogram(MuD);
h5.BinWidth = 0.5;              %Currently the Bin Width is set to 0.5MuSeconds thick
h5.FaceColor = [0.5 0.5 0.5];   
h5.EdgeColor = 'k';             %'k'-Blue outline.
hold on  
BinCnt5 = histc(MuD,BinR5);
        %[counts,centers] = histogram-bin counts of data set(MuD),
                    %with range BinRD
BinR5 = BinR5(:);               %(:) converts Bin Range to column
BinCnt5 = BinCnt5(:);             %(:) converts Bin Count into a column
CentrBn5 = (BinR5+0.5);             %Center of Each Bin

hold on
title(' No. of Muon Decay event N(t) between 0.12 and 20 microsec ')
xlabel('Time in Micro-seconds (X 10^-6)')
grid
hold on

Xq5 = CentrBn5(:,1); %Interval the DecayEventData set down collum one
Yq5 = BinCnt5(:,1); %Interval the Bin Range set down collum one
format long
LLambdaq5 = Xq5\Yq5;       %llambda is slope
yCalcq5 = (LLambdaq5)*Xq5;
q5 =plot(Xq5,Yq5,':b','LineWidth',2,'MarkerSize',3);
hold on
q6 =plot(Xq5,yCalcq5,'r','LineWidth',2,'MarkerSize',3);
hold on
Yq5LN =(1/LLambdaq5);
qln = plot(Yq5LN, DecayDomain,'g','LineWidth',2,'MarkerSize',3);

%-----------------------------------------------------------------------------------
%-----------------------------------------BackGround---------------------------------
%-----------------------------------------------------------------------------------
BGCntr5 = [12.29;12.79;13.29;13.79;14.29;14.79;15.29;15.79;16.29;16.79;17.29;17.79;18.29;18.79;19.29;19.79];
BGRng5 = [20;20;25;19;22;20;16;26;14;16;16;20;20;33;19;0];
BinR6 = 0.012:0.5:19.96;        %start, interval (0.2microseconds), end
[BinC6] = histc(MuD,BinR6); 

ylabel('Number of Decays events')
q7 = bar(BGCntr5,BGRng5);
q7.FaceColor = [0 0 1];   
q7.EdgeColor = [0 0 1];

legend([h5 q5 q6 q7],'Histogram of MuonDecay Events',' Center of Bin vs. BinCound', 'BinRange vs. BinCount', 'BackGround Range')
hold on

title('Muon Decay Life vs. Number of Decays (width 0.5)')
xlabel('Time in Micro-seconds (X 10^-6)')
ylabel('Number of Decay Events')

% Figure 6
figure('Name','Non-Linear Least Square Fit')
LeastSqrtFitTot5 = lsqnonneg(CentrBn5,BinCnt5);   %(4.0319)
LeastSqrtFitBG5 = lsqnonneg(BGCntr5,BGRng5);             %(2.4711)

FittMLS = times(MuD,LeastSqrtFitTot5);
FittMLSBG = times(MuD,LeastSqrtFitBG5);

q4 = histogram(FittMLS);
q4.FaceColor = [0 1 0.7];   
q4.EdgeColor = [0 0 0];
hold on

qq5 = histogram(FittMLSBG);
qq5.FaceColor = [0.5 0.9 0.2 ];   
qq5.EdgeColor = [0 0 0];
hold on


[hight5,CentrBn5] = hist(MuD,BinR5);   %BinCnt,CntrD
    %    F(R) - F(L)  =  h*(R - L).
hold on
nh5 = length(CentrBn5);          %number of bars (Equals 20)
wh5 = CentrBn5(2)-CentrBn5(1);     %Muon Decay split into bits by a factor of 1.
                                    %^^^^^^^^Defined by BinRange^^^^(Equals 1)
th5 = linspace(CentrBn5(1)-(wh5/2), CentrBn5(40)+wh5/2, nh5+1);

ph5 = fix(nh5/2);       %"Fix" rounds value to nearest whole value
ln5 = plot(CentrBn5([ph5 ph5]),[0 hight5(ph5)],'r:','LineWidth', 2, 'MarkerSize',2);
h5p = text(CentrBn5(ph5)-.2,hight5(ph5)/2,'     h');
h51 = text(CentrBn5(1)-.2,hight5(1)/2,'     hi');
h52 = text(CentrBn5(2)-.2,hight5(2)/2,'     hi+1');
h53 = text(CentrBn5(3)-.2,hight5(3)/2,'     hi+2');
%h(i) its height 
hold on
dt5 = diff(th5); %dt(i its width
Fvals5 = cumsum([0.04,hight5.*dt5]); % |F| is complete cubic spline interpolant.
F5 = spline(th5, [0.04, Fvals5, 19.96]); %Bounds From 0.04 to 19.96
DivF5 = fnder(F5);  % computes its first derivative
hold on
fnplt(DivF5, 'r', 2)   %Function Plot
legend('NonLinear Square Fit','Subtracted BackGround','Fit');
title('Non-Linear Least Square Fit')
xlabel('Time in Micro-seconds (X 10^-6)')
ylabel('Number of Decay Events')
%--------------------------------------------------------------------------------------------------------------------------------------------------------------
MuBC5 = mean(BinCnt5);                   %Average Bin Height
sigmaBC5 = std(BinCnt5);                 %StandardD
maximumDBC5 = max(BinCnt5);              %Max
minimumDBC5 = min(BinCnt5);              %Min
%--------------------------------------------------------------------------------------------------------------------------------------------------------------
disp(d1)
BinCount05wdth = table(MuBC5,sigmaBC5,minimumDBC5,maximumDBC5,'VariableNames',{'AveragDcay','StndardDviation','MinimumDcay','MaximiumDcay'})
disp(d1)
%--------------------------------------------------------------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------------
%-------------------------------------Normalized, Width 0.2-------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------------------------------------------------------

% Variables
DeltaN2 = NoD/100;                     %Change in N = NumberOfCounts/BinCount (28.6200)
SrtB2 = sqrt(BinCnt2);                 %SquareRoot of the Number of Decay for uncertainty
BGh22 = 129.332;                      %Background = (Sum DecayEvent from 12 to 20MuSec)/(Number of Bins 12 to 20 MuS)
beta02 = 0.2;                         % Intervaled at 0.2microseconds
BinR2 = 0.012:0.2:19.96;             %start, interval (0.2microseconds), end
BinR2 = BinR2(:);                    % (:) creates a vertical array
%-----------------------------------------------------------------------------------------------------------------------
%Figure 7)
figure
X2 = CentrBn2(:,1);  %BinCnt2
Y2 = BinCnt2(:,1); %CentrBn2
format long
b1 = X2\Y2
y2Calc1 = b1*X2;
scatter(X2,Y2)
hold on

scatter(Y2,y2Calc1,'r')
grid on
hold on

XX2 = [ones(length(X2),1) X2];
slp2 = XX2\Y2;
yCalc2 = XX2*slp2;
plot(X2,yCalc2,'ok')
legend('DecayEvents','Linear','Slope & Intercept','Location','best');

RmnS02 = 1 - sum((Y2 - y2Calc1).^2)/sum((Y2 - mean(Y2)).^2)
RmnS021 = 1 - sum((Y2 - yCalc2).^2)/sum((Y2 - mean(Y2)).^2)

%Figure 8)
figure('Name','Trial_ Two, Banwidth of 0.2MicroSeconds')              % Figure 2
ht2 = histogram(MuD);  %#ok<*NOPTS>
ht2.BinWidth = 0.2;              %Currently the Bin Width is set to 0.2MuSeconds thick
ht2.FaceColor = [0 0.7 0.5];   
ht2.EdgeColor = [0 0.7 0.5];             
hold on

BinCnt2 = histc(MuD,BinR2);
%[counts,centers] = hist(___)

title(' No. of Muon Decay event N(t) vs. time width 0.2MicroSecond Width ')
xlabel('Time in Micro-seconds (X 10^-6)'),
ylabel('Number of Decays events (BinWidth 0.2MicroSecond)')
grid on
hold on
BinCnt2 = BinCnt2(:);                    %(:) converts Bin Count into a column

p020 = plot(CentrBn2,BinCnt2, 'b');
hold on

p021 = plot(BinR2,BinCnt2,'-.r');
hold on

%-----------------------------------------------------------------------------------
%-----------------------------BackGround, Width 0.2---------------------------------
%-----------------------------------------------------------------------------------
p022 = bar(BGCenr2,BGR2);   %Plot of Back Ground
p022.FaceColor = [0 0 1];   
p022.EdgeColor = [0 0 1];
hold on
legend([ht2 p020 p021 p022 ],'Histogram of MuonDecay Events','Center vs, BinCount','BinRange vs. BinCount', 'BackGound Range')
hold off

%Figure 9)
figure('Name','Non-Linear Least Square Fit')
LeastSqrtFitTot2 = lsqnonneg(BinR2,BinCnt2);   %    (4.0319)
LeastSqrtFitBG2 = lsqnonneg(BGCenr2,BGR2);     %    (2.4711)

FittMLS2 = times(DecayDomain,LeastSqrtFitTot2);
FittMLSBG2 = times(DecayDomain,LeastSqrtFitBG2);

p23 = histogram(FittMLS2);
p23.FaceColor = [0.5 0.9 0.2];   
p23.EdgeColor = [0 0 0];
hold on
p24 = histogram(FittMLSBG2);
p24.FaceColor = [0 0.1 0.7];   
p24.EdgeColor = [0 0 0];
hold on
legend([p23 p24],'NonLinearSquareFit','BackGroundRemoved')
hold on
%--------------------------------------------------------------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------Removing BackGround From 0.2width-------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------------------------------------------------------
[X2pos,resnorm2,residual2] = lsqnonneg(CentrBn2,BinCnt2); %Obtain the solution and residual information.
  NormTest2 = norm(residual2)^2;   %^LineSquareNonNegativeFunction
        %^Test Residual Value for Validity        
%Nonnegative Least Squares with Nondefault Options
options2 = optimset('Display','final');

%Figure 10)
figure;
h02w = histogram(MuD);
h02w.BinWidth = 0.2; 
hold on

h2wr = bar(CentrBn2, residual2);     
h2wr.FaceColor = [0 1 1];   
h2wr.EdgeColor = [0 0 0];  

title(' (Chi^2) : residual Muon Decay event N(t) between 0.12 and 20 microsec ')
xlabel('Time in Micro-seconds (X 10^-6)')
ylabel('Number of Decays events (BinWidth 0.2MicroSecond)')
grid on
legend([h02w h2wr], 'Length of Muon Decay Events vs. Binwidth 0.2', 'Chi^2 residual data vs. time')
hold on

%--------------------------------------------------------------------------------------------------------------------------------------------------------------

%Figure 11)
figure('Name','MuonDecayEvents Fitted to Line'); 

h2wFL = histogram(MuD); %#ok<*NOPTS> 
h2wFL.BinWidth = 0.2;              %Currently the Bin Width is set to 0.2MuSeconds thick
h2wFL.FaceColor = [0.5 0.5 0.5];   
h2wFL.EdgeColor = 'k';             %'k'-black outline.
hold on

[hght,CntrBn002] = hist(DecayDomain,BinR2);   %BinCnt,CntrD
    %    F(R) - F(L)  =  h*(R - L).
hold on
nh2 = length(CntrBn002);          %number of bars (Equals 20)
wh2 = CntrBn002(2)-CntrBn002(1);     %Muon Decay split into bits by a factor of 1.
                                    %^^^^^^^^Defined by BinRange^^^^(Equals 1)
th2 = linspace(CntrBn002(1)-(wh2/2), CntrBn002(100)+wh2/2, nh2+1);

ph2 = fix(nh2/2);       %"Fix" rounds value to nearest whole value
lnP2 = plot(CntrBn002([ph2 ph2]),[0 hght(ph2)],'b:','LineWidth', 2, 'MarkerSize',2);
hlt = text(CntrBn002(ph2)-.2,hght(ph2)/2,'     h');
hlt = text(CntrBn002(1)-.2,hght(1)/2,'     hi');
hlt = text(CntrBn002(2)-.2,hght(2)/2,'     hi+1');
hlt = text(CntrBn002(3)-.2,hght(3)/2,'     hi+2');
%h(i) its height
hold on
dtt = diff(th2); %dt(i its width
Fvals2 = cumsum([0.04,hght.*dtt]); % |F| is complete cubic spline interpolant.
F2 = spline(th2, [0.04, Fvals2, 19.96]); %Bounds From 0.04 to 19.96
DivF2 = fnder(F2);  % computes its first derivative
hold on
fnplt(DivF2, 'r', 2)   %Function Plot

title('Muon Life time vs. No. of Muon Decay event N(t)')
xlabel('Time in Micro-seconds (X 10^-6)')
ylabel('Number of Decays events')
grid on
legend('Muon Decay Events','h center', 'SplineFit')
%--------------------------------------------------------------------------------------------------------------------------------------------------------------
MuBC2 = mean(BinCnt2);                   %Average Bin Height
sigmaBC2 = std(BinCnt2);                 %StandardD
maximumDBC2 = max(BinCnt2);              %Max
minimumDBC2 = min(BinCnt2);              %Min
%--------------------------------------------------------------------------------------------------------------------------------------------------------------
disp(d1)
BinCount2wdth = table(MuBC2,sigmaBC2,minimumDBC2,maximumDBC2,'VariableNames',{'AverageDecay','StandardDeviation','MinimumDecay','MaximiumDecay'})
disp(d1)
%--------------------------------------------------------------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------------------------------------------------------
disp(d1)

%--------------------------------------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------Data Fitting-------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------------------------------------------------------
%NonLinear Least Square Fit for Ln(N(t)) vs. t
avgYln = mean(yln)
NLSFlnnt = lsqnonneg(yln,MuD);              %(3.8468)
lnntUnc = yln\MuD;                          %(3.8468)
constrtNm = norm(yln*NLSFlnnt - MuD);        % ConstrainedNormalization (168.8596)
UncnstrtNorm = norm(yln*lnntUnc - MuD);      % UnConstrainedNormalization (168.8596)

%NonLinear Least Square Fit for Ln(N(t)) vs. t(T is organized low to high)
NLSFlnnt2 = lsqnonneg(yln,DecayDomain);              %(1.297269198379082)
lnntUnct2 = yln\DecayDomain;                          %(1.297269198379082)
constrtNm2 = norm(yln*NLSFlnnt2 - DecayDomain);   % ConstrainedNormalization (3.331691199826960e+02)
UncnstrtNorm2 = norm(yln*lnntUnct2 - DecayDomain); % UnConstrainedNormalization (3.331691199826960e+02)

dvdFrtyln = avgYln\40;
chsqdd = ((dvdFrtyln - Expected)^2)\Expected;
%--------------------------------------------------------------------------------------------------------------------------------------------------------------
% ln(P(llambda, t)) = ln(llambda -llambda*time)
dlnp = (1/llambda) - MuD;               %div[ln(P)]/divllambda = (1/llambda) - Time;
dlnp2 = (1/LLambda) - MuD;              %div[ln(P)]/div LLambda = (1/Slope (^Above^)) - Time; 

dlnppos = ((1/llambda) - MuD)+0.5;      %div[ln(P)]/divllambda = (1/llambda) - Time;
dlnppos2 = ((1/LLambda) - MuD)+0.5;     %div[ln(P)]/div LLambda = (1/LLambda (^Above^)) - Time;

dlnpneg = ((1/llambda) - MuD)-0.5;      %div[ln(P)]/divllambda = (1/llambda) - Time;
dlnpneg2 = ((1/LLambda) - MuD)-0.5;     %div[ln(P)]/div LLambda = (1/LLambda (^Above^)) - Time;

Avgdlnp = mean(dlnp);                   %Average of div[ln(P)]/divllambda
Avgdlnp2 = mean(dlnp2);                 %Average of div[ln(P)]/div LLambda

sigmaLl = std(dlnp);                    %StandardD (Should be about +/-0.5)
sigmaLl2 = std(dlnp2);                  %StandardD LLambda(about +/-0.5)
LlPos = Avgdlnp+sigmaLl;                % SD +0.5
LlPos2 = Avgdlnp2+sigmaLl2;             % SDLLambda +0.5
LlNeg = Avgdlnp-+sigmaLl;               % SD -0.5
LlNeg2 = Avgdlnp2-+sigmaLl2;            % SDLLambda -0.5

dvdFrtyln2 = (1\llambda)\40;
chsqdd2 = ((dvdFrtyln2 - Expected)^2)\Expected;
%--------------------------------------------------------------------------------------------------------------------------------------------------------------
%SlopeFit(NumberofDecays/Sum) = Equations with LLambda (2)
%SlopeFit2(Max Likely Hood) = Equations with LLAMBDA - BackGround
%--------------------------------------------------------------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------------------------------------------------------
SlopeFit = table(llambda,Avgdlnp,sigmaLl,LlPos,LlNeg,'VariableNames',{'ApproxLLAMBDA','ApproxAverage','ApproxStandardDev','StndrdDevPos','StndrdDevNeg'}) %#ok<NASGU>
disp(d1)
SlopeFit2 = table(LLambda, Avgdlnp2,sigmaLl2,LlPos2,LlNeg2,'VariableNames',{'SlpofLLAMBDAminsBckGrnd','SlpApproxAvge','SlpStandardDev','StndrdDevPos','StndrdDevNeg'}) %#ok<NASGU>
disp(d1)
FitEvntData = table(NLSFlnnt,constrtNm,NLSFlnnt2,constrtNm2,'VariableNames',{'LnMuD','Normlized','LnDomain','Normlzed'}) %#ok<NASGU>
disp(d1)
%-------------------------------------------------------------------
%Figure 12)
yyaxis left
figure('Name','ln(MuonDecay Events) vs. MuonDecay Events')
hDf = histogram(MuD);           
hDf.FaceColor = [0 0.5 0.5];
hDf.EdgeColor = 'k';
ylabel('Muon Decay Count')
hold on

%Linear Fit
yyaxis right
Xln = CentrBn2(:,1); %Interval the DecayEventData-BackGround set down collum one
Yln = MudMinsBG(:,1); %Interval the Bin Range set down collum one
format long
LLambda = Xln\Yln;       %llambda is slope
yCalc1 = times(Xln,-LLambda);   %Slope multiplied by x position 
plot(Xln,yCalc1,'y','LineWidth', 3);
hold on


legend('Histogram','BackGround'); %'LinearFit minus Background',
xlabel('Linear Event Distribution from 0.04MicroSeconds to 19.96MicroSeconds')
title('Linear Regression  of Muon Decays')
grid on
%--------------------------------------------------------------------------------------------------------------------------------------------------------------
%---------------------------------Chi Squared------------------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------------------------------------------------------
dvdFrtyln3 = (1\LLambda);
chsqdd3 = ((dvdFrtyln3 - Expected)^2)\Expected;



%Figure 13)
figure('Name','MuonDecay Lifespan vs. Events Chi with probability distribution')
ProbDisthD2 = pdf('Chisquare',MuD,2,sigma);
ProbDisthD22 = pdf('Chisquare',DecayDomain,2,sigma); 
ChiMu = times(ProbDisthD2,MuD);         %Average probability density was 0.212997742826911
ChiMu2 = times(ProbDisthD22,DecayDomain);

yyaxis left
hChi = histogram(MuD);            
hChi.FaceColor = [0.8 0.8 1];   
hChi.EdgeColor = 'k';
ylabel('Number of Decays events (BinWidth 1MicroSecond)')
hold on

yyaxis right
scatter(MuD,ProbDisthD2,'b');
ylabel('Decays event distribution (BinWidth 1MicroSecond)')
hold on

scatter(MuD,ChiMu,'g','LineWidth', 0.5);
hold on

plot(DecayDomain,ProbDisthD22,'r','LineWidth', 3);
hold on
plot(DecayDomain,ChiMu2,'y','LineWidth', 3);
hold on

title('Decay Event Probablity distribution')
xlabel('Time in Micro-seconds (X 10^-6)')
grid on
legend('DecayLifetime vs. number of Events','Positive Distribution Curve', 'Chi^2 Probability Density')

%Figure 14)
figure('Name','Cubic Spline Interpolant of DecayEvents')
hCUbe = histogram(MuD);          
hCUbe.FaceColor = [0.4 0.4 0.5];   
hCUbe.EdgeColor = 'k';             
hold on
XX = linspace(minimumD,maximumD,2862);
%CentrBn2,BinR2
%CentrBn,BinCounthD
plot(XX,csapi(CentrBn,BinCounthD,XX),'k:',CentrBn,BinCounthD,'r')
title('Cubic Spline Interpolantion of Muon Event Data')
xlabel('Time in Micro-seconds (X 10^-6)')
ylabel('Number of Decays events (BinWidth 1MicroSecond)')
grid on
legend('Muon Life time vs. Number of Decay Events.','Cubic Splin CenterBin vs.BinRange','BinRange 0.2wdth')

%Figure 15)
figure('Name','MuonDecay Lifespan vs. Events Chi with probability distribution')
error = (BinCnt2).^(1\2);
plot(CentrBn2, residual2,'r','LineWidth', 3);
SLopeOFBackGroundRemoved = residual2\CentrBn2;
OneOver = 1\SLopeOFBackGroundRemoved
hold on
plot(CentrBn2,NonLinearLeastSquareFit,'y','LineWidth', 3);
hold on
legend('show');
title('Decay Event Probablity distribution')
xlabel('Time in Micro-seconds (X 10^-6)')
grid on
legend('Positive Distribution Curve', 'Chi^2 Probability Density')
%--------------------------------------------------------------------------------------------------------------------------------------------------------------
ChiSquared = table(chsqdd,chsqdd2,chsqdd3,'VariableNames',{'chiLNMuonDecayData','chiApproxOneOverLLAMBDA','chiSqOneOverLlambdaMnusBckGrnd'}) %#ok<NASGU>
disp(d1)
%--------------------------------------------------------------------------------------------------------------------------------------------------------------
disp(d1)
disp(lab)
disp(win)
disp(d1)
%-----------------------------------------------------------------
diary ('MuonHistoT2v9')                               %Save File
%whos                                                 %FileSize
%-----------------------------------------------------------------