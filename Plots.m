%% This is a script to load data and create plots for EARS 202 final project
addpath("Matlab functions/")
%% Alignments

AI = readmatrix("Output/AlignedIndices.csv");
AI = AI(:,2);
aligned = readmatrix("Output/AlignedSignals.csv");


%% All cores
cores = {'DEN13A','DEN13B','DEN19A','DEN22A','DEN22B'};
Ncore = length(cores); % number of cores
ab 
%% Loads isotope and CFA data
for c = 1:Ncore
    core = cores{c};
    ISO.(core) = readtable(['Data/Denali/Isotopes/' core '_WithTime.csv']);
    if c == 2
        CFA.(core) = readtable(['Data/Denali/Dust/' core '_16_Dust.csv']);
    elseif c < 5
        CFA.(core) = readtable(['Data/Denali/Dust/' core '_22_Dust.csv']);
    end
end

%% Loads IC / ICPMS Data
CHEM.DEN13A = readtable("Data/Denali/DIC1_IC_ICPMS.csv");
CHEM.DEN13B = readtable("Data/Denali/DIC2_IC_ICPMS.csv");
CHEM.DEN19A = readtable("Data/Denali/Denali_2019_ICPMS_AllData.csv");
CHEM.DEN22A = readtable("Data/Denali/Denali_2022_ICPMS.csv");

%% Loads agescales
TS.DEN13A_2 = readtable("Data/Denali/DICfire_Composite_AgeScale_VolcAdjusted.csv");
TS.DEN13A = readtable("Data/Denali/DIC1_Timescale_20230620.csv");
TS.DEN13B = readtable("Data/Denali/DIC2_Timescale_20230620.csv");
TS.DEN19A = readtable("Data/Denali/DEN19A_AgeScale_20250223.csv");
TS.DEN22A = readtable("Data/Denali/DEN22A_AgeScale_20250223.csv");

% adds drill date to each agescale, approximated as 0.5 years past year start
TS.DEN13A_2 = [{2013.5,0};TS.DEN13A_2];
TS.DEN13B = [{2013.5,0};TS.DEN13B];
TS.DEN19A = [{0, 2019.5};TS.DEN19A];
TS.DEN22A = [{0, 2022.5};TS.DEN22A];

%% Reinterpolates most recent age scale to chemistry data
% DEN13A
ISO.DEN13A.TopTime = interp1(TS.DEN13A.DepthJan_1_m_,TS.DEN13A.Year_CE_,ISO.DEN13A.TopDepth);
ISO.DEN13A.MidTime = interp1(TS.DEN13A.DepthJan_1_m_,TS.DEN13A.Year_CE_,ISO.DEN13A.MidDepth);
ISO.DEN13A.BotTime = interp1(TS.DEN13A.DepthJan_1_m_,TS.DEN13A.Year_CE_,ISO.DEN13A.BotDepth);
CHEM.DEN13A.TopTime = interp1(TS.DEN13A.DepthJan_1_m_,TS.DEN13A.Year_CE_,CHEM.DEN13A.TopDepth);
CHEM.DEN13A.MidTime = interp1(TS.DEN13A.DepthJan_1_m_,TS.DEN13A.Year_CE_,CHEM.DEN13A.MidDepth);
CHEM.DEN13A.BotTime = interp1(TS.DEN13A.DepthJan_1_m_,TS.DEN13A.Year_CE_,CHEM.DEN13A.BotDepth);
CFA.DEN13A.Time = interp1(TS.DEN13A_2.Depth_m_,TS.DEN13A_2.Year,CFA.DEN13A.Depth_m_);


% DEN13B
ISO.DEN13B.TopTime = interp1(TS.DEN13B.DepthJan_1_m_,TS.DEN13B.Year_CE_,ISO.DEN13B.TopDepth);
ISO.DEN13B.MidTime = interp1(TS.DEN13B.DepthJan_1_m_,TS.DEN13B.Year_CE_,ISO.DEN13B.MidDepth);
ISO.DEN13B.BotTime = interp1(TS.DEN13B.DepthJan_1_m_,TS.DEN13B.Year_CE_,ISO.DEN13B.BotDepth);
CHEM.DEN13B.TopTime = interp1(TS.DEN13B.DepthJan_1_m_,TS.DEN13B.Year_CE_,CHEM.DEN13B.TopDepth);
CHEM.DEN13B.MidTime = interp1(TS.DEN13B.DepthJan_1_m_,TS.DEN13B.Year_CE_,CHEM.DEN13B.MidDepth);
CHEM.DEN13B.BotTime = interp1(TS.DEN13B.DepthJan_1_m_,TS.DEN13B.Year_CE_,CHEM.DEN13B.BotDepth);
CFA.DEN13B.Time = interp1(TS.DEN13B.DepthJan_1_m_,TS.DEN13B.Year_CE_,CFA.DEN13B.Depth_m_);

% DEN19A
ISO.DEN19A.TopTime = interp1(TS.DEN19A.Depth_m_,TS.DEN19A.Year,ISO.DEN19A.TopDepth);
ISO.DEN19A.MidTime = interp1(TS.DEN19A.Depth_m_,TS.DEN19A.Year,ISO.DEN19A.MidDepth);
ISO.DEN19A.BotTime = interp1(TS.DEN19A.Depth_m_,TS.DEN19A.Year,ISO.DEN19A.BotDepth);
CHEM.DEN19A.TopTime = interp1(TS.DEN19A.Depth_m_,TS.DEN19A.Year,CHEM.DEN19A.depth_top);
CHEM.DEN19A.MidTime = interp1(TS.DEN19A.Depth_m_,TS.DEN19A.Year,CHEM.DEN19A.depth_avg);
CHEM.DEN19A.BotTime = interp1(TS.DEN19A.Depth_m_,TS.DEN19A.Year,CHEM.DEN19A.depth_bottom);
CFA.DEN19A.Time = interp1(TS.DEN19A.Depth_m_,TS.DEN19A.Year,CFA.DEN19A.Depth_m_);

% DEN22A
ISO.DEN22A.TopTime = interp1(TS.DEN22A.Depth_m_,TS.DEN22A.Year,ISO.DEN22A.TopDepth);
ISO.DEN22A.MidTime = interp1(TS.DEN22A.Depth_m_,TS.DEN22A.Year,ISO.DEN22A.MidDepth);
ISO.DEN22A.BotTime = interp1(TS.DEN22A.Depth_m_,TS.DEN22A.Year,ISO.DEN22A.BotDepth);
CHEM.DEN22A.TopTime = interp1(TS.DEN22A.Depth_m_,TS.DEN22A.Year,CHEM.DEN22A.depth_top);
CHEM.DEN22A.MidTime = interp1(TS.DEN22A.Depth_m_,TS.DEN22A.Year,CHEM.DEN22A.depth_avg);
CHEM.DEN22A.BotTime = interp1(TS.DEN22A.Depth_m_,TS.DEN22A.Year,CHEM.DEN22A.depth_bottom);
CFA.DEN22A.Time = interp1(TS.DEN22A.Depth_m_,TS.DEN22A.Year,CFA.DEN22A.Depth_m_);

% DEN22B
ISO.DEN22B.TopTime = interp1(TS.DEN22A.Depth_m_,TS.DEN22A.Year,ISO.DEN22B.TopDepth);
ISO.DEN22B.MidTime = interp1(TS.DEN22A.Depth_m_,TS.DEN22A.Year,ISO.DEN22B.MidDepth);
ISO.DEN22B.BotTime = interp1(TS.DEN22A.Depth_m_,TS.DEN22A.Year,ISO.DEN22B.BotDepth);


%% Resampling isotope data to annual and 3-year averages
ISO.DEN13A.Deuterium(ISO.DEN13A.Deuterium == 0) = NaN;
ISO.DEN13A.Oxygen18(ISO.DEN13A.Deuterium == 0) = NaN;
ISO.DEN13B.Deuterium(ISO.DEN13B.Deuterium == 0) = NaN;
ISO.DEN13B.Oxygen18(ISO.DEN13B.Deuterium == 0) = NaN;
years1 = 1700:2012;
years3 = 1201:3:2011;
interpX = 1200:0.01:2013;
isoInterp.DEN13A = nan(length(interpX),2);
isoInterp.DEN13B = nan(length(interpX),2);
for i = 1:size(ISO.DEN13B,1)
    top = ISO.DEN13B.TopTime(i);
    bot = ISO.DEN13B.BotTime(i);
    [~,i1] = min(abs(top-interpX));
    [~,i2] = min(abs(bot-interpX));
    isoInterp.DEN13B(i2:i1,1) = ISO.DEN13B.Deuterium(i);
    isoInterp.DEN13B(i2:i1,2) = ISO.DEN13B.Oxygen18(i);
end
for i = 1:size(ISO.DEN13A,1)
    top = ISO.DEN13A.TopTime(i);
    bot = ISO.DEN13A.BotTime(i);
    [~,i1] = min(abs(top-interpX));
    [~,i2] = min(abs(bot-interpX));
    isoInterp.DEN13A(i2:i1,1) = ISO.DEN13A.Deuterium(i);
    isoInterp.DEN13A(i2:i1,2) = ISO.DEN13A.Oxygen18(i);
end
isoInterp.DEN13A(:,1) = naninterp(isoInterp.DEN13A(:,1));
isoInterp.DEN13A(:,2) = naninterp(isoInterp.DEN13A(:,2));
isoInterp.DEN13B(:,1) = naninterp(isoInterp.DEN13B(:,1));
isoInterp.DEN13B(:,2) = naninterp(isoInterp.DEN13B(:,2));

iso1.DEN13A = nan(length(years1),2);
iso1.DEN13B = nan(length(years1),2);
iso3.DEN13A = nan(length(years3),2);
iso3.DEN13B = nan(length(years3),2);
for y = 1:length(years1)
    yr = years1(y);
    i = floor(interpX) == yr;
    iso1.DEN13A(y,:) = mean(isoInterp.DEN13A(i,:),1);
    iso1.DEN13B(y,:) = mean(isoInterp.DEN13B(i,:),1);
end
for y = 1:length(years3)
    yr = years3(y);
    i = interpX >= yr - 1 & interpX < yr + 2;
    iso3.DEN13A(y,:) = mean(isoInterp.DEN13A(i,:),1);
    iso3.DEN13B(y,:) = mean(isoInterp.DEN13B(i,:),1);
end

iso3.DEN13A(years3 < 1865,:) = NaN;
iso1.DEN13A(years1 < 1865,:) = NaN;



%% load reanalyses
REAN = {'NCEP','NCEP2','NARR','CR20'};
%REAN = {'NCEP','NCEP2','NARR','JRA55','CR20'};

NREAN = length(REAN);

%% PWT
yrs = 1980:2010; % yrs for analysis
for r = 1:NREAN
    dtaID = REAN{r};
    DTA.(dtaID) = readtable(['Data/Reanalysis/' dtaID '.csv']);

    % data wrangling
    LVL.(dtaID) = DTA.(dtaID).Properties.VariableNames(6:end); % atmospheric level
    for l = 1:length(LVL.(dtaID))
        lvl = LVL.(dtaID){l};
        LVL.(dtaID){l} = lvl(6:end);
    end
    LVL.(dtaID) = double(string(LVL.(dtaID)));
    TEMP.(dtaID) = table2array(DTA.(dtaID)(:,6:end)) - 273.15;
    PRECIP.(dtaID) = table2array(DTA.(dtaID)(:,5));
    for i = 1:size(DTA.(dtaID),1)
        y = DTA.(dtaID).Year(i);
        m = DTA.(dtaID).Month(i);
        d = DTA.(dtaID).Day(i);
        dm = findDayIndex(y,m,d,y);
        ndays = isLeap(y) + 365;
        dy = y + dm/ndays;
        DECYR.(dtaID)(i) = dy;
    end

    % calculates annual T and PWT
    startYr = DTA.(dtaID).Year(1);
    endYr = DTA.(dtaID).Year(end);
    YEARS.(dtaID) = (startYr:endYr)';
    n = endYr - startYr + 1;
    Tann.(dtaID) = nan(n,size(TEMP.(dtaID),2));
    PWTann.(dtaID) = nan(n,size(TEMP.(dtaID),2));
    for y = startYr:endYr
        i = DTA.(dtaID).Year == y;
        t = TEMP.(dtaID)(i,:);
        p = PRECIP.(dtaID)(i);
        pwt = sum(t .* p) / sum(p);
        Tann.(dtaID)(y-startYr+1,:) = mean(t,1);
        PWTann.(dtaID)(y-startYr+1,:) = mean(pwt,1);
    end
    yi = ismember(YEARS.(dtaID),yrs);
    Tann.(dtaID) = Tann.(dtaID) - mean(Tann.(dtaID)(yi,:),1);
    PWTann.(dtaID) = PWTann.(dtaID) - mean(PWTann.(dtaID)(yi,:),1);

    % correlates annual T and PWT to ice core data
    PWT_R.(dtaID) = nan(size(TEMP.(dtaID),2),2);
    T_R.(dtaID) = nan(size(TEMP.(dtaID),2));
    PWT_R.(dtaID)(:,1) = corr(PWTann.(dtaID)(yi,:),iso1.DEN13B(ismember(years1,yrs),1));
    T_R.(dtaID)(:,1) = corr(Tann.(dtaID)(yi,:),iso1.DEN13B(ismember(years1,yrs),1));
    PWT_R.(dtaID)(:,2) = corr(PWTann.(dtaID)(yi,:),iso1.DEN13B(ismember(years1,yrs),2));
    T_R.(dtaID)(:,2) = corr(Tann.(dtaID)(yi,:),iso1.DEN13B(ismember(years1,yrs),2));

    % calculates biweekly T and PWT
    PWT_BIWEEK.(dtaID) = nan(n*26,size(TEMP.(dtaID),2));
    PRECIP_BIWEEK.(dtaID) = nan(n*26,1);
    T_BIWEEK.(dtaID) = nan(n*26,size(TEMP.(dtaID),2));
    for yr = 1:n
        yI = find(DTA.(dtaID).Year == YEARS.(dtaID)(yr),1);
        for w = 1:26
            i = yI(1)-1 + ((w-1)*14+1:(w*14));
            PRECIP_BIWEEK.(dtaID)((yr-1)*26+w) = sum(PRECIP.(dtaID)(i));
            PWT_BIWEEK.(dtaID)((yr-1)*26+w,:) = sum(PRECIP.(dtaID)(i) .* TEMP.(dtaID)(i,:)) / sum(PRECIP.(dtaID)(i));
            T_BIWEEK.(dtaID)((yr-1)*26+w,:) = mean(TEMP.(dtaID)(i,:));
        end
    end
    dates.(dtaID) = nan(n*26,3);
    for i = 1:n
        dates.(dtaID)((i-1)*26+1:i*26,1) = YEARS.(dtaID)(i);
        dates.(dtaID)((i-1)*26+1:i*26,2) = 1:26;
        dates.(dtaID)((i-1)*26+1:i*26,3) = YEARS.(dtaID)(i) + ((1:26)-1)/26;    
    end
end


%% Plotting parameters
font = 'Calibri';
labFS = 16;
axFS = 12;
titFS = 18;
legFS = 14;
cols = colororder;
pos1 = [584 528 700 400];
posv = [20 420 800 500];
res = 400;

%% 3-year water isotopes, 1200-2013, both cores
close all; figure('Position',pos1); hold on
plot(years3,iso3.DEN13A(:,1),'Color',cols(1,:),'LineWidth',2)
plot(years3,iso3.DEN13B(:,1),'Color',cols(2,:),'LineWidth',2)
set(gca,'FontSize',axFS,'FontName',font)
ylabel("δD (‰)",'FontSize',labFS)
xlabel("Year",'FontSize',labFS)
ylim([-240 -180])
xlim([1200 2012])
legend(cores,'location','northwest','FontSize',legFS,'NumColumns',2)
exportgraphics(gcf,'Figures/Isotopes_3Year.png','Resolution',res)

%% annual water isotopes, 1700-2012, both cores
close all; figure('Position',pos1); hold on
plot(years1,iso1.DEN13A(:,1),'Color',cols(1,:),'LineWidth',2)
plot(years1,iso1.DEN13B(:,1),'Color',cols(2,:),'LineWidth',2)
set(gca,'FontSize',axFS,'FontName',font)
ylabel("δD (‰)",'FontSize',labFS)
xlabel("Year",'FontSize',labFS)
ylim([-240 -180])
xlim([1700 2012])
legend(cores,'location','northwest','FontSize',legFS,'NumColumns',2)
exportgraphics(gcf,'Figures/Isotopes_Annual.png','Resolution',res)

%% annual water isotopes, 1980-2010, both cores
close all; figure('Position',pos1); hold on
plot(years1,iso1.DEN13A(:,1),'Color',cols(1,:),'LineWidth',2)
plot(years1,iso1.DEN13B(:,1),'Color',cols(2,:),'LineWidth',2)
set(gca,'FontSize',axFS,'FontName',font)
ylabel("δD (‰)",'FontSize',labFS)
xlabel("Year",'FontSize',labFS)
ylim([-240 -180])
xlim([1980 2010])
legend(cores,'location','northwest','FontSize',legFS,'NumColumns',2)
exportgraphics(gcf,'Figures/Isotopes_1980_2010.png','Resolution',res)

%% annual water isotopes, 1980-2010, ONLY CORE B
close all; figure('Position',pos1); hold on
plot(years1,iso1.DEN13B(:,1),'Color','k','LineWidth',2)
set(gca,'FontSize',axFS,'FontName',font)
ylabel("δD (‰)",'FontSize',labFS)
xlabel("Year",'FontSize',labFS)
ylim([-240 -180])
xlim([1980 2010])
legend([cores(2)],'location','northoutside','FontSize',legFS,'NumColumns',6)
exportgraphics(gcf,'Figures/Isotopes_1980_2010_CoreB.png','Resolution',res)

%% Reanalysis 1000 hPa T w/ annual water isotopes, 1980-2010
close all; figure('Position',pos1); hold on
plot(years1,iso1.DEN13B(:,1),'Color','k','LineWidth',2)
set(gca,'FontSize',axFS,'FontName',font)
ylabel("δD (‰)",'FontSize',labFS)
xlabel("Year",'FontSize',labFS)
ylim([-240 -180])
xlim([1980 2010])
yyaxis right
for r = 1:NREAN
    d = REAN{r};
    [~,lvli] = min(abs(LVL.(dtaID) - 1000));
    plot(YEARS.(d),Tann.(d)(:,lvli),'LineStyle','-','LineWidth',1.5,'Marker','none','Color',cols(r,:))
end
set(gca,'YColor','k')
ylabel("T_{anom} (°C)",'FontSize',labFS)
legend([cores(2) REAN],'location','northoutside','FontSize',legFS,'NumColumns',6)
exportgraphics(gcf,'Figures/Isotopes_Surface.png','Resolution',res)

%% Reanalysis 350 hPa T w/ annual water isotopes, 1980-2010
close all; figure('Position',pos1); hold on
plot(years1,iso1.DEN13B(:,1),'Color','k','LineWidth',2)
set(gca,'FontSize',axFS,'FontName',font)
ylabel("δD (‰)",'FontSize',labFS)
xlabel("Year",'FontSize',labFS)
ylim([-240 -180])
xlim([1980 2010])
yyaxis right
for r = 1:NREAN
    d = REAN{r};
    [~,lvli] = min(abs(LVL.(dtaID) - 350));
    plot(YEARS.(d),Tann.(d)(:,lvli),'LineStyle','-','LineWidth',2,'Marker','none','Color',cols(r,:))
end
set(gca,'YColor','k')
ylabel("T_{anom} (°C)",'FontSize',labFS)
legend([cores(2) REAN],'location','northoutside','FontSize',legFS,'NumColumns',6)
exportgraphics(gcf,'Figures/Isotopes_350hPa_T.png','Resolution',res)

%% Reanalysis 500 hPa PWT w/ annual water isotopes, 1980-2010
close all; figure('Position',pos1); hold on
plot(years1,iso1.DEN13B(:,1),'Color','k','LineWidth',2)
set(gca,'FontSize',axFS,'FontName',font)
ylabel("δD (‰)",'FontSize',labFS)
xlabel("Year",'FontSize',labFS)
ylim([-240 -180])
xlim([1980 2010])
yyaxis right
for r = 1:NREAN
    d = REAN{r};
    [~,lvli] = min(abs(LVL.(dtaID) - 500));
    plot(YEARS.(d),PWTann.(d)(:,lvli),'LineStyle','-','LineWidth',2,'Marker','none','Color',cols(r,:))
end
set(gca,'YColor','k')
ylabel("PWT_{anom} (°C)",'FontSize',labFS)
legend([cores(2) REAN],'location','northoutside','FontSize',legFS,'NumColumns',6)
exportgraphics(gcf,'Figures/Isotopes_500hPa_PWT.png','Resolution',res)

%% Reanalysis 500 hPa PWT (NARR emphasized) w/ annual water isotopes, 1980-2010
close all; figure('Position',pos1); hold on
plot(years1,iso1.DEN13B(:,1),'Color','k','LineWidth',3)
set(gca,'FontSize',axFS,'FontName',font)
ylabel("δD (‰)",'FontSize',labFS)
xlabel("Year",'FontSize',labFS)
ylim([-240 -180])
xlim([1980 2010])
yyaxis right
for r = 1:NREAN
    d = REAN{r};
    [~,lvli] = min(abs(LVL.(dtaID) - 500));
    if string(d) == "NARR"
        plot(YEARS.(d),PWTann.(d)(:,lvli),'LineStyle','-','LineWidth',3,'Marker','none','Color',cols(r,:))
    else
        plot(YEARS.(d),PWTann.(d)(:,lvli),'LineStyle','-','LineWidth',1,'Marker','none','Color',0.9*ones(3,1))
    end
end
ylabel("PWT_{anom} (°C)",'FontSize',labFS)
legend([cores(2) REAN],'location','northoutside','FontSize',legFS,'NumColumns',6)
exportgraphics(gcf,'Figures/Isotopes_500hPa_PWT_NARR.png','Resolution',res)


%% T and PWT / core correlation by pressure levels
figure('Position',posv)
for r = 1:NREAN
    dtaID = REAN{r};
    x = PWT_R.(dtaID)(:,1);
    subplot(1,2,1)
    plot(PWT_R.(dtaID)(:,1),LVL.(dtaID),'LineWidth',2)
    title("PWT Relationship")
    hold on
    subplot(1,2,2)
    plot(T_R.(dtaID)(:,1),LVL.(dtaID),'LineWidth',2)
    hold on
end
for i=1:2
    subplot(1,2,i)
    set(gca,'YDir','reverse','FontSize',axFS,'FontName',font)
    xlabel("r","FontSize",labFS)
    legend(REAN,'FontSize',legFS)
    grid on
    xlim([0 0.8])
    ylim([50 1000])
end
subplot(1,2,1)
title("PWT / isotope correlation",'FontSize',titFS)
ylabel("Level (hPa)",'FontSize',labFS)
subplot(1,2,2)
title("Temperature / isotope correlation",'FontSize',labFS)
exportgraphics(gcf,'Figures/Correlations.png','Resolution',res)


%% Raw water isotopes, all cores
close all; figure('Position',pos1); hold on
for c = 1:length(cores)
    core = cores{c};
    plot(ISO.(core).MidTime,ISO.(core).Deuterium,'LineWidth',2)
end
set(gca,'FontSize',axFS,'FontName',font)
ylabel("δD (‰)",'FontSize',labFS)
xlabel("Year",'FontSize',labFS)
xlim([2005 2020])
grid on
legend(cores,'location','northoutside','FontSize',legFS,'NumColumns',5)
exportgraphics(gcf,'Figures/Isotopes_Seasonality.png','Resolution',res)


%% PWT versus T, 1980-2010
close all; figure('Position',pos1);
s1=subplot(2,1,1);hold on
for r = 1:NREAN
    d = REAN{r};
    plot(YEARS.(d),Tann.(d)(:,1),'LineStyle','-','LineWidth',2,'Marker','none','Color',cols(r,:))
end
xlim([1980 2010])
set(gca,'FontSize',axFS,'FontName',font)
ylabel("T_{anom} (°C)",'FontSize',labFS)
legend(REAN,'location','northoutside','FontSize',legFS,'NumColumns',6)

s2=subplot(2,1,2);hold on
for r = 1:NREAN
    d = REAN{r};
    plot(YEARS.(d),PWTann.(d)(:,1),'LineStyle','-','LineWidth',2,'Marker','none','Color',cols(r,:))
end
xlim([1980 2010])
set(gca,'FontSize',axFS,'FontName',font)
xlabel("Year",'FontSize',labFS)
ylabel("PWT_{anom} (°C)",'FontSize',labFS)

s1.Position(2)=0.54;
s1.Position(4)=0.35;
s2.Position(4)=0.35;

exportgraphics(gcf,'Figures/PWT.png','Resolution',res)


%% NARR annual cycle T (also with precip)
close all; figure('Position',pos1);
plot(DECYR.NARR,TEMP.NARR(:,1),'LineStyle','-','LineWidth',0.5,'Color','k')
xlim([2000 2010])
set(gca,'FontSize',axFS,'FontName',font)
ylabel("T (°C)",'FontSize',labFS)
ylim([-10 25])
xlabel("Year",'FontSize',labFS)
exportgraphics(gcf,'Figures/NARR_1980_2010.png','Resolution',res)

xlim([2000 2010])
yyaxis right
plot(DECYR.NARR,PRECIP.NARR,'LineStyle','-','LineWidth',0.5,'Color','b')
ylabel("Precip (kg m^{–2})",'FontSize',labFS)
set(gca,'YColor','b')
exportgraphics(gcf,'Figures/NARR_1980_2010_Precip.png','Resolution',res)


%% NARR annual cycles - biweekly
close all; figure('Position',pos1);
[~,lvli] = min(abs(LVL.NARR - 500));
plot(dates.NARR(:,3),T_BIWEEK.NARR(:,lvli),'LineStyle','-','LineWidth',1,'Color','k')
xlim([2000 2010])
set(gca,'FontSize',axFS,'FontName',font)
ylabel("T (°C)",'FontSize',labFS)
ylim([-45 -10])
xlabel("Year",'FontSize',labFS)
exportgraphics(gcf,'Figures/NARR_1980_2010_BIWEEK.png','Resolution',res)

yyaxis right
plot(dates.NARR(:,3),PRECIP_BIWEEK.NARR,'LineStyle','-','LineWidth',1,'Color','b')
ylabel("Precip (kg m^{–2})",'FontSize',labFS)
set(gca,'YColor','b')
ylim([0 350])
exportgraphics(gcf,'Figures/NARR_1980_2010_Precip_BIWEEK.png','Resolution',res)

close all; figure('Position',pos1);
[~,lvli] = min(abs(LVL.NARR - 500));
plot(dates.NARR(:,3),PWT_BIWEEK.NARR(:,lvli),'LineStyle','-','LineWidth',1,'Color','k')
xlim([2000 2010])
set(gca,'FontSize',axFS,'FontName',font)
ylabel("PWT (°C)",'FontSize',labFS)
ylim([-45 -10])
xlabel("Year",'FontSize',labFS)
exportgraphics(gcf,'Figures/NARR_1980_2010_BIWEEK_PWT.png','Resolution',res)

yyaxis right
[x,y] = getStepped(ISO.DEN13B.BotTime,ISO.DEN13B.TopTime,ISO.DEN13B.Deuterium);
plot(x,y,'LineStyle','-','LineWidth',1,'Color','r')
ylabel("δD (‰)",'FontSize',labFS)
set(gca,'YColor','r')
exportgraphics(gcf,'Figures/NARR_1980_2010_BIWEEK_PWT_ISO.png','Resolution',res)



%% NARR and water isotopes - biweekly normalized
close all; figure('Position',pos1); hold on
plot(aligned(:,1),aligned(:,2),'LineStyle','-','LineWidth',1,'Color','k')
xlim([2000 2010])
set(gca,'FontSize',axFS,'FontName',font)
ylabel("Normalized values (unitless)",'FontSize',labFS)
xlabel("Year",'FontSize',labFS)
plot(aligned(:,1),aligned(:,3),'LineStyle','-','LineWidth',1,'Color','r')
legend(["PWT","δD"])
box on
exportgraphics(gcf,'Figures/NARR_1980_2010_BIWEEK_PWT_ISO_NORM.png','Resolution',res)


%% NARR and water isotopes - biweekly ALIGNED
close all; figure('Position',pos1); hold on
plot(aligned(:,1),aligned(:,4),'LineStyle','-','LineWidth',1,'Color','k')
xlim([2000 2010])
set(gca,'FontSize',axFS,'FontName',font)
ylabel("Normalized values (unitless)",'FontSize',labFS)
xlabel("Year",'FontSize',labFS)
plot(aligned(:,1),aligned(:,5),'LineStyle','-','LineWidth',1,'Color','r')
legend(["PWT","δD Aligned"])
box on
exportgraphics(gcf,'Figures/NARR_1980_2010_BIWEEK_ALIGNMENT.png','Resolution',res)

%% NARR and water isotopes - biweekly UNALIGNED and ALIGNED
close all; figure('Position',[584 528 600 500]); 
subplot(2,1,1); hold on
plot(aligned(:,1),aligned(:,2),'LineStyle','-','LineWidth',1,'Color','k')
xlim([2000 2010])
set(gca,'FontSize',axFS,'FontName',font)
ylabel("Norm. (unitless)",'FontSize',labFS)
xlabel("Year",'FontSize',labFS)
plot(aligned(:,1),aligned(:,3),'LineStyle','-','LineWidth',1,'Color','#bd4f00')
legend(["PWT","δD"],'NumColumns',2,'location','southwest')
box on
text(2000.3,2.5,"a",'FontName',font,'FontSize',titFS,'FontWeight','bold')

subplot(2,1,2); hold on
plot(aligned(:,1),aligned(:,4),'LineStyle','-','LineWidth',1,'Color','k')
xlim([2000 2010])
set(gca,'FontSize',axFS,'FontName',font)
ylabel("Norm. (unitless)",'FontSize',labFS)
xlabel("Year",'FontSize',labFS)
plot(aligned(:,1),aligned(:,5),'LineStyle','-','LineWidth',1,'Color','#2447b3')
legend(["PWT","δD Aligned"],'NumColumns',2,'location','southwest')
box on
text(2000.3,2.5,"b",'FontName',font,'FontSize',titFS,'FontWeight','bold')
exportgraphics(gcf,'Figures/NARR_1980_2010_BIWEEK_ALIGNMENT_SUBPLOT.png','Resolution',res)


%% Water isotopes - biweekly UNALIGNED and ALIGNED
close all; figure('Position',pos1); hold on
plot(aligned(:,1),aligned(:,3),'LineStyle','-','LineWidth',1,'Color','k')
xlim([2000 2010])
set(gca,'FontSize',axFS,'FontName',font)
ylabel("Normalized values (unitless)",'FontSize',labFS)
xlabel("Year",'FontSize',labFS)
plot(aligned(:,1),aligned(:,5),'LineStyle','-','LineWidth',1,'Color','r')
legend(["δD","δD Aligned"])
box on
exportgraphics(gcf,'Figures/NARR_1980_2010_BIWEEK_ALIGNMENT_ISO.png','Resolution',res)
xlim([2002 2005])
exportgraphics(gcf,'Figures/NARR_1980_2010_BIWEEK_ALIGNMENT_ISO_ZOOM.png','Resolution',res)


%% Chemistry aligned figures
YRS = 1980:2010;
NYRS=length(YRS);
% CFA
CFA_BIWEEK = nan(NYRS*26,1);
t = CFA.DEN13A.Year;
d = CFA.DEN13A.DustConcentation_ppb_;
for i = 1:NYRS
    y = YRS(i);
    decY = y+1/52:1/26:y+51/52;
    for j = 1:26
        startT = y+(j-1)/26;
        endT = y + j/26;
        CFA_BIWEEK((i-1)*26+j) = mean(d(t >= startT & t < endT),'omitnan');
    end
end

% Mg and MSA
Mg_BIWEEK = nan(NYRS*26,1);
MSA_BIWEEK = nan(NYRS*26,1);
t = CHEM.DEN13B.MidTime;
d = CHEM.DEN13B.ICMg2__ppb_;
e = CHEM.DEN13B.ICMSA_ppb_;
i = ~isnan(t);
t = t(i);
d = d(i);
e = e(i);
[~,i] = unique(t);
t = t(i);
d = d(i);
e = e(i);
for i = 1:NYRS
    y = YRS(i);
    decY = y+1/52:1/26:y+51/52;
    Mg_BIWEEK((i-1)*26+1:i*26) = interp1(t,d,decY);
    MSA_BIWEEK((i-1)*26+1:i*26) = interp1(t,e,decY);
end

% ISO
ISO_BIWEEK = nan(NYRS*26,1);
t = ISO.DEN13B.MidTime;
d = ISO.DEN13B.Deuterium;
i = ~isnan(t);
t = t(i);
d = d(i);
[~,i] = unique(t);
t = t(i);
d = d(i);
for i = 1:NYRS
    y = YRS(i);
    decY = y+1/52:1/26:y+51/52;
    ISO_BIWEEK((i-1)*26+1:i*26) = interp1(t,d,decY);
end


t = aligned(:,1);

XL = [2000 2010];
close all; figure('Position',[584 528 800 600])
subplot(4,1,1)
plot(t,ISO_BIWEEK,'LineStyle','-','LineWidth',2,'Color',cols(1,:))
xlim(XL)
set(gca,'FontSize',axFS,'FontName',font)
ylabel("δD (‰)",'FontSize',labFS)
ylim([-280 -150])

subplot(4,1,2)
plot(t,Mg_BIWEEK,'LineStyle','-','LineWidth',2,'Color',cols(2,:))
xlim(XL)
set(gca,'FontSize',axFS,'FontName',font)
ylabel("Mg (ppb)",'FontSize',labFS)
ylim([0 10])

subplot(4,1,3)
plot(t,MSA_BIWEEK,'LineStyle','-','LineWidth',2,'Color',cols(3,:))
xlim(XL)
set(gca,'FontSize',axFS,'FontName',font)
ylabel("MSA (ppb)",'FontSize',labFS)
ylim([0 10])

subplot(4,1,4)
plot(t,CFA_BIWEEK,'LineStyle','-','LineWidth',2,'Color',cols(4,:))
xlim(XL)
set(gca,'FontSize',axFS,'FontName',font)
ylabel("Dust (ppb)",'FontSize',labFS)
xlabel("Year",'FontSize',labFS)
ylim([0 1000])
for i = 1:4
    subplot(4,1,i)
    box on
end
exportgraphics(gcf,'Figures/SeasonalityPreAlignment.png','Resolution',res)



close all; figure('Position',[584 528 800 600])
subplot(4,1,1)
plot(t,ISO_BIWEEK(AI),'LineStyle','-','LineWidth',2,'Color',cols(1,:))
%yyaxis right
xlim(XL)
set(gca,'FontSize',axFS,'FontName',font)
ylabel("δD (‰)",'FontSize',labFS)
ylim([-280 -150])

subplot(4,1,2)
plot(t,Mg_BIWEEK(AI),'LineStyle','-','LineWidth',2,'Color',cols(2,:))
xlim(XL)
set(gca,'FontSize',axFS,'FontName',font)
ylabel("Mg (ppb)",'FontSize',labFS)
ylim([0 10])

subplot(4,1,3)
plot(t,MSA_BIWEEK(AI),'LineStyle','-','LineWidth',2,'Color',cols(3,:))
xlim(XL)
set(gca,'FontSize',axFS,'FontName',font)
ylabel("MSA (ppb)",'FontSize',labFS)
ylim([0 10])

subplot(4,1,4)
plot(t,CFA_BIWEEK(AI),'LineStyle','-','LineWidth',2,'Color',cols(4,:))
xlim(XL)
set(gca,'FontSize',axFS,'FontName',font)
ylabel("Dust (ppb)",'FontSize',labFS)
xlabel("Year",'FontSize',labFS)
ylim([0 1000])
for i = 1:4
    subplot(4,1,i)
    box on
end
exportgraphics(gcf,'Figures/SeasonalityPostAlignment.png','Resolution',res)



t = aligned(:,1);
col1 = "#c99206";
col2 = "#167a2a";

XL = [2000 2010];
close all; figure('Position',[584 528 800 600])
subplot(4,1,1); hold on
plot(t,ISO_BIWEEK,'LineStyle','-','LineWidth',1,'Color',col1)
xlim(XL)
set(gca,'FontSize',axFS,'FontName',font)
ylabel("δD (‰)",'FontSize',labFS)
ylim([-280 -150])

subplot(4,1,2); hold on
plot(t,Mg_BIWEEK,'LineStyle','-','LineWidth',1,'Color',col1)
xlim(XL)
set(gca,'FontSize',axFS,'FontName',font)
ylabel("Mg (ppb)",'FontSize',labFS)
ylim([0 10])

subplot(4,1,3); hold on
plot(t,MSA_BIWEEK,'LineStyle','-','LineWidth',1,'Color',col1)
xlim(XL)
set(gca,'FontSize',axFS,'FontName',font)
ylabel("MSA (ppb)",'FontSize',labFS)
ylim([0 10])

subplot(4,1,4); hold on
plot(t,CFA_BIWEEK,'LineStyle','-','LineWidth',1,'Color',col1)
xlim(XL)
set(gca,'FontSize',axFS,'FontName',font)
ylabel("Dust (ppb)",'FontSize',labFS)
xlabel("Year",'FontSize',labFS)
ylim([0 1000])
for i = 1:4
    subplot(4,1,i)
    box on
end


subplot(4,1,1)
plot(t,ISO_BIWEEK(AI),'LineStyle','-','LineWidth',2,'Color',col2)
%yyaxis right
xlim(XL)
set(gca,'FontSize',axFS,'FontName',font)
ylabel("δD (‰)",'FontSize',labFS)
ylim([-280 -150])

subplot(4,1,2)
plot(t,Mg_BIWEEK(AI),'LineStyle','-','LineWidth',2,'Color',col2)
xlim(XL)
set(gca,'FontSize',axFS,'FontName',font)
ylabel("Mg (ppb)",'FontSize',labFS)
ylim([0 10])

subplot(4,1,3)
plot(t,MSA_BIWEEK(AI),'LineStyle','-','LineWidth',2,'Color',col2)
xlim(XL)
set(gca,'FontSize',axFS,'FontName',font)
ylabel("MSA (ppb)",'FontSize',labFS)
ylim([0 10])

subplot(4,1,4)
plot(t,CFA_BIWEEK(AI),'LineStyle','-','LineWidth',2,'Color',col2)
xlim(XL)
set(gca,'FontSize',axFS,'FontName',font)
ylabel("Dust (ppb)",'FontSize',labFS)
xlabel("Year",'FontSize',labFS)
ylim([0 1000])
for i = 1:4
    subplot(4,1,i)
    box on
    set(gca,'XGrid','on')
end
subplot(4,1,2)
legend(["Original","Aligned"],'Location','northwest','FontSize',legFS)
exportgraphics(gcf,'Figures/SeasonalityAlignmentSubplot.png','Resolution',res)


%% Isotope top/bottom time vs middle time
close all; figure('Position',[525 564 1000 400])
plot(ISO.DEN13B.MidTime,ISO.DEN13B.Deuterium,'LineWidth',2)
hold on
[x,y] = getStepped(ISO.DEN13B.BotTime,ISO.DEN13B.TopTime,ISO.DEN13B.Deuterium);
plot(x,y,'LineWidth',2)
xlim([1995 2000])
set(gca,'FontSize',axFS,'FontName',font)
legend(["Middle time only","Top/bottom times"],'Location','northeast','FontSize',legFS)
ylabel("δD (‰)",'FontSize',labFS)
xlabel("Year",'FontSize',labFS)
xticks(1980:2010)
exportgraphics(gcf,'Figures/StepVsMid.png','Resolution',res)

