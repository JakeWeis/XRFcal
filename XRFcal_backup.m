function XRFcal(varargin)
% Calibration of element concentrations determined using a portable XRF
% analyser.
% 
% The calibration is based on linear regression models fitted through
% semi-quantitative pXRF data and corresponding high-accuracy ICP-MS data.
% Measurements were conducted on a wide range of marine sediment samples
% collected in the South Indian Ocean and Tasman Sea. Each sample was
% analysed both using a portable XRF analyser (Olympus Vanta Series) and an
% inductively-coupled plasma mass spectrometer (Thermo Fischer Scientific
% Element 2). The calibration data is stored in the "CalData" folder. In
% order to improve the accuracy of the calibration (especially of elements
% present in low quantities in the samples the calibration is currently
% based on) the calibration data sets should be extended with corresponding
% pXRF and ICP-MS data of different sediment types.
%
% As pXRF and ICP-MS target slightly different sets of elements, two
% different calibrations are being performed.
%
% Element-specific calibration:
% For each element that was analysed both by pXRF and ICP-MS, an element-
% specific calibration (ESC) line is being determined by fitting a linear
% regression model through the corresponding pXRF and ICP-MS element
% concentrations. The pXRF concentrations are then corrected for slope (m)
% and y-intercept (b) of the corresponding ESC.
%
% General calibration:
% Where an element was analysed only by pXRF, a general calibration (GC) is
% applied. The GC line is based on a linear regression model fitted through
% the entire data set of corresponding pXRF and ICP-MS concentrations (i.e.
% an average of all ESCs). The pXRF concentrations are then corrected for
% slope and y-intercept of the GC.
%
% PXRF analyses are conducted in two different modes, "geoChem" and "Soil",
% each targeting slightly different sets of elements. Therefore,
% ESC and GC are calculated separately for each pXRF mode.
%
% How to use the script:
% - Replace input data (pXRF concentrations) in the "input" folder. Make
%   sure to leave the structure of the excel file as is. Replace only the
%   the data (concentrations and, optionally, the sample-IDs, remove or add
%   rows if necessary). Make sure, that the pasted concentrations align
%   with the concentrations of the input file (row 1).
% - Set Matlab working directory to corresponding path of the "XRFcal"
%   folder
% - Run script by typing "XRFcal" into the Command Window (now input
%   variables required)
% - When prompted, select pXRF data to be calibrated (geoChem, Soil or
%   both)
% - Calibration output will be stored in a date/time subfolder within the
%   the "output" folder.
%
% Calibration Info:
% GC and ESC coefficients (slope, y-intercept and R2 error) can be found
% in the "Calibration_Info" subfolder (CalibrationInfo.txt).
% "CalibrationFigure_GEOCHEM" and "CalibrationFigure_SOIL" display
% regression plots and lines (ICP-MS vs. pXRF concentrations; green plots -
% R2>.9, orange - .5<R2<.9, red - R2<.5).

%% Import calibration raw data
ICPcaldat = readmatrix('CalData/ICPMS_CalibrationData.csv','Range',[1,2]);
XRFcaldat(:,:,1) = readmatrix('CalData/XRF_geoChem_CalibrationData.csv','Range',[1,2]);
XRFcaldat(:,:,2) = readmatrix('CalData/XRF_soil_CalibrationData.csv','Range',[1,2]);
fid = fopen('CalData/XRF_geoChem_CalibrationData.csv');
CalElementsXRF = strsplit(fgetl(fid), ';');
CalElementsXRF(1) = [];
fclose(fid);
fid = fopen('CalData/ICPMS_CalibrationData.csv');
CalElementsICPMS = strsplit(fgetl(fid), ';');
CalElementsICPMS(1) = [];
fclose(fid);
warning('off','curvefit:fit:noStartPoint');
warning('off','curvefit:cfit:subsasgn:coeffsClearingConfBounds');

for k = 1:2
    %% General Calibration (GC)
    %%% linear regression based on complete datasets
    GCfit = fittype('m*x+b',...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'m','b'});
    
    XRFvec = XRFcaldat(:,sum(isnan(XRFcaldat(:,:,k)))==0,k);
    XRFvec(:,13) = [];
    XRFvec = XRFvec(1:end)';
    ICPvec = ICPcaldat(:,sum(isnan(XRFcaldat(:,:,k)))==0);
    ICPvec(:,13) = [];
    ICPvec = ICPvec(1:end)';
    
    GC.coeffs = fit(ICPvec,XRFvec,GCfit);
    GC.XRFcal = GC.coeffs.m * ICPvec + GC.coeffs.b;
    GC.Rsq = 1 - sum((XRFvec - GC.XRFcal).^2)/sum((XRFvec - mean(XRFvec)).^2);
    GC.lineX = linspace(min(ICPvec),max(ICPvec),100);
    GC.lineY = GC.coeffs.m * GC.lineX+GC.coeffs.b;
    
    CalInfo_val = [GC.coeffs.m,GC.coeffs.b,GC.Rsq];
    CalInfo_str = {'General Calibration'};
    
%     figure(k)
%     p = subplot(5,5,1);
%     plot(ICPvec,XRFvec,'.k')
%     hold on
%     plot(GC.lineX,GC.lineY)
%     hold off
%     title('GC')
%     axis square
%     if GC.Rsq >= 0.9
%         p.Color = [0 1 0];
%     elseif GC.Rsq >= 0.5 && ESC.Rsq(n,1) < 0.9
%         p.Color = [0.75 0.75 0.3];
%     elseif GC.Rsq < 0.5
%         p.Color = [1 0.4 0.4];
%     elseif isnan(GC.Rsq)
%         p.Color = [1 1 1];
%         xticks([]);
%         yticks([]);
%     end
    
    %% Element-Specific Calibration (ESC)
    ESCfit = fittype('m*x + b',...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'m','b'});
    
    ESC.XRFcal = zeros(size(ICPcaldat));
    ESC.Rsq = zeros(size(ICPcaldat,2),1);
    ESC.lineX = zeros(size(ICPcaldat,2),2);
    ESC.lineY = zeros(size(ICPcaldat,2),2);
    
    for i = 1 : size(ICPcaldat,2)
        if sum(isnan(XRFcaldat(:,i,k))) == 0
            if sum(XRFcaldat(:,i,k)) > 0
                ESC.coeffs{i} = fit(ICPcaldat(:,i),XRFcaldat(:,i,k),ESCfit);
            elseif sum(XRFcaldat(:,i,k)) == 0
                ESC.coeffs{i}.m = 0;
                ESC.coeffs{i}.b = 0;
            end
        else
            ESC.coeffs{i}.m = NaN;
            ESC.coeffs{i}.b = NaN;
        end
        ESC.XRFcal(:,i) = ESC.coeffs{i}.m * ICPcaldat(:,i) + ESC.coeffs{i}.b;
        if ~(ESC.coeffs{i}.m == 0)
            ESC.Rsq(i,1) = 1 - sum((XRFcaldat(:,i,k) - ESC.XRFcal(:,i)).^2)/sum((XRFcaldat(:,i,k) - mean(XRFcaldat(:,i,k))).^2);
        else
            ESC.Rsq(i,1) = 0;
        end
        ESC.lineX(i,:) = [min(ICPcaldat(:,i)), max(ICPcaldat(:,i))];
        ESC.lineY(i,:) = ESC.coeffs{i}.m * ESC.lineX(i,:) + ESC.coeffs{i}.b;
        
        CalInfo_val(i+1,:) = [ESC.coeffs{i}.m,ESC.coeffs{i}.b,ESC.Rsq(i)];
        CalInfo_str{i+1,1} = [CalElementsXRF{i}(1:end-14),' (pXRF) vs. ',...
            CalElementsICPMS{i},' (ICP-MS)'];
        
%         figure(k)
%         p = subplot(5,5,n+1);
%         plot(ICPcaldat(:,n),XRFcaldat(:,n,k),'.k')
%         hold on
%         plot(ESC.lineX(n,:),ESC.lineY(n,:))
%         hold off
%         title(CalElementsXRF{n}(1:end-14))
%         axis square
%         if ESC.Rsq(n,1) >= 0.9
%             p.Color = [0 1 0];
%         elseif ESC.Rsq(n,1) >= 0.5 && ESC.Rsq(n,1) < 0.9
%             p.Color = [0.75 0.75 0.3];
%         elseif ESC.Rsq(n,1) < 0.5
%             p.Color = [1 0.4 0.4];
%         elseif isnan(ESC.Rsq(n,1))
%             p.Color = [1 1 1];
%             xticks([]);
%             yticks([]);
%         end
    end
    
    if k == 1
        GCvals_G = CalInfo_val(1,:);
        ESCvals_G = CalInfo_val(2:end,:);
    elseif k == 2
        GCvals_S = CalInfo_val(1,:);
        ESCvals_S = CalInfo_val(2:end,:);
    end
       
end
infoG_m = [GCvals_G(:,1);ESCvals_G(:,1)];
infoG_b = [GCvals_G(:,2);ESCvals_G(:,2)];
infoG_R = [GCvals_G(:,3);ESCvals_G(:,3)];
infoS_m = [GCvals_S(:,1);ESCvals_S(:,1)];
infoS_b = [GCvals_S(:,2);ESCvals_S(:,2)];
infoS_R = [GCvals_S(:,3);ESCvals_S(:,3)];
infoG_m = char({' ';'slope';' ';num2str(infoG_m,4)});
infoG_b = char({'GeoChem';'y-int';' ';num2str(infoG_b,4)});
infoG_R = char({' ';'Rsq';' ';num2str(infoG_R,4)});
infoS_m = char({' ';'slope';' ';num2str(infoS_m,4)});
infoS_b = char({'Soil';'y-int';' ';num2str(infoS_b,4)});
infoS_R = char({' ';'Rsq';' ';num2str(infoS_R,4)});
for i = 1 : size(infoG_m,1)
    colsep(i,1:3) = ' | ';
    colseps(i,1:2) = '  ';
end
CalInfo_str = char({'Calibration';' ';' ';char(CalInfo_str)});

info = [CalInfo_str,colsep,...
    infoG_m,colseps,infoG_b,colseps,infoG_R,colsep...
    infoS_m,colseps,infoS_b,colseps,infoS_R];
writematrix(info,'Calibration_Info/CalibrationInfo.txt','Delimiter','tab','FileType','text')


%% DATA CORRECTION
modes = varargin;
if isempty(modes)
    answer = questdlg('Choose pXRF data for calibration', ...
        'pXRF Calibration', ...
        'GeoChem & Soil','GeoChem','Soil',...
        'GeoChem & Soil');
    switch answer
        case 'GeoChem & Soil'
            modes = 'GS';
        case 'GeoChem'
            modes = 'G';
        case 'Soil'
            modes = 'S';
        case 'none (generate calibration info file only)'
            modes = 'info';
    end
end

if sum(strcmp(modes, 'GS')) == 1 || sum(strcmp(modes, 'SG')) == 1
    RawData_G(:,:) = readmatrix('input/INPUT_GEOCHEM.csv','Range',[1,2]);
    RawData_S(:,:) = readmatrix('input/INPUT_SOIL.csv','Range',[1,2]);
elseif sum(strcmp(modes, 'G')) == 1
    RawData_G = readmatrix('input/INPUT_GEOCHEM.csv','Range',[1,2]);
    RawData_S = [];
elseif sum(strcmp(modes, 'S')) == 1
    RawData_S = readmatrix('input/INPUT_SOIL.csv','Range',[1,2]);
    RawData_G = [];
elseif sum(strcmp(modes, 'info')) == 1
    RawData_G = [];
    RawData_S = [];
elseif isempty(modes)
    warndlg('No input data specified. Only calibration information file generated.')
    RawData_G = [];
    RawData_S = [];
end
fid = fopen('input/INPUT_GEOCHEM.csv');
Elements = strsplit(fgetl(fid), ';');
Elements(1) = [];
fclose(fid);

CalOutput = readtable('input/INPUT_GEOCHEM.csv');
ColumnNames = CalOutput.Properties.VariableNames;
SampleIDs = table2cell(CalOutput(:,1));
folder = join(string(round(clock)),'_');
mkdir(join(['output/',folder],''))

if ~(isempty(RawData_G))
    CalData_G = zeros(size(RawData_G));
    for i = 1 : size(RawData_G,2)
        if isempty(find(strcmp(Elements(i),CalElementsXRF),1))
            CalData_G(:,i) = (RawData_G(:,i) - GCvals_G(2)) / GCvals_G(1);
        else
            n = find(strcmp(Elements(i),CalElementsXRF),1);
            CalData_G(:,i) = (RawData_G(:,i) - ESCvals_G(n,2)) / ESCvals_G(n,1);
        end
    end
    CalTab_G = array2table(CalData_G,...
        'VariableNames',ColumnNames(2:end),'RowNames',SampleIDs)
    CalTab_G.Properties.DimensionNames(1) = {'SampleID'};
    writetable(CalTab_G,join(['output/',folder,'/OUTPUT_GEOCHEM.csv'],''),'WriteRowNames',true)
end

if ~(isempty(RawData_S))
    CalData_S = zeros(size(RawData_S));
    for i = 1 : size(RawData_S,2)
        if isempty(find(strcmp(Elements(i),CalElementsXRF),1))
            CalData_S(:,i) = (RawData_S(:,i) - GCvals_S(2)) / GCvals_S(1);
        else
            n = find(strcmp(Elements(i),CalElementsXRF),1);
            CalData_S(:,i) = (RawData_S(:,i) - ESCvals_S(n,2)) / ESCvals_S(n,1);
        end
    end
    CalTab_S = array2table(CalData_S,...
        'VariableNames',ColumnNames(2:end),'RowNames',SampleIDs)
    CalTab_S.Properties.DimensionNames(1) = {'SampleID'};
    writetable(CalTab_S,join(['output/',folder,'/OUTPUT_SOIL.csv'],''),'WriteRowNames',true)
end

end