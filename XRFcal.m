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

CalFig = 0;

for k = 1:2
    %% LINEAR REGRESSION MODEL: General Calibration (GC)
    %%% linear regression based on complete datasets
    XRFvec = XRFcaldat(:,sum(isnan(XRFcaldat(:,:,k)))==0,k);
    XRFvec(:,13) = [];
    XRFvec = XRFvec(1:end)';
    ICPvec = ICPcaldat(:,sum(isnan(XRFcaldat(:,:,k)))==0);
    ICPvec(:,13) = [];
    ICPvec = ICPvec(1:end)';
    
    GCmodel = fitlm(ICPvec,XRFvec);                     % GC CALIBRATION MODEL OUTPUT
    GC.coeffs.m = GCmodel.Coefficients.Estimate(2);     % slope
    GC.coeffs.b = GCmodel.Coefficients.Estimate(1);     % y-intercept
    GC.coeffs.mSE = GCmodel.Coefficients.SE(2);         % standard error of slope
    GC.coeffs.bSE = GCmodel.Coefficients.SE(1);         % standard error of y-intercept
    GC.Rsq = GCmodel.Rsquared.Adjusted;                 % R-squared of model
    
    CalInfo_val = [GC.coeffs.m,GC.coeffs.b,GC.coeffs.mSE,GC.coeffs.bSE,GC.Rsq];
    CalInfo_str = {'General Calibration'};
    
    % Calibration figure
    if CalFig == 1
        figure(k);clf
        GC.lineX = linspace(min(ICPvec),max(ICPvec),100);
        GC.lineY = GC.coeffs.m * GC.lineX + GC.coeffs.b;
        
        % Subplot setup: [gap-h gap-v],[mar-b mar-t],[mar-l mar-r]
        p = tight_subplot(5,5,[.05 .05],[.05 .05],[.05 .05]);
        axes(p(1))
        plot(ICPvec,XRFvec,'.k')
        hold on
        plot(GC.lineX,GC.lineY,'k')
        hold off
        title('GC')
        axis square
        if GC.Rsq >= 0.9
            p(1).Color = [0 1 0];
        elseif GC.Rsq >= 0.5 && ESC.Rsq(nobcGn,1) < 0.9
            p(1).Color = [0.75 0.75 0.3];
        elseif GC.Rsq < 0.5
            p(1).Color = [1 0.4 0.4];
        elseif isnan(GC.Rsq)
            p(1).Color = [1 1 1];
            xticks([]);
            yticks([]);
        end
    end
    
    %% LINEAR REGRESSION MODEL: Element-Specific Calibration (ESC)
    ESC.Rsq = zeros(size(ICPcaldat,2),1);
    ESC.lineX = zeros(size(ICPcaldat,2),2);
    ESC.lineY = zeros(size(ICPcaldat,2),2);
    for i = 1 : size(ICPcaldat,2)
        if sum(isnan(XRFcaldat(:,i,k))) == 0
            if sum(XRFcaldat(:,i,k)) > 0
                ESCmodel = fitlm(ICPcaldat(:,i),XRFcaldat(:,i,k));      % ESC CALIBRATION MODEL OUTPUT
                ESC.coeffs{i}.m = ESCmodel.Coefficients.Estimate(2);    % slope
                ESC.coeffs{i}.b = ESCmodel.Coefficients.Estimate(1);    % y-intercept
                ESC.coeffs{i}.mSE = ESCmodel.Coefficients.SE(2);        % standard error of slope
                ESC.coeffs{i}.bSE = ESCmodel.Coefficients.SE(1);        % standard error of y-intercept
                ESC.Rsq(i,1) = ESCmodel.Rsquared.Adjusted;              % R-squared of model
            elseif sum(XRFcaldat(:,i,k)) == 0
                ESC.coeffs{i}.m = 0;
                ESC.coeffs{i}.b = 0;
                ESC.coeffs{i}.mSE = 0;
                ESC.coeffs{i}.bSE = 0;
                ESC.Rsq(i,1) = 0;
            end
        else
            ESC.coeffs{i}.m = NaN;
            ESC.coeffs{i}.b = NaN;
            ESC.coeffs{i}.mSE = NaN;
            ESC.coeffs{i}.bSE = NaN;
            ESC.Rsq(i,1) = NaN;
        end

        CalInfo_val(i+1,:) = [ESC.coeffs{i}.m,ESC.coeffs{i}.b,ESC.coeffs{i}.mSE,ESC.coeffs{i}.bSE,ESC.Rsq(i)];
        CalInfo_str{i+1,1} = [CalElementsXRF{i}(1:end-14),' (pXRF) vs. ',CalElementsICPMS{i},' (ICP-MS)'];
        
        % Calibration figure
        if CalFig == 1
            figure(k)
            ESC.lineX(i,:) = [min(ICPcaldat(:,i)), max(ICPcaldat(:,i))];
            ESC.lineY(i,:) = ESC.coeffs{i}.m * ESC.lineX(i,:) + ESC.coeffs{i}.b;
            axes(p(i+1))
            plot(ICPcaldat(:,i),XRFcaldat(:,i,k),'.k')
            hold on
            plot(ESC.lineX(i,:),ESC.lineY(i,:),'k')
            hold off
            title(CalElementsXRF{i}(1:end-14))
            axis square
            if ESC.Rsq(i,1) >= 0.9
                p(i+1).Color = [0 1 0];
            elseif ESC.Rsq(i,1) >= 0.5 && ESC.Rsq(i,1) < 0.9
                p(i+1).Color = [0.75 0.75 0.3];
            elseif ESC.Rsq(i,1) < 0.5
                p(i+1).Color = [1 0.4 0.4];
            elseif isnan(ESC.Rsq(i,1))
                p(i+1).Color = [1 1 1];
                xticks([]);
                yticks([]);
            end
        end
    end
    
    if CalFig == 1
        FigFontResize(14,1.1,1.2)
        set(gcf, 'InvertHardcopy', 'off')
        if k == 1
            print(gcf,'Calibration_Info/CalibrationFigure_GEOCHEM','-dpng')
        else
            print(gcf,'Calibration_Info/CalibrationFigure_SOIL','-dpng')
        end
    end
    
    
    if k == 1
        GCvals_G = CalInfo_val(1,:);
        ESCvals_G = CalInfo_val(2:end,:);
    elseif k == 2
        GCvals_S = CalInfo_val(1,:);
        ESCvals_S = CalInfo_val(2:end,:);
    end
       
end

%% Calibration info text file
infoG_m = [GCvals_G(:,1);ESCvals_G(:,1)];
infoG_mSE = [GCvals_G(:,3);ESCvals_G(:,3)];
infoG_m = char({' ';' ';'Slope';' ';[strjust(num2str(infoG_m,'%1.4f'),'center'),repmat(' (+-',length(infoG_m),1),strjust(num2str(infoG_mSE,'%1.4f'),'center'),repmat(')',length(infoG_m),1)]});
infoG_m(1:3,:) = strjust(infoG_m(1:3,:),'center');
infoG_m([2,4],:) = '-';
infoG_m(find(isnan(ESCvals_G(:,1)))+5,:) = repmat(pad('--',size(infoG_m,2),'both'),2,1);
infoG_m(find(ESCvals_G(:,1)==0)+5,:) = repmat(pad('0 (+- 0)',size(infoG_m,2),'both'),length(find(ESCvals_G(:,1)==0)),1);

infoG_b = [GCvals_G(:,2);ESCvals_G(:,2)];
infoG_bSE = [GCvals_G(:,4);ESCvals_G(:,4)];
infoG_b = char({'GeoChem';' ';'Intercept';' ';[strjust(num2str(infoG_b,'%1.4f'),'center'),repmat(' (+-',length(infoG_b),1),strjust(num2str(infoG_bSE,'%1.4f'),'center'),repmat(')',length(infoG_b),1)]});
infoG_b(1:3,:) = strjust(infoG_b(1:3,:),'center');
infoG_b([2,4],:) = '-';
infoG_b(find(isnan(ESCvals_G(:,1)))+5,:) = repmat(pad('--',size(infoG_b,2),'both'),2,1);
infoG_b(find(ESCvals_G(:,1)==0)+5,:) = repmat(pad('0 (+- 0)',size(infoG_b,2),'both'),length(find(ESCvals_G(:,1)==0)),1);

infoG_R = [GCvals_G(:,5);ESCvals_G(:,5)];
infoG_R = strjust(char({' ';' ';'Rsq';' ';num2str(infoG_R,'%1.4f')}),'center');
infoG_R([2,4],:) = '-';
infoG_R(find(isnan(ESCvals_G(:,1)))+5,:) = repmat(pad('--',size(infoG_R,2),'both'),2,1);
infoG_R(find(isnan(ESCvals_G(:,1)))+5,:) = repmat(pad('--',size(infoG_R,2),'both'),2,1);
infoG_R(find(ESCvals_G(:,1)==0)+5,:) = repmat(pad('0',size(infoG_R,2),'both'),length(find(ESCvals_G(:,1)==0)),1);

infoS_m = [GCvals_S(:,1);ESCvals_S(:,1)];
infoS_mSE = [GCvals_S(:,3);ESCvals_S(:,3)];
infoS_m = char({' ';' ';'Slope';' ';[strjust(num2str(infoS_m,'%1.4f'),'center'),repmat(' (+-',length(infoS_m),1),strjust(num2str(infoS_mSE,'%1.4f'),'center'),repmat(')',length(infoS_m),1)]});
infoS_m(1:3,:) = strjust(infoS_m(1:3,:),'center');
infoS_m([2,4],:) = '-';
infoS_m(find(isnan(ESCvals_S(:,1)))+5,:) = repmat(pad('--',size(infoS_m,2),'both'),2,1);
infoS_m(find(ESCvals_S(:,1)==0)+5,:) = repmat(pad('0 (+- 0)',size(infoS_m,2),'both'),length(find(ESCvals_S(:,1)==0)),1);

infoS_b = [GCvals_S(:,2);ESCvals_S(:,2)];
infoS_bSE = [GCvals_S(:,4);ESCvals_S(:,4)];
infoS_b = char({'Soil';' ';'Intercept';' ';[strjust(num2str(infoS_b,'%1.4f'),'center'),repmat(' (+-',length(infoS_b),1),strjust(num2str(infoS_bSE,'%1.4f'),'center'),repmat(')',length(infoS_b),1)]});
infoS_b(1:3,:) = strjust(infoS_b(1:3,:),'center');
infoS_b([2,4],:) = '-';
infoS_b(find(isnan(ESCvals_S(:,1)))+5,:) = repmat(pad('--',size(infoS_b,2),'both'),2,1);
infoS_b(find(ESCvals_S(:,1)==0)+5,:) = repmat(pad('0 (+- 0)',size(infoS_b,2),'both'),length(find(ESCvals_S(:,1)==0)),1);

infoS_R = [GCvals_S(:,5);ESCvals_S(:,5)];
infoS_R = strjust(char({' ';' ';'Rsq';' ';num2str(infoS_R,'%1.4f')}),'center');
infoS_R([2,4],:) = '-';
infoS_R(find(isnan(ESCvals_S(:,1)))+5,:) = repmat(pad('--',size(infoS_R,2),'both'),2,1);
infoS_R(find(ESCvals_S(:,1)==0)+5,:) = repmat(pad('0',size(infoS_R,2),'both'),length(find(ESCvals_S(:,1)==0)),1);

Acolsep = repmat('  ||  ',size(infoG_m,1),1);
Bcolsep = repmat('    ',size(infoG_m,1),1);

CalInfo_str = strjust(char({' ';'Calibration';' ';' ';char(CalInfo_str)}),'center');CalInfo_str(4,:) = '-';

info = [CalInfo_str,Acolsep,...
    infoG_m,Bcolsep,infoG_b,Bcolsep,infoG_R,Acolsep...
    infoS_m,Bcolsep,infoS_b,Bcolsep,infoS_R];
writematrix(info,'Calibration_Info/CalibrationInfo.txt','Delimiter','tab','FileType','text')

%% DATA CORRECTION
mode = varargin;
if isempty(mode)
    % Choose calibration dlg
    fn = {'Best case','GeoChem & Soil','GeoChem','Soil','none (generate calibration info file)'};
    mode = listdlg('PromptString','Select data to be calibrated:',...
        'SelectionMode','single','ListString',fn,'ListSize',[200,70]);
    % Specify R2 threshold dlg
    if mode ~= 5
        prompt = {'In order to ignore calibrations of insufficient quality, specify a threshold R^2 value (between 0 and 1) below which to use the raw data instead of the calibrated data.'};
        dlgtitle = 'R squared threshold';
        definput = {'0.5'};
        opts.Interpreter = 'tex';
        opts.Resize = 'on';
        threshold = str2double(inputdlg(prompt,dlgtitle,[1 40],definput,opts));
    end
end

switch mode
    case {1,2}
    RawData_G(:,:) = readmatrix('input/INPUT_GEOCHEM.csv','Range',[1,2]);
    RawData_S(:,:) = readmatrix('input/INPUT_SOIL.csv','Range',[1,2]);
    case 3
    RawData_G = readmatrix('input/INPUT_GEOCHEM.csv','Range',[1,2]);
    RawData_S = [];
    case 4
    RawData_S = readmatrix('input/INPUT_SOIL.csv','Range',[1,2]);
    RawData_G = [];
    case 5
    warndlg('No input data specified. Only calibration information file generated.')
    RawData_G = [];
    RawData_S = [];
end

fid = fopen('input/INPUT_GEOCHEM.csv');
Elements = strsplit(fgetl(fid), ',');
Elements(1) = [];
fclose(fid);

CalOutput = readtable('input/INPUT_GEOCHEM.csv','PreserveVariableNames',1);
ColumnNames = CalOutput.Properties.VariableNames;
SampleIDs = table2cell(CalOutput(:,1));
[~,i,~] = unique(SampleIDs,'first');
IDdupl = find(~ismember(1:numel(SampleIDs),i));
if length(IDdupl) > 1
    for i = 1 : length(IDdupl)
        SampleIDs(IDdupl(i)) = {[SampleIDs{IDdupl(i)}(:)','(2)']};
    end
end

folder = join(string(round(clock)),'_');
if mode ~= 5
    mkdir(join(['output/',folder],''))
end

%% GeoChem calibration & error propagation
if mode == 1 || mode == 2 || mode == 3
    CalData_G = zeros(size(RawData_G));
    for i = find(contains(Elements,'Concentration'))
        if isempty(find(strcmp(Elements(i),CalElementsXRF),1))
            % GC
            CalData_G(:,i) = (RawData_G(:,i) - GCvals_G(2)) / GCvals_G(1);
            % Error propagation
            CalData_G(:,i+1) = abs(CalData_G(:,i)) .* sqrt(((sqrt(RawData_G(:,i+1).^2 + GCvals_G(4).^2)) ./ ((RawData_G(:,i) - GCvals_G(2)))).^2 + (GCvals_G(3) ./ GCvals_G(1)).^2);
        else
            nobcGn = find(strcmp(Elements(i),CalElementsXRF),1);
            % ESC
            CalData_G(:,i) = (RawData_G(:,i) - ESCvals_G(nobcGn,2)) / ESCvals_G(nobcGn,1);
            % Error propagation
            CalData_G(:,i+1) = abs(CalData_G(:,i)) .* sqrt(((sqrt(RawData_G(:,i+1).^2 + ESCvals_G(nobcGn,4).^2)) ./ ((RawData_G(:,i) - ESCvals_G(nobcGn,2)))).^2 + (ESCvals_G(nobcGn,3)./ESCvals_G(nobcGn,1)).^2);
        end
    end
    % Export calibrated dataset as csv file
    if mode ~=1
        CalTab_G = array2table(CalData_G,'VariableNames',Elements,'RowNames',SampleIDs);
        CalTab_G.Properties.DimensionNames(1) = {'SampleID'};
        writetable(CalTab_G,join(['output/',folder,'/OUTPUT_GEOCHEM.csv'],''),'WriteRowNames',true)
    end
end

%% Soil calibration & error propagation
if mode == 1 || mode == 2 || mode == 4
    CalData_S = zeros(size(RawData_S));
    for i = find(contains(Elements,'Concentration'))
        if isempty(find(strcmp(Elements(i),CalElementsXRF),1))
            % GC
            CalData_S(:,i) = (RawData_S(:,i) - GCvals_S(2)) / GCvals_S(1);
            % Error propagation
            CalData_S(:,i+1) = CalData_S(:,i) .* sqrt(((sqrt(RawData_S(:,i+1).^2 + GCvals_S(4).^2)) ./ ((RawData_S(:,i) - GCvals_S(2)))).^2 + (GCvals_S(3) ./ GCvals_S(1)).^2);
        else
            nobcGn = find(strcmp(Elements(i),CalElementsXRF),1);
            % ESC
            CalData_S(:,i) = (RawData_S(:,i) - ESCvals_S(nobcGn,2)) / ESCvals_S(nobcGn,1);
            % Error propagation
            CalData_S(:,i+1) = CalData_S(:,i) .* sqrt(((sqrt(RawData_S(:,i+1).^2 + ESCvals_S(nobcGn,4).^2)) ./ ((RawData_S(:,i) - ESCvals_S(nobcGn,2)))).^2 + (ESCvals_S(nobcGn,3)./ESCvals_S(nobcGn,1)).^2);
        end
    end
    % Export calibrated dataset as csv file
    if mode ~= 1
        CalTab_S = array2table(CalData_S,'VariableNames',Elements,'RowNames',SampleIDs);
        CalTab_S.Properties.DimensionNames(1) = {'SampleID'};
        writetable(CalTab_S,join(['output/',folder,'/OUTPUT_SOIL.csv'],''),'WriteRowNames',true)
    end
end

%% Compile best case calibrations & errors
if mode == 1
                                                                      % GeoChem best case where ...
    bcG = ESCvals_G(:,5) > threshold & ...                            % GeoChem R2 is larger than the threshold value ...
        (ESCvals_G(:,5) > ESCvals_S(:,5) | isnan(ESCvals_S(:,5)));    % ... [AND] GeoChem R2 is larger than Soil R2 [OR] Soil R2 is NaN.
    
                                                                      % Soil best case where ...
    bcS = ESCvals_S(:,5) > threshold & ...                            % ... Soil R2 is larger than the threshold value ...
        (ESCvals_S(:,5) > ESCvals_G(:,5) | isnan(ESCvals_G(:,5)));    % ... [AND] Soil R2 is larger than GeoChem R2 [OR] GeoChem R2 is NaN.
    
                                                                      % uncalibrated GeoChem data where ...
    nobcG = CalElementsXRF((bcG == 0 & bcS == 0) & ...                % ... no best case was applied ...
        (ESCvals_G(:,5) > ESCvals_S(:,5) | isnan(ESCvals_S(:,5))))    % ... [AND] GeoChem R2 is larger than Soil R2 [OR] Soil R2 is NaN.
    
                                                                      % uncalibrated Soil data where ...
    nobcS = CalElementsXRF((bcG == 0 & bcS == 0) & ...                % ... no best case was applied ...
        (isnan(ESCvals_G(:,5)) | ESCvals_S(:,5) > ESCvals_G(:,5) |... % [AND] GeoChem R2 is larger than Soil R2 [OR] Soil R2 is NaN ...
        (ESCvals_S(:,5) == 0 & ESCvals_G(:,5) == 0)))                 % [OR] either R2 is zero (applies to Ag and Cd)
    
    bcG = CalElementsXRF(bcG)
    bcS = CalElementsXRF(bcS)
    
    CalData_BC = NaN(size(CalData_S));
    BCinfo = cell(1,size(CalData_S,2));
    for i = 1 : length(bcG)
        bcGn(i) = find(strcmp(bcG(i),ColumnNames(2:end)),1);
        CalData_BC(:,bcGn(i)) = CalData_G(:,bcGn(i));
        CalData_BC(:,bcGn(i)+1) = CalData_G(:,bcGn(i)+1);
        BCinfo{bcGn(i)} = 'GeoChem ESC';
        BCinfo{bcGn(i)+1} = 'GeoChem ESC';
    end
    for i = 1 : length(bcS)
        bcSn(i) = find(strcmp(bcS(i),ColumnNames(2:end)),1);
        CalData_BC(:,bcSn(i)) = CalData_S(:,bcSn(i));
        CalData_BC(:,bcSn(i)+1) = CalData_S(:,bcSn(i)+1);
        BCinfo{bcSn(i)} = 'Soil ESC';
        BCinfo{bcSn(i)+1} = 'Soil ESC';
    end
    for i = 1 : length(nobcG)
        nobcGn(i) = find(strcmp(nobcG(i),ColumnNames(2:end)),1);
        CalData_BC(:,nobcGn(i)) = RawData_G(:,nobcGn(i));
        CalData_BC(:,nobcGn(i)+1) = RawData_G(:,nobcGn(i)+1);
        BCinfo{nobcGn(i)} = 'GeoChem uncalibrated';
        BCinfo{nobcGn(i)+1} = 'GeoChem uncalibrated';
    end
    for i = 1 : length(nobcS)
        nobcSn(i) = find(strcmp(nobcS(i),ColumnNames(2:end)),1);
        CalData_BC(:,nobcSn(i)) = RawData_S(:,nobcSn(i));
        CalData_BC(:,nobcSn(i)+1) = RawData_S(:,nobcSn(i)+1);
        BCinfo{nobcSn(i)} = 'Soil uncalibrated';
        BCinfo{nobcSn(i)+1} = 'Soil uncalibrated';
    end

    GCn = sum(isnan(CalData_BC)) == size(CalData_BC,1);
    GCn(47:end-2) = 0;
    CalData_BC(:,GCn) = CalData_G(:,GCn);
    BCinfo(GCn) = repmat({'GeoChem GC'},size(BCinfo(GCn)));
    
    GCn = sum(isnan(CalData_BC)) == size(CalData_BC,1);
    GCn([1:46,end-1:end]) = 0;
    CalData_BC(:,GCn) = CalData_S(:,GCn);
    BCinfo(GCn) = repmat({'Soil GC'},size(BCinfo(GCn)));
    
    % Export calibrated dataset as csv file
    CalTab_BC = array2table(CalData_BC,'VariableNames',Elements,'RowNames',SampleIDs');
    CalTab_BC.Properties.DimensionNames(1) = {'SampleID'};
    writetable(CalTab_BC,join(['output/',folder,'/OUTPUT_BEST-CASE.csv'],''),'WriteRowNames',true)
    
    CalTab_BCinfo = array2table(BCinfo,'VariableNames',Elements);
    writetable(CalTab_BCinfo,join(['output/',folder,'/OUTPUT_BEST-CASE-INFO.csv'],''),'WriteRowNames',true)
    
end

end