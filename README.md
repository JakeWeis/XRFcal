# XRFcal

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
%   variables required
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
