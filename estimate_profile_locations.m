% ------------------------------------------------------------------------------
% Estimation of not located profile positions using a method based on the
% "Terrain-following" method (Kaihe Yamazaki et al. paper
% https://doi.org/10.1029/2019JC015406).
%
% The implemented method is based on Kaihe Yamazaki's method that we improved:
%  - to avoid "dead end" issues
%  - to consider float in situ measurements.
%
% SYNTAX :
%   estimate_profile_locations(6902899)
%
% INPUT PARAMETERS :
%   varargin : WMO number of float to process
%
% OUTPUT PARAMETERS :
%
% EXAMPLES :
%
% SEE ALSO :
% AUTHORS  : Jean-Philippe Rannou (Altran)(jean-philippe.rannou@altran.com)
% ------------------------------------------------------------------------------
% RELEASES :
%   03/01/2022 - RNU - creation
% ------------------------------------------------------------------------------
function estimate_profile_locations(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONFIGURATION - START

% default list of floats to process
FLOAT_LIST_FILE_NAME = 'C:\Users\jprannou\_RNU\DecArgo_soft\lists\tmp.txt';
% FLOAT_LIST_FILE_NAME = 'C:\Users\jprannou\_RNU\DecArgo_soft\lists\liste_snapshot_202202.txt';
FLOAT_LIST_FILE_NAME = 'C:\Users\jprannou\_RNU\DecArgo_soft\lists\liste_snapshot_202202_prof_qc_8_9.txt';

% top directory of the NetCDF files
DIR_INPUT_NC_FILES = 'C:\Users\jprannou\_DATA\DATA_UNDER_ICE\IN\';
% DIR_INPUT_NC_FILES = 'D:\202202-ArgoData\coriolis\';

% directory of output files
DIR_OUTPUT_FILES = 'C:\Users\jprannou\_DATA\DATA_UNDER_ICE\OUT6\';

% directory to store the log file
DIR_LOG_FILE = 'C:\Users\jprannou\_RNU\DecArgo_soft\work\log\';

% GEBCO bathymetric file
GEBCO_FILE = 'C:\Users\jprannou\_RNU\_ressources\GEBCO_2021\GEBCO_2021.nc';

% max difference (in meters) between sea bottom and float parking drift to start
DIFF_DEPTH_TO_START = 100000;

% tolerance (in meters) used to compare float pressure and GEBCO depth
FLOAT_VS_BATHY_TOLERANCE = 10;

% tolerance (in meters) used to compare float pressure and GEBCO depth when the
% float grounded
FLOAT_VS_BATHY_TOLERANCE_FOR_GRD = 170;

% first range value (in kilometers)
FIRST_RANGE = 20;

% last range value (in kilometers)
LAST_RANGE = 20;

% range period (in kilometers)
RANGE_PERIOD = 5;

% CONFIGURATION - END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% global values initialization
init_global_values;

% specific global variables
global g_estProfLoc_version;
global g_estProfLoc_diffDepthToStart;
global g_estProfLoc_floatVsbathyTolerance;
global g_estProfLoc_floatVsbathyToleranceForGrd;
global g_estProfLoc_firstRange;
global g_estProfLoc_lastRange;
global g_estProfLoc_rangePeriod;

g_estProfLoc_version = '1.0';
g_estProfLoc_diffDepthToStart = DIFF_DEPTH_TO_START;
g_estProfLoc_floatVsbathyTolerance = FLOAT_VS_BATHY_TOLERANCE;
g_estProfLoc_floatVsbathyToleranceForGrd = FLOAT_VS_BATHY_TOLERANCE_FOR_GRD;
g_estProfLoc_firstRange = FIRST_RANGE;
g_estProfLoc_lastRange = LAST_RANGE;
g_estProfLoc_rangePeriod = RANGE_PERIOD;


% check inputs
if (nargin == 0)
   if ~(exist(FLOAT_LIST_FILE_NAME, 'file') == 2)
      fprintf('ERROR: File not found: %s\n', FLOAT_LIST_FILE_NAME);
      return
   end
end
if ~(exist(DIR_INPUT_NC_FILES, 'dir') == 7)
   fprintf('ERROR: Directory not found: %s\n', DIR_INPUT_NC_FILES);
   return
end
if ~(exist(GEBCO_FILE, 'file') == 2)
   fprintf('ERROR: File not found: %s\n', GEBCO_FILE);
   return
end

% get floats to process
if (nargin == 0)
   % floats to process come from default list
   fprintf('Floats from list: %s\n', FLOAT_LIST_FILE_NAME);
   floatList = load(FLOAT_LIST_FILE_NAME);
else
   % floats to process come from input parameters
   floatList = cell2mat(varargin);
end

% create and start log file recording
if (nargin == 0)
   [~, name, ~] = fileparts(FLOAT_LIST_FILE_NAME);
   name = ['_' name];
else
   name = sprintf('_%d', floatList);
end

% store the start time of the run
currentTime = datestr(now, 'yyyymmddTHHMMSSZ');

% create log directory
if ~(exist(DIR_LOG_FILE, 'dir') == 7)
   mkdir(DIR_LOG_FILE);
end

logFile = [DIR_LOG_FILE '/' 'estimate_profile_locations' name '_' currentTime '.log'];
diary(logFile);
tic;

% process floats of the list
nbFloats = length(floatList);
for idFloat = 1:nbFloats

   floatNum = floatList(idFloat);
   floatNumStr = num2str(floatNum);
   fprintf('%03d/%03d %s\n', idFloat, nbFloats, floatNumStr);

   % retrieve float data from NetCDF files
   floatData = get_float_data(floatNum, [DIR_INPUT_NC_FILES '/' floatNumStr '/']);
   if (isempty(floatData))
      fprintf('No profile location to estimate\n');
      continue
   end

   % process float data
   process_float_data(floatNum, floatData, [DIR_OUTPUT_FILES '/' floatNumStr '/'], GEBCO_FILE);

end

ellapsedTime = toc;
fprintf('done (Elapsed time is %.1f seconds)\n', ellapsedTime);

diary off;

return

% ------------------------------------------------------------------------------
% Process one float data.
%
% SYNTAX :
%  process_float_data(a_floatNum, a_floatData, a_outputDir, a_gebcoFilePathName)
%
% INPUT PARAMETERS :
%   a_floatNum          : float WMO number
%   a_floatData         : float data
%   a_outputDir         : output directory
%   a_gebcoFilePathName : GEBCO file path name
%
% OUTPUT PARAMETERS :
%   o_ncData : nc data
%
% EXAMPLES :
%
% SEE ALSO :
% AUTHORS  : Jean-Philippe Rannou (Altran)(jean-philippe.rannou@altran.com)
% ------------------------------------------------------------------------------
% RELEASES :
%   03/01/2022 - RNU - creation
% ------------------------------------------------------------------------------
function process_float_data(a_floatNum, a_floatData, a_outputDir, a_gebcoFilePathName)

% specific global variables
global g_estProfLoc_diffDepthToStart;

% QC flag values (numerical)
global g_decArgo_qcInterpolated;
global g_decArgo_qcMissing;


% define the sets of cycles to process
pos = ones(size(a_floatData.positionQc));
idF = find((a_floatData.positionQc == g_decArgo_qcInterpolated) | (a_floatData.positionQc == g_decArgo_qcMissing));
pos(idF) = 0;
startIdList = find(diff(pos) == -1);
stopIdList = find(diff(pos) == 1) + 1;

if (isempty(startIdList))
   return
end

if (length(startIdList) ~= length(stopIdList))
   if ((a_floatData.positionQc(end) == g_decArgo_qcInterpolated) || (a_floatData.positionQc(end) == g_decArgo_qcMissing))
      fprintf('ERROR: Float %d: inconsistent data (last profile location has QC = 8 or 9) - ignored\n', a_floatNum);
      return
   else
      fprintf('ERROR: Float %d: unknown reason (TO BE CHECKED) - ignored\n', a_floatNum);
      return
   end
end

% TEMP - START / look for QC = 9 loc or grounded cycles
% fprintf('@@FLOAT@@%d', a_floatNum);
%
% if (any((a_floatData.positionQc == g_decArgo_qcMissing)))
%    fprintf('@@QC=9 (%d)', length(find(a_floatData.positionQc == g_decArgo_qcMissing)));
% end
% if (any(((a_floatData.positionQc == g_decArgo_qcInterpolated) | ...
%       (a_floatData.positionQc == g_decArgo_qcMissing)) & ...
%       (a_floatData.grounded == 1)))
%    fprintf('@@GRD (%d)', length(find(((a_floatData.positionQc == g_decArgo_qcInterpolated) | ...
%       (a_floatData.positionQc == g_decArgo_qcMissing)) & ...
%       (a_floatData.grounded == 1))));
% end
%
% fprintf('\n');
%
% return
% TEMP - END

% create output directory
if ~(exist(a_outputDir, 'dir') == 7)
   mkdir(a_outputDir);
end

% interpolate anew profile locations (because some of them are missing
% (positionQc == g_decArgo_qcMissing) or are badly interpolated (5906033_074-084
% and 084-101))
paramJuld = get_netcdf_param_attributes('JULD');
idFv = find(a_floatData.juldLocation == paramJuld.fillValue);
a_floatData.juldLocation(idFv) = a_floatData.juld(idFv);

for idS = 1:length(startIdList)
   idStart = startIdList(idS);
   idStop = stopIdList(idS);

   % interpolate the locations
   [lonInter, latInter] = interpolate_between_2_locations(...
      a_floatData.juldLocation(idStart), a_floatData.longitude(idStart), a_floatData.latitude(idStart), ...
      a_floatData.juldLocation(idStop), a_floatData.longitude(idStop), a_floatData.latitude(idStop), ...
      a_floatData.juldLocation(idStart+1:idStop-1)');
   a_floatData.longitude(idStart+1:idStop-1) = lonInter';
   a_floatData.latitude(idStart+1:idStop-1) = latInter';
   a_floatData.positionQc(idStart+1:idStop-1) = g_decArgo_qcInterpolated;
end

% compute speeds
speed = nan(size(a_floatData.juldLocation));
for idC = 2:length(a_floatData.juldLocation)
   speed(idC) = ...
      100*distance_lpo([a_floatData.latitude(idC-1) a_floatData.latitude(idC)], ...
      [a_floatData.longitude(idC-1) a_floatData.longitude(idC)]) / ...
      ((a_floatData.juldLocation(idC)-a_floatData.juldLocation(idC-1))*86400);
end
a_floatData.speed = speed;

% retrieve GEBCO depth
a_floatData.gebcoDepth = get_gebco_depth(a_floatData.longitude, a_floatData.latitude, a_gebcoFilePathName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process the sets of cycles
for idS = 1:length(startIdList)
   idStart = startIdList(idS);
   idStop = stopIdList(idS);

   fprintf('   Processing one set of cycles: %d %03d-%03d\n', ...
      a_floatNum, ...
      a_floatData.cycleNumber(idStart), ...
      a_floatData.cycleNumber(idStop));

   % be sure RPP is not to far from bottom depth
   reducedFlag = 0;
   if (idStart > 1)
      for idC = idStart:idStop
         if (~isnan(a_floatData.rpp(idC)) ...
               && (a_floatData.gebcoDepth(idC) - a_floatData.rpp(idC) < g_estProfLoc_diffDepthToStart))
            break
         end
         idStart = idC;
         reducedFlag = 1;
      end
   end
   for idC = idStop:-1:idStart
      if (~isnan(a_floatData.rpp(idC)) ...
            && (a_floatData.gebcoDepth(idC) - a_floatData.rpp(idC) < g_estProfLoc_diffDepthToStart))
         break
      end
      idStop = idC;
      reducedFlag = 1;
   end
   if (idStart == idStop)
      fprintf('   Drifting depth too far from bathymetry - nothing done\n');
      continue
   elseif (reducedFlag == 1)
      fprintf('   Set of cycles reduced to: %d %03d-%03d\n', ...
         a_floatNum, ...
         a_floatData.cycleNumber(idStart), ...
         a_floatData.cycleNumber(idStop));
   end

   % process the set of cycles
   a_floatData = process_set_of_cycles(a_floatNum, a_floatData, idS, idStart, idStop, a_outputDir, a_gebcoFilePathName);

end

% write CSV report
print_csv_report(a_floatNum, a_floatData, a_outputDir);

return

% ------------------------------------------------------------------------------
% Process one set of cycles.
%
% SYNTAX :
%  [o_floatData] = process_set_of_cycles(a_floatNum, a_floatData, ...
%    a_setNum, a_idStart, a_idStop, a_outputDir, a_gebcoFilePathName)
%
% INPUT PARAMETERS :
%   a_floatNum          : float WMO number
%   a_floatData         : input float data
%   a_setNum            : number of the set of cycles
%   a_idStart           : first Id of the set of cycles
%   a_idStop            : last Id of the set of cycles
%   a_outputDir         : output directory
%   a_gebcoFilePathName : GEBCO file path name
%
% OUTPUT PARAMETERS :
%   o_ncData : nc data
%
% EXAMPLES :
%
% SEE ALSO :
% AUTHORS  : Jean-Philippe Rannou (Altran)(jean-philippe.rannou@altran.com)
% ------------------------------------------------------------------------------
% RELEASES :
%   03/01/2022 - RNU - creation
% ------------------------------------------------------------------------------
function [o_floatData] = process_set_of_cycles(a_floatNum, a_floatData, ...
   a_setNum, a_idStart, a_idStop, a_outputDir, a_gebcoFilePathName)

% output parameters initialization
o_floatData = a_floatData;

% specific global variables
global g_estProfLoc_floatVsbathyTolerance;
global g_estProfLoc_floatVsbathyToleranceForGrd;
global g_estProfLoc_firstRange;
global g_estProfLoc_lastRange;
global g_estProfLoc_rangePeriod;


% consider only data of the set
cycleNumber = a_floatData.cycleNumber(a_idStart:a_idStop);
juld = a_floatData.juld(a_idStart:a_idStop);
juldLocation = a_floatData.juldLocation(a_idStart:a_idStop);
latitude = a_floatData.latitude(a_idStart:a_idStop);
longitude = a_floatData.longitude(a_idStart:a_idStop);
profPresMax = a_floatData.profPresMax(a_idStart:a_idStop);
grounded = a_floatData.grounded(a_idStart:a_idStop);
groundedPres = a_floatData.groundedPres(a_idStart:a_idStop);
gebcoDepth = a_floatData.gebcoDepth(a_idStart:a_idStop);

% define depth constraint for each cycle
depthConstraint = interp1q([juld(1); juld(end)], [gebcoDepth(1); gebcoDepth(end)], juld')';
lastId = 1;
for idC = 2:length(cycleNumber)-1
   if (~isnan(profPresMax(idC)) && ((profPresMax(idC) > depthConstraint(idC)) || (grounded(idC) == 1)))
      depthConstraint(idC) = profPresMax(idC);
      depthConstraint(lastId:idC) = interp1q([juld(lastId); juld(idC)], [depthConstraint(lastId); depthConstraint(idC)], juld(lastId:idC)')';
      lastId = idC;
   end
end
if (lastId ~= 1)
   depthConstraint(lastId:end) = interp1q([juld(lastId); juld(end)], [depthConstraint(lastId); depthConstraint(end)], juld(lastId:end)')';
end

o_floatData.setNumber(a_idStart:a_idStop) = a_setNum;
o_floatData.depthConstraint(a_idStart:a_idStop) = depthConstraint;

% create the figure
close(findobj('Name', 'Estimate profile locations'));
warning off;

screenSize = get(0, 'ScreenSize');
figure('Name', 'Estimate profile locations', ...
   'Position', [1 screenSize(4)*(1/3) screenSize(3) screenSize(4)*(2/3)-90], ...
   'Color', 'w');

longitudeOri = longitude;
latitudeOri = latitude;

if (any(abs(diff(longitude)) > 180))
   id = find(longitude < 0);
   longitude(id) = longitude(id) + 360;
end

% angle of the normal to linear trajectory
[x, y] = latLon_2_xy(longitude([1 end]), latitude([1 end]));
x1 = x(1);
x2 = x(end);
y1 = y(1);
y2 = y(end);
den = x2 - x1 + eps;
tetaRad = atan((y2-y1)/den) + pi*(x2<x1) - pi/2;
tetaDeg = tetaRad*180/pi;
if (tetaDeg < 0)
   tetaDeg = tetaDeg + 360;
end
% fprintf('Teta: %.1f deg\n', tetaDeg);

% result arrays for forward and backward tries
resultF = nan(length(longitude), 2);
resultB = nan(length(longitude), 2);
resultF(1, 1) = longitude(1);
resultF(1, 2) = latitude(1);
resultB(1, 1) = longitude(end);
resultB(1, 2) = latitude(end);

% 2 loops: one for foreward, one for backward
for idLoop = 1:2

   if (idLoop == 2)
      longitude = fliplr(longitude);
      latitude = fliplr(latitude);
      depthConstraint = fliplr(depthConstraint);
   end

   % increase range until the path is found
   for range = g_estProfLoc_firstRange:g_estProfLoc_rangePeriod:g_estProfLoc_lastRange

      if (idLoop == 1)
         fprintf('   Trying forward with RANGE = %d', range);
      else
         fprintf('   Trying backward with RANGE = %d', range);
      end

      % create the map of locations to check
      nbCol = length(longitude);
      nbLig = (nbCol-1)*2*range+1;
      depthTabVal = nan(nbLig, nbCol);
      diffTabVal = nan(nbLig, nbCol);
      devTabFlag = ones(nbLig, nbCol);
      lonTabAll = nan(nbLig, nbCol);
      latTabAll = nan(nbLig, nbCol);
      for idC = 1:length(longitude)-1

         % create the set of locations on the search segment
         [lonTab, latTab] = get_loc_on_search_range(longitude([idC idC+1]), latitude([idC idC+1]), idC*range, tetaDeg);

         % retrieve location depth
         depthVal = get_gebco_depth(lonTab, latTab, a_gebcoFilePathName);

         depthFlag = ones(size(depthVal));
         diffVal = depthVal - depthConstraint(idC+1);
         if (grounded(idC+1) == 0)
            idOk = find(diffVal >= -g_estProfLoc_floatVsbathyTolerance);
         else
            idOk = find((diffVal >= -g_estProfLoc_floatVsbathyToleranceForGrd) & ...
               (diffVal <= g_estProfLoc_floatVsbathyToleranceForGrd));
         end
         depthFlag(idOk) = 0;

         depthTabVal((nbCol-(idC+1))*range+(1:length(depthVal)), idC+1) = depthVal;
         diffTabVal((nbCol-(idC+1))*range+(1:length(depthVal)), idC+1) = diffVal;
         devTabFlag((nbCol-(idC+1))*range+(1:length(depthVal)), idC+1) = depthFlag;
         lonTabAll((nbCol-(idC+1))*range+(1:length(depthVal)), idC+1) = lonTab;
         latTabAll((nbCol-(idC+1))*range+(1:length(depthVal)), idC+1) = latTab;
      end

      % try to find a path
      result = nan(length(longitude), 1);
      curId = (nbCol-1)*range + 1;
      idC = 1;
      done = 1;
      while (idC < length(longitude))
         searchId = curId-range:curId+range;
         idToCheck = find(devTabFlag(searchId, idC+1) == 0);
         if (~isempty(idToCheck))
            idToCheck = searchId(idToCheck);
            [~, minId] = min(abs(diffTabVal(idToCheck, idC+1)));
            curId = idToCheck(minId);
            devTabFlag((devTabFlag(:, idC+1) == 2), idC+1) = 3;
            devTabFlag(curId, idC+1) = 2;
            result(idC+1) = curId;
            idC = idC + 1;
         else
            if (idC <= 2)
               done = 0;
               break
            end
            idC = idC - 1;
            curId = result(idC);
         end
      end

      if (idLoop == 1)
         dir = 'Foreward';
         dir2 = '1_foreward';
      else
         dir = 'Backward';
         dir2 = '2_backward';
      end
      if (done == 1)
         koOk = 'OK';
      else
         koOk = 'KO';
      end
      fprintf(' - %s\n', koOk);
      label = sprintf('Float: %d - Cycles: %03d to %03d - %s - Range %d km - %s', ...
         a_floatNum, ...
         cycleNumber(1), ...
         cycleNumber(end), ...
         dir, ...
         range, ...
         koOk);
      pngFileName = sprintf('%d_%03d-%03d_%s_range_%d_%s.png', ...
         a_floatNum, ...
         cycleNumber(1), ...
         cycleNumber(end), ...
         dir2, ...
         range, ...
         koOk);

      % display result of the try

      % arrays to store legend information
      legendPlots = [];
      legendLabels = [];

      idDone = find(devTabFlag == 2);
      [lonMin, lonMax, latMin, latMax] = compute_geo_extrema( ...
         [], [longitudeOri lonTabAll(idDone)'], [latitudeOri latTabAll(idDone)'], 0);
      [elevC, lonC , latC] = get_gebco_elev_zone(lonMin, lonMax, latMin, latMax, '');

      cla;

      m_proj('mercator', 'latitudes', [latMin latMax], 'longitudes', [lonMin lonMax]);
      m_grid('box', 'fancy', 'tickdir', 'out', 'linestyle', 'none');
      hold on;

      isobath = -unique(round(depthConstraint));
      isobath = min(isobath):100:max(isobath);
      if (length(isobath) == 1)
         isobath = [isobath isobath];
      end
      [contourMatrix, contourHdl] = m_contour(lonC, latC, elevC, isobath, 'c');
      if (~isempty(contourMatrix))
         legendPlots = [legendPlots contourHdl];
         legendLabels = [legendLabels {'depth constraint isobath'}];
      end

      m_line([longitude(1) longitude(end)], [latitude(1) latitude(end)], 'linestyle', '-', 'visible', 'on');

      title(label, 'FontSize', 14);

      plotHdl = m_plot(longitude(1), latitude(1), 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'Markersize', 4);
      if (~isempty(plotHdl))
         legendPlots = [legendPlots plotHdl];
         legendLabels = [legendLabels {'starting location'}];
      end

      for idC = 1:length(longitude)-1
         lonT = lonTabAll(:, idC+1);
         latT = latTabAll(:, idC+1);

         lonL = lonT(~isnan(lonT));
         latL = latT(~isnan(latT));
         m_line([lonL(1) lonL(end)], [latL(1) latL(end)], 'linestyle', '-', 'visible', 'on');

         idNotChecked = find(devTabFlag(:, idC+1) == 0);
         plotHdl = m_plot(lonT(idNotChecked), latT(idNotChecked), 'h', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'Markersize',  4);
         if (~isempty(plotHdl))
            if (~any(strcmp(legendLabels, 'eligible locations')))
               legendPlots = [legendPlots plotHdl];
               legendLabels = [legendLabels {'eligible locations'}];
            end
         end

         idFailed = find(devTabFlag(:, idC+1) == 3);
         plotHdl = m_plot(lonT(idFailed), latT(idFailed), 'h', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'Markersize',  4);
         if (~isempty(plotHdl))
            if (~any(strcmp(legendLabels, 'failed locations')))
               legendPlots = [legendPlots plotHdl];
               legendLabels = [legendLabels {'failed locations'}];
            end
         end

         idDone = find(devTabFlag(:, idC+1) == 2);
         plotHdl = m_plot(lonT(idDone), latT(idDone), 'h', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'Markersize',  4);
         if (~isempty(plotHdl))
            if (~any(strcmp(legendLabels, 'final locations')))
               legendPlots = [legendPlots plotHdl];
               legendLabels = [legendLabels {'final locations'}];
            end
         end

         if (done)
            if (idLoop == 1)
               resultF(idC+1, 1) = lonT(idDone);
               resultF(idC+1, 2) = latT(idDone);
            else
               resultB(idC+1, 1) = lonT(idDone);
               resultB(idC+1, 2) = latT(idDone);
            end
         end
      end

      % plot legend
      legend(legendPlots, legendLabels, 'Location', 'NorthEastOutside', 'Tag', 'Legend');

      if (done || (range == g_estProfLoc_lastRange))
         print('-dpng', [a_outputDir '/' pngFileName]);
      end

      %       fprintf('Press any key ...');
      %       pause
      %       fprintf('\n');

      if (done)
         break
      end
   end
end

% plot final trajectory
if (done)

   label = sprintf('Float: %d - Cycles: %03d to %03d - Final trajectory\n', ...
      a_floatNum, ...
      cycleNumber(1), ...
      cycleNumber(end));

   pngFileName = sprintf('%d_%03d-%03d_3_final.png', ...
      a_floatNum, ...
      cycleNumber(1), ...
      cycleNumber(end));

   % arrays to store legend information
   legendPlots = [];
   legendLabels = [];

   [lonMin, lonMax, latMin, latMax] = compute_geo_extrema( ...
      [], [longitudeOri resultF(:, 1)' resultB(:, 1)'], [latitudeOri resultF(:, 2)' resultB(:, 2)'], 0);
   [elevC, lonC , latC] = get_gebco_elev_zone(lonMin, lonMax, latMin, latMax, '');

   cla;

   m_proj('mercator', 'latitudes', [latMin latMax], 'longitudes', [lonMin lonMax]);
   m_grid('box', 'fancy', 'tickdir', 'out', 'linestyle', 'none');
   hold on;

   isobath = -unique(round(depthConstraint));
   isobath = min(isobath):100:max(isobath);
   if (length(isobath) == 1)
      isobath = [isobath isobath];
   end
   m_contour(lonC, latC, elevC, isobath, 'c');

   m_line([longitude(1) longitude(end)], [latitude(1) latitude(end)], 'linestyle', '-', 'visible', 'on');

   title(label, 'FontSize', 14);

   for idC = 1:length(longitude)-1
      plotHdl = m_plot(resultF(idC+1, 1), resultF(idC+1, 2), 'o', 'Markersize', 3, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
      if (~isempty(plotHdl))
         if (~any(strcmp(legendLabels, 'foreward locations')))
            legendPlots = [legendPlots plotHdl];
            legendLabels = [legendLabels {'foreward locations'}];
         end
      end

      plotHdl = m_plot(resultB(idC+1, 1), resultB(idC+1, 2), 'o', 'Markersize',  6, 'MarkerEdgeColor', 'r');
      if (~isempty(plotHdl))
         if (~any(strcmp(legendLabels, 'backward locations')))
            legendPlots = [legendPlots plotHdl];
            legendLabels = [legendLabels {'backward locations'}];
         end
      end
   end

   trajLon = nan(length(longitude), 1);
   trajLat = nan(length(longitude), 1);
   nb1 = ceil(length(longitude)/2);
   nb2 = nb1;
   if (mod(length(longitude), 2) ~= 0)
      nb2 = nb2 - 1;
   end
   trajLon(1:nb1) = resultF(1:nb1, 1);
   trajLat(1:nb1) = resultF(1:nb1, 2);
   trajLon(nb1+1:end) = flipud(resultB(1:nb2, 1));
   trajLat(nb1+1:end) = flipud(resultB(1:nb2, 2));
   lineHdl = m_line(trajLon, trajLat, 'linestyle', '-', 'color', 'r', 'visible', 'on');
   if (~isempty(lineHdl))
      legendPlots = [legendPlots lineHdl];
      legendLabels = [legendLabels {'merged trajectory'}];
   end

   % plot legend
   legend(legendPlots, legendLabels, 'Location', 'NorthEastOutside', 'Tag', 'Legend');

   print('-dpng', [a_outputDir '/' pngFileName]);


   %    fprintf('Press any key ...');
   %    pause
   %    fprintf('\n');
end

% store output parameters

if (done)

   speed = nan(size(juld));
   for idC = 2:length(juld)
      speed(idC) = ...
         100*distance_lpo([trajLat(idC-1) trajLat(idC)], [trajLon(idC-1) trajLon(idC)]) / ...
         ((juld(idC)-juld(idC-1))*86400);
   end

   forwardLat = resultF(:, 2);
   forwardLon = resultF(:, 1);
   if (any(forwardLon > 180))
      id = find(forwardLon > 180);
      forwardLon(id) = forwardLon(id) - 360;
   end
   backwardLat = resultB(:, 2);
   backwardLon = resultB(:, 1);
   backwardLat = flipud(backwardLat);
   backwardLon = flipud(backwardLon);
   if (any(backwardLon > 180))
      id = find(backwardLon > 180);
      backwardLon(id) = backwardLon(id) - 360;
   end
   if (any(trajLon > 180))
      id = find(trajLon > 180);
      trajLon(id) = trajLon(id) - 360;
   end

   o_floatData.forwardLat(a_idStart:a_idStop) = forwardLat;
   o_floatData.forwardLon(a_idStart:a_idStop) = forwardLon;
   o_floatData.forwardGebcoDepth(a_idStart:a_idStop) = get_gebco_depth(forwardLon, forwardLat, a_gebcoFilePathName);
   o_floatData.backwardLat(a_idStart:a_idStop) = backwardLat;
   o_floatData.backwardLon(a_idStart:a_idStop) = backwardLon;
   o_floatData.backwardGebcoDepth(a_idStart:a_idStop) = get_gebco_depth(backwardLon, backwardLat, a_gebcoFilePathName);
   o_floatData.trajLat(a_idStart:a_idStop) = trajLat;
   o_floatData.trajLon(a_idStart:a_idStop) = trajLon;
   o_floatData.speedEst(a_idStart:a_idStop) = speed;
end

return

% ------------------------------------------------------------------------------
% Retrieve relevent data from float NetCDf files.
%
% SYNTAX :
%  [o_ncData] = get_float_data(a_floatNum, a_ncFileDir)
%
% INPUT PARAMETERS :
%   a_floatNum  : float WMO number
%   a_ncFileDir : float nc files directory
%
% OUTPUT PARAMETERS :
%   o_ncData : nc data
%
% EXAMPLES :
%
% SEE ALSO :
% AUTHORS  : Jean-Philippe Rannou (Altran)(jean-philippe.rannou@altran.com)
% ------------------------------------------------------------------------------
% RELEASES :
%   03/01/2022 - RNU - creation
% ------------------------------------------------------------------------------
function [o_floatData] = get_float_data(a_floatNum, a_ncFileDir)

% output parameters initialization
o_floatData = [];

% global measurement codes
global g_MC_Grounded;

% QC flag values (char)
global g_decArgo_qcStrGood;
global g_decArgo_qcStrProbablyGood;

% QC flag values (numerical)
global g_decArgo_qcInterpolated;
global g_decArgo_qcMissing;


paramJuld = get_netcdf_param_attributes('JULD');
paramPres = get_netcdf_param_attributes('PRES');
paramLat = get_netcdf_param_attributes('LATITUDE');
paramLon = get_netcdf_param_attributes('LONGITUDE');
floatData = get_data_init_struct;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% retrieve information from mono PROF files

profDirName = [a_ncFileDir '/profiles/'];

floatFiles = dir([profDirName '/' sprintf('*%d_*.nc', a_floatNum)]);
for idFile = 1:length(floatFiles)
   floatFileName = floatFiles(idFile).name;
   if (floatFileName(1) == 'B')
      continue
   end
   floatFilePathName = [profDirName '/' floatFileName];

   % retrieve data from file
   wantedVars = [ ...
      {'FORMAT_VERSION'} ...
      {'CYCLE_NUMBER'} ...
      {'DIRECTION'} ...
      {'DATA_MODE'} ...
      {'JULD'} ...
      {'JULD_QC'} ...
      {'JULD_LOCATION'} ...
      {'LATITUDE'} ...
      {'LONGITUDE'} ...
      {'POSITION_QC'} ...
      {'PRES'} ...
      {'PRES_ADJUSTED'} ...
      {'CONFIG_MISSION_NUMBER'} ...
      ];
   ncData = get_data_from_nc_file(floatFilePathName, wantedVars);

   formatVersion = get_data_from_name('FORMAT_VERSION', ncData)';
   formatVersion = strtrim(formatVersion);
   cycleNumber = get_data_from_name('CYCLE_NUMBER', ncData);
   direction = get_data_from_name('DIRECTION', ncData);
   dataMode = get_data_from_name('DATA_MODE', ncData);
   juld = get_data_from_name('JULD', ncData);
   juldQc = get_data_from_name('JULD_QC', ncData);
   juldLocation = get_data_from_name('JULD_LOCATION', ncData);
   latitude = get_data_from_name('LATITUDE', ncData);
   longitude = get_data_from_name('LONGITUDE', ncData);
   positionQc = get_data_from_name('POSITION_QC', ncData);
   pres = get_data_from_name('PRES', ncData);
   presAdjusted = get_data_from_name('PRES_ADJUSTED', ncData);
   configMissionNumber = get_data_from_name('CONFIG_MISSION_NUMBER', ncData);

   % check the file format version
   if (~strcmp(formatVersion, '3.1'))
      fprintf('ERROR: Input mono prof file (%s) is expected to be of 3.1 format version (but FORMAT_VERSION = %s) - ignored\n', ...
         floatFileName, formatVersion);
      continue
   end

   % check data consistency
   if ((length(unique(cycleNumber)) > 1) || (length(unique(direction)) > 1) || ...
         (length(unique(juld)) > 1) || (length(unique(juldQc)) > 1) || ...
         (length(unique(juldLocation)) > 1) || (length(unique(latitude)) > 1) || ...
         (length(unique(longitude)) > 1) || (length(unique(positionQc)) > 1) || ...
         (length(unique(configMissionNumber)) > 1))

      fprintf('ERROR: Inconsistent data in file: %s - ignored\n', floatFileName);
      continue
   end

   if (all(juld == paramJuld.fillValue))
      fprintf('WARNING: Not dated profile in file: %s - ignored\n', floatFileName);
      continue
   end

   floatData.cycleNumber = [floatData.cycleNumber unique(cycleNumber)];
   if (unique(direction) == 'D')
      direct = 1;
   else
      direct = 2;
   end
   floatData.direction = [floatData.direction direct];
   floatData.juld = [floatData.juld unique(juld)];
   floatData.juldQc = [floatData.juldQc str2double(unique(juldQc))];
   floatData.juldLocation = [floatData.juldLocation unique(juldLocation)];
   floatData.latitude = [floatData.latitude unique(latitude)];
   floatData.longitude = [floatData.longitude unique(longitude)];
   floatData.positionQc = [floatData.positionQc str2double(unique(positionQc))];
   floatData.configMissionNumber = [floatData.configMissionNumber unique(configMissionNumber)];

   % compute profile max pressure
   presMax = -1;
   for idProf = 1:length(dataMode)
      if (dataMode(idProf) == 'R')
         presVal = pres(:, idProf);
      else
         presVal = presAdjusted(:, idProf);
      end
      presVal(presVal == paramPres.fillValue) = [];
      if (max(presVal) > presMax)
         presMax = max(presVal);
      end
   end
   floatData.profPresMax = [floatData.profPresMax presMax];
end

if ~(any((floatData.positionQc == g_decArgo_qcInterpolated) | (floatData.positionQc == g_decArgo_qcMissing)))
   return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% retrieve information from META file

floatFiles = dir([a_ncFileDir '/' sprintf('%d_meta.nc', a_floatNum)]);
if (isempty(floatFiles))
   fprintf('ERROR: Meta-data file not found - ignored\n');
   return
end

floatFilePathName = [a_ncFileDir '/' floatFiles(1).name];

% retrieve information from file
wantedVars = [ ...
   {'FORMAT_VERSION'} ...
   {'LAUNCH_CONFIG_PARAMETER_NAME'} ...
   {'LAUNCH_CONFIG_PARAMETER_VALUE'} ...
   {'CONFIG_PARAMETER_NAME'} ...
   {'CONFIG_PARAMETER_VALUE'} ...
   {'CONFIG_MISSION_NUMBER'} ...
   ];
ncData = get_data_from_nc_file(floatFilePathName, wantedVars);

formatVersion = get_data_from_name('FORMAT_VERSION', ncData)';
formatVersion = strtrim(formatVersion);
launchConfigParamName = get_data_from_name('LAUNCH_CONFIG_PARAMETER_NAME', ncData);
launchConfigValue = get_data_from_name('LAUNCH_CONFIG_PARAMETER_VALUE', ncData);
configParamName = get_data_from_name('CONFIG_PARAMETER_NAME', ncData);
configValue = get_data_from_name('CONFIG_PARAMETER_VALUE', ncData);
configMissionNumberMeta = get_data_from_name('CONFIG_MISSION_NUMBER', ncData);

% check the file format version
if (~strcmp(formatVersion, '3.1'))
   fprintf('ERROR: Input meta file (%s) is expected to be of 3.1 format version (but FORMAT_VERSION = %s) - ignored\n', ...
      floatFiles(1).name, formatVersion);
   return
end

% retrieve the needed configuration parameters
[~, nParam] = size(launchConfigParamName);
launchConfigName = [];
for idParam = 1:nParam
   launchConfigName{end+1} = deblank(launchConfigParamName(:, idParam)');
end
[~, nParam] = size(configParamName);
configName = [];
for idParam = 1:nParam
   configName{end+1} = deblank(configParamName(:, idParam)');
end

% process retrieved data
parkP = -1;
idF = find(strcmp('CONFIG_ParkPressure_dbar', launchConfigName(:)) == 1, 1);
if (~isempty(idF) && (launchConfigValue(idF) ~= paramPres.fillValue))
   parkP = launchConfigValue(idF);
end
parkPres = ones(size(configMissionNumberMeta))*parkP;
idF = find(strcmp('CONFIG_ParkPressure_dbar', configName(:)) == 1, 1);
if (~isempty(idF))
   parkPres = configValue(idF, :);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% retrieve information from TRAJ file

floatFileName = '';
floatFiles = dir([a_ncFileDir '/' sprintf('%d_*traj.nc', a_floatNum)]);
for idFile = 1:length(floatFiles)
   if (any(floatFiles(idFile).name == 'B'))
      continue
   end
   floatFileName = floatFiles(idFile).name;
end

if (isempty(floatFileName))
   fprintf('ERROR: Trajectory file not found - ignored\n');
   return
end

floatFilePathName = [a_ncFileDir '/' floatFileName];

% retrieve information from file
wantedVars = [ ...
   {'FORMAT_VERSION'} ...
   {'TRAJECTORY_PARAMETERS'} ...
   {'JULD'} ...
   {'JULD_QC'} ...
   {'LATITUDE'} ...
   {'LONGITUDE'} ...
   {'POSITION_ACCURACY'} ...
   {'POSITION_QC'} ...
   {'CYCLE_NUMBER'} ...
   {'MEASUREMENT_CODE'} ...
   {'PRES'} ...
   {'PRES_ADJUSTED'} ...
   {'TRAJECTORY_PARAMETER_DATA_MODE'} ...
   {'GROUNDED'} ...
   {'REPRESENTATIVE_PARK_PRESSURE'} ...
   {'CONFIG_MISSION_NUMBER'} ...
   {'CYCLE_NUMBER_INDEX'}, ...
   {'DATA_MODE'} ...
   ];
ncData = get_data_from_nc_file(floatFilePathName, wantedVars);
formatVersion = get_data_from_name('FORMAT_VERSION', ncData)';
formatVersion = strtrim(formatVersion);
trajParam = get_data_from_name('TRAJECTORY_PARAMETERS', ncData);
juld = get_data_from_name('JULD', ncData);
juldQc = get_data_from_name('JULD_QC', ncData);
latitude = get_data_from_name('LATITUDE', ncData);
longitude = get_data_from_name('LONGITUDE', ncData);
positionAccuracy = get_data_from_name('POSITION_ACCURACY', ncData);
positionQc = get_data_from_name('POSITION_QC', ncData);
cycleNumber = get_data_from_name('CYCLE_NUMBER', ncData);
measCode = get_data_from_name('MEASUREMENT_CODE', ncData);
pres = get_data_from_name('PRES', ncData);
presAdjusted = get_data_from_name('PRES_ADJUSTED', ncData);
trajParamDataMode = get_data_from_name('TRAJECTORY_PARAMETER_DATA_MODE', ncData);
grounded = get_data_from_name('GROUNDED', ncData);
rpp = get_data_from_name('REPRESENTATIVE_PARK_PRESSURE', ncData);
configMissionNumber = get_data_from_name('CONFIG_MISSION_NUMBER', ncData);
cycleNumberIndex = get_data_from_name('CYCLE_NUMBER_INDEX', ncData);
dataMode = get_data_from_name('DATA_MODE', ncData);

% check the file format version
if (~ismember(formatVersion, [{'3.1'} {'3.2'}]))
   fprintf('ERROR: Input trajectory file (%s) is expected to be of 3.1 format version (but FORMAT_VERSION = %s)\n', ...
      floatFileName, formatVersion);
   return
end

% add cycle #0 location if any or launch location otherwise
idLoc0 = find((cycleNumber == 0) & ...
   (latitude ~= paramLat.fillValue) & ...
   (longitude ~= paramLon.fillValue) & ...
   (positionAccuracy ~= 'I') & ...
   ((positionQc == g_decArgo_qcStrGood) | (positionQc == g_decArgo_qcStrProbablyGood)));
if (~isempty(idLoc0))
   [~, idMax] = max(juld(idLoc0));
   juld0 = juld(idLoc0(idMax));
   juldQc0 = juldQc(idLoc0(idMax));
   longitude0 = longitude(idLoc0(idMax));
   latitude0 = latitude(idLoc0(idMax));
   positionQc0 = positionQc(idLoc0(idMax));
else
   idLoc0 = find(cycleNumber == -1);
   juld0 = juld(idLoc0);
   juldQc0 = juldQc(idLoc0);
   longitude0 = longitude(idLoc0);
   latitude0 = latitude(idLoc0);
   positionQc0 = positionQc(idLoc0);
end

floatData.cycleNumber = [0 floatData.cycleNumber];
floatData.direction = [0 floatData.direction];
floatData.juld = [juld0 floatData.juld];
floatData.juldQc = [str2double(juldQc0) floatData.juldQc];
floatData.juldLocation = [juld0 floatData.juldLocation];
floatData.latitude = [latitude0 floatData.latitude];
floatData.longitude = [longitude0 floatData.longitude];
floatData.positionQc = [str2double(positionQc0) floatData.positionQc];
floatData.configMissionNumber = [-1 floatData.configMissionNumber];
floatData.profPresMax = [0 floatData.profPresMax];

floatData.speed = nan(size(floatData.cycleNumber));
floatData.rpp = nan(size(floatData.cycleNumber));
floatData.grounded = zeros(size(floatData.cycleNumber));
floatData.groundedPres = nan(size(floatData.cycleNumber));
floatData.gebcoDepth = nan(size(floatData.cycleNumber));
floatData.setNumber = nan(size(floatData.cycleNumber));
floatData.depthConstraint = nan(size(floatData.cycleNumber));
floatData.forwardLat = nan(size(floatData.cycleNumber));
floatData.forwardLon = nan(size(floatData.cycleNumber));
floatData.forwardGebcoDepth = nan(size(floatData.cycleNumber));
floatData.backwardLat = nan(size(floatData.cycleNumber));
floatData.backwardLon = nan(size(floatData.cycleNumber));
floatData.backwardGebcoDepth = nan(size(floatData.cycleNumber));
floatData.trajLat = nan(size(floatData.cycleNumber));
floatData.trajLon = nan(size(floatData.cycleNumber));
floatData.speedEst = nan(size(floatData.cycleNumber));

if (strcmp(formatVersion, '3.2') && any(grounded == 'Y'))
   for idP = 1:size(trajParam, 2)
      paramName = strtrim(trajParam(:, idP)');
      if (strcmp(paramName, 'PRES'))
         presParamId = idP;
         break
      end
   end
end

% get RPP (Representative Parking Pressure) and GROUNDED (flag and pressure if
% any)
cycleNumberList = unique(floatData.cycleNumber);
for idCy = 1:length(cycleNumberList)
   cyNum = cycleNumberList(idCy);
   if (cyNum == 0)
      continue
   end
   idForCy = find(floatData.cycleNumber == cyNum);

   rppVal = rpp(cycleNumberIndex == cyNum);
   if (isempty(rppVal))
      floatData.rpp(idForCy) = nan;
   elseif (rppVal ~= paramPres.fillValue)
      floatData.rpp(idForCy) = rppVal;
   else
      rppVal = parkPres(configMissionNumberMeta == configMissionNumber(cycleNumberIndex == cyNum));
      floatData.rpp(idForCy) = rppVal;
   end

   groundedFlag = grounded(cycleNumberIndex == cyNum);
   if (groundedFlag == 'Y')
      floatData.grounded(idForCy) = 1;

      idGrd = find((cycleNumber == cyNum) & (measCode == g_MC_Grounded));
      if (~isempty(idGrd))
         if (strcmp(formatVersion, '3.1'))
            if (dataMode(cycleNumberIndex == cyNum) == 'R')
               presVal = pres;
            else
               presVal = presAdjusted;
            end
            presGrd = presVal(idGrd);
         elseif (strcmp(formatVersion, '3.2'))
            presGrd = [];
            for idG = 1:length(idGrd)
               if (trajParamDataMode(presParamId, idGrd(idG)) == 'R')
                  presGrd = [presGrd pres(idGrd(idG))];
               else
                  presGrd = [presGrd presAdjusted(idGrd(idG))];
               end
            end
         end
         presGrd(presGrd == paramPres.fillValue) = [];
         if (~isempty(presGrd))
            floatData.groundedPres(idForCy) = min(presGrd);
         end
      end
   end
end

% specific to float anomalies
if (ismember(a_floatNum, [6901880]))
   switch a_floatNum
      case 6901880
         floatData.profPresMax(floatData.cycleNumber == 14) = nan;
         floatData.rpp(floatData.cycleNumber == 14) = nan;
         floatData.grounded(floatData.cycleNumber == 14) = nan;
         floatData.groundedPres(floatData.cycleNumber == 14) = nan;
   end
end

% sort the data in chronological order
[~, idSort] = sort(floatData.juld);
floatData.cycleNumber = floatData.cycleNumber(idSort);
floatData.direction = floatData.direction(idSort);
floatData.juld = floatData.juld(idSort);
floatData.juldQc = floatData.juldQc(idSort);
floatData.juldLocation = floatData.juldLocation(idSort);
floatData.latitude = floatData.latitude(idSort);
floatData.longitude = floatData.longitude(idSort);
floatData.positionQc = floatData.positionQc(idSort);
floatData.speed = floatData.speed(idSort);
floatData.configMissionNumber = floatData.configMissionNumber(idSort);
floatData.profPresMax = floatData.profPresMax(idSort);

floatData.rpp = floatData.rpp(idSort);
floatData.grounded = floatData.grounded(idSort);
floatData.groundedPres = floatData.groundedPres(idSort);
floatData.gebcoDepth = floatData.gebcoDepth(idSort);
floatData.setNumber = floatData.setNumber(idSort);
floatData.depthConstraint = floatData.depthConstraint(idSort);
floatData.forwardLat = floatData.forwardLat(idSort);
floatData.forwardLon = floatData.forwardLon(idSort);
floatData.forwardGebcoDepth = floatData.forwardGebcoDepth(idSort);
floatData.backwardLat = floatData.backwardLat(idSort);
floatData.backwardLon = floatData.backwardLon(idSort);
floatData.backwardGebcoDepth = floatData.backwardGebcoDepth(idSort);
floatData.trajLat = floatData.trajLat(idSort);
floatData.trajLon = floatData.trajLon(idSort);
floatData.speedEst = floatData.speedEst(idSort);

% output data
o_floatData = floatData;

return

% ------------------------------------------------------------------------------
% Print estimated profile locations in CSV file report.
%
% SYNTAX :
%  print_csv_report(a_floatNum, a_floatData, a_outputDir)
%
% INPUT PARAMETERS :
%   a_floatNum  : float WMO number
%   a_floatData : input float data
%   a_outputDir : output directory
%
% OUTPUT PARAMETERS :
%
% EXAMPLES :
%
% SEE ALSO :
% AUTHORS  : Jean-Philippe Rannou (Altran)(jean-philippe.rannou@altran.com)
% ------------------------------------------------------------------------------
% RELEASES :
%   03/01/2022 - RNU - creation
% ------------------------------------------------------------------------------
function print_csv_report(a_floatNum, a_floatData, a_outputDir)

% specific global variables
global g_estProfLoc_version;
global g_estProfLoc_diffDepthToStart;
global g_estProfLoc_floatVsbathyTolerance;
global g_estProfLoc_floatVsbathyToleranceForGrd;


% CSV file creation
outputFileName = [a_outputDir '/estimate_profile_locations_' num2str(a_floatNum) '.csv'];
fidOut = fopen(outputFileName, 'wt');
if (fidOut == -1)
   fprintf('ERROR: Unable to create CSV output file: %s\n', outputFileName);
   return
end

% print file header
header = ['WMO;CyNum;Dir;Juld;JuldQC;JuldLoc;Lat;Lon;PosQC;Speed;ProfPresMax;' ...
   'Rpp;Grd;GrdPres;GebcoDepth;SetNum;DepthConstraint;' ...
   'ForwLat;ForwLon;ForwGebcoDepth;ForwDiffDepth;' ...
   'BackwLat;BackLon;BackwGebcoDepth;BackDiffDepth;TrajLat;TrajLon;SpeedEst;' ...
   ';DIFF_DEPTH_TO_START;FLOAT_VS_BATHY_TOLERANCE;FLOAT_VS_BATHY_TOLERANCE_FOR_GRD;TOOL_VERSION'];
fprintf(fidOut, '%s\n', header);

for idC = 1:length(a_floatData.cycleNumber)
   fprintf(fidOut, ...
      '%d;%d;%d;%s;%d;%s;%.3f;%.3f;%d;%.3f;%.1f;%.1f;%d;%.1f;%.1f;%d;%.1f;%.3f;%.3f;%.1f;%.1f;%.3f;%.3f;%.1f;%.1f;%.3f;%.3f;%.3f;;%d;%d;%d;%s\n', ...
      a_floatNum, ...
      a_floatData.cycleNumber(idC), ...
      a_floatData.direction(idC), ...
      julian_2_gregorian_dec_argo(a_floatData.juld(idC)), ...
      a_floatData.juldQc(idC), ...
      julian_2_gregorian_dec_argo(a_floatData.juldLocation(idC)), ...
      a_floatData.latitude(idC), ...
      a_floatData.longitude(idC), ...
      a_floatData.positionQc(idC), ...
      a_floatData.speed(idC), ...
      a_floatData.profPresMax(idC), ...
      a_floatData.rpp(idC), ...
      a_floatData.grounded(idC), ...
      a_floatData.groundedPres(idC), ...
      a_floatData.gebcoDepth(idC), ...
      a_floatData.setNumber(idC), ...
      a_floatData.depthConstraint(idC), ...
      a_floatData.forwardLat(idC), ...
      a_floatData.forwardLon(idC), ...
      a_floatData.forwardGebcoDepth(idC), ...
      a_floatData.forwardGebcoDepth(idC)-a_floatData.depthConstraint(idC), ...
      a_floatData.backwardLat(idC), ...
      a_floatData.backwardLon(idC), ...
      a_floatData.backwardGebcoDepth(idC), ...
      a_floatData.backwardGebcoDepth(idC)-a_floatData.depthConstraint(idC), ...
      a_floatData.trajLat(idC), ...
      a_floatData.trajLon(idC), ...
      a_floatData.speedEst(idC), ...
      g_estProfLoc_diffDepthToStart, ...
      g_estProfLoc_floatVsbathyTolerance, ...
      g_estProfLoc_floatVsbathyToleranceForGrd, ...
      g_estProfLoc_version ...
      );
end

fclose(fidOut);

return

% ------------------------------------------------------------------------------
% Get data from name in a {var_name}/{var_data} list.
%
% SYNTAX :
%  [o_dataValues] = get_data_from_name(a_dataName, a_dataList)
%
% INPUT PARAMETERS :
%   a_dataName : name of the data to retrieve
%   a_dataList : {var_name}/{var_data} list
%
% OUTPUT PARAMETERS :
%   o_dataValues : concerned data
%
% EXAMPLES :
%
% SEE ALSO :
% AUTHORS  : Jean-Philippe Rannou (Altran)(jean-philippe.rannou@altran.com)
% ------------------------------------------------------------------------------
% RELEASES :
%   06/12/2018 - RNU - creation
% ------------------------------------------------------------------------------
function [o_dataValues] = get_data_from_name(a_dataName, a_dataList)

% output parameters initialization
o_dataValues = [];

idVal = find(strcmp(a_dataName, a_dataList(1:2:end)) == 1, 1);
if (~isempty(idVal))
   o_dataValues = a_dataList{2*idVal};
end

return

% ------------------------------------------------------------------------------
% Get the basic structure to store a profile information.
%
% SYNTAX :
%  [o_profStruct] = get_profile_init_struct( ...
%    a_cycleNum, a_profNum, a_phaseNum, a_PrimarySamplingProfileFlag)
%
% INPUT PARAMETERS :
%   a_cycleNum                    : cycle number
%   a_profNum                     : profile number
%   a_phaseNum                    : phase number
%   a_PrimarySamplingProfileFlag  : 1 if it is a primary sampling profile,
%                                   0 otherwise
%
% OUTPUT PARAMETERS :
%   o_profStruct : profile initialized structure
%
% EXAMPLES :
%
% SEE ALSO :
% AUTHORS  : Jean-Philippe Rannou (Altran)(jean-philippe.rannou@altran.com)
% ------------------------------------------------------------------------------
% RELEASES :
%   02/25/2013 - RNU - creation
% ------------------------------------------------------------------------------
function [o_dataStruct] = get_data_init_struct

o_dataStruct = struct( ...
   'cycleNumber', [], ...
   'direction', [], ...
   'juld', [], ...
   'juldQc', [], ...
   'juldLocation', [], ...
   'latitude', [], ...
   'longitude', [], ...
   'positionQc', [], ...
   'speed', [], ...
   'configMissionNumber', [], ...
   'profPresMax', [], ...
   'rpp', [], ...
   'grounded', [], ...
   'groundedPres', [], ...
   'gebcoDepth', [], ...
   'setNumber', [], ...
   'depthConstraint', [], ...
   'forwardLat', [], ...
   'forwardLon', [], ...
   'forwardGebcoDepth', [], ...
   'backwardLat', [], ...
   'backwardLon', [], ...
   'backwardGebcoDepth', [], ...
   'trajLat', [], ...
   'trajLon', [], ...
   'speedEst', [] ...
   );

return

% ------------------------------------------------------------------------------
% Get GEBCO depth associated to a list of locations.
%
% SYNTAX :
%  [o_depth] = get_gebco_depth(a_lon, a_lat, a_gebcoPathFileName)
%
% INPUT PARAMETERS :
%   a_lon               : latitudes
%   a_lat               : longitudes
%   a_gebcoPathFileName : GEBCO path file name
%
% OUTPUT PARAMETERS :
%   o_depth : depths
%
% EXAMPLES :
%
% SEE ALSO :
% AUTHORS  : Jean-Philippe Rannou (Altran)(jean-philippe.rannou@altran.com)
% ------------------------------------------------------------------------------
% RELEASES :
%   03/01/2022 - RNU - creation
% ------------------------------------------------------------------------------
function [o_depth] = get_gebco_depth(a_lon, a_lat, a_gebcoPathFileName)

[elevOri] = get_gebco_elev_point(a_lon, a_lat, a_gebcoPathFileName);
elev = mean(elevOri, 2);
if (any(isnan(elev)))
   idNan = find(isnan(elev));
   for idL = idNan'
      elev(idL) = mean(elevOri(idL, ~isnan(elevOri(idL, :))));
   end
end
o_depth = -elev';

return

% ------------------------------------------------------------------------------
% Retrieve the surounding elevations of a list of locations from the GEBCO 2019
% file.
%
% SYNTAX :
%  [o_elev] = get_gebco_elev_point(a_lon, a_lat, a_gebcoFileName)
%
% INPUT PARAMETERS :
%   a_lon           : list of location longitudes
%   a_lat           : list of location atitudes
%   a_gebcoFileName : GEBCO 2019 file path name
%
% OUTPUT PARAMETERS :
%   o_elev : surrounding elevations of each location
%            (size(o_elev) = [length(a_lon) 4]
%             4 elevations are generally provided [elevSW elevNW elevSE elevNE]
%             when only 1 or 2 are provided other ones are set to NaN)
%
% EXAMPLES :
%
% SEE ALSO :
% AUTHORS  : Jean-Philippe Rannou (Altran)(jean-philippe.rannou@altran.com)
% ------------------------------------------------------------------------------
% RELEASES :
%   04/29/2020 - RNU - creation
% ------------------------------------------------------------------------------
function [o_elev] = get_gebco_elev_point(a_lon, a_lat, a_gebcoFileName)

% output parameters initialization
o_elev = nan(length(a_lon), 4);


% check inputs
if (a_lon < -180)
   fprintf('ERROR: get_gebco_elev_point: input lon < -180\n');
   return
end
if (a_lon >= 360)
   fprintf('ERROR: get_gebco_elev_point: input lat >= 360\n');
   return
end
if (a_lat < -90)
   fprintf('ERROR: get_gebco_elev_point: input lat < -90\n');
   return
elseif (a_lat > 90)
   fprintf('ERROR: get_gebco_elev_point: input lat > 90\n');
   return
end

if (a_lon >= 180)
   a_lon = a_lon - 360;
end

% check GEBCO file exists
if ~(exist(a_gebcoFileName, 'file') == 2)
   fprintf('ERROR: GEBCO file not found (%s)\n', a_gebcoFileName);
   return
end

% open NetCDF file
fCdf = netcdf.open(a_gebcoFileName, 'NC_NOWRITE');
if (isempty(fCdf))
   fprintf('RTQC_ERROR: Unable to open NetCDF input file: %s\n', a_gebcoFileName);
   return
end

lonVarId = netcdf.inqVarID(fCdf, 'lon');
latVarId = netcdf.inqVarID(fCdf, 'lat');
elevVarId = netcdf.inqVarID(fCdf, 'elevation');

lon = netcdf.getVar(fCdf, lonVarId);
lat = netcdf.getVar(fCdf, latVarId);
minLon = min(lon);
maxLon = max(lon);

for idP = 1:length(a_lat)

   if (isnan(a_lat(idP)) || isnan(a_lon(idP)))
      continue
   end

   idLigStart = find(lat <= a_lat(idP), 1, 'last');
   if (isempty(idLigStart))
      idLigStart = 1;
   end
   idLigEnd = find(lat >= a_lat(idP), 1, 'first');
   if (isempty(idLigEnd))
      idLigEnd = length(lat);
   end
   %    latVal = lat(fliplr(idLigStart:idLigEnd));

   % a_lon(idP) is in the [-180, 180[ interval
   % it can be in 3 zones:
   % case 1: [-180, minLon[
   % case 2: [minLon, maxLon]
   % case 3: ]maxLon, -180[
   if ((a_lon(idP) >= minLon) && (a_lon(idP) <= maxLon))
      % case 2
      idColStart = find(lon <= a_lon(idP), 1, 'last');
      idColEnd = find(lon >= a_lon(idP), 1, 'first');

      elev = nan(length(idLigStart:idLigEnd), length(idColStart:idColEnd));
      for idL = idLigStart:idLigEnd
         elev(end-(idL-idLigStart), :) = netcdf.getVar(fCdf, elevVarId, fliplr([idL-1 idColStart-1]), fliplr([1 length(idColStart:idColEnd)]))';
      end

      %       lonVal = lon(idColStart:idColEnd);
   elseif (a_lon(idP) < minLon)
      % case 1
      elev1 = nan(length(idLigStart:idLigEnd), 1);
      for idL = idLigStart:idLigEnd
         elev1(end-(idL-idLigStart), :) = netcdf.getVar(fCdf, elevVarId, fliplr([idL-1 length(lon)-1]), fliplr([1 1]))';
      end

      %       lonVal1 = lon(end);

      elev2 = nan(length(idLigStart:idLigEnd), 1);
      for idL = idLigStart:idLigEnd
         elev2(end-(idL-idLigStart), :) = netcdf.getVar(fCdf, elevVarId, fliplr([idL-1 0]), fliplr([1 1]))';
      end

      %       lonVal2 = lon(1) + 360;

      elev = cat(2, elev1, elev2);
      %       lonVal = cat(1, lonVal1, lonVal2);
      clear elev1 elev2
   elseif (a_lon(idP) > maxLon)
      % case 3
      elev1 = nan(length(idLigStart:idLigEnd), 1);
      for idL = idLigStart:idLigEnd
         elev1(end-(idL-idLigStart), :) = netcdf.getVar(fCdf, elevVarId, fliplr([idL-1 length(lon)-1]), fliplr([1 1]))';
      end

      %       lonVal1 = lon(end);

      elev2 = nan(length(idLigStart:idLigEnd), 1);
      for idL = idLigStart:idLigEnd
         elev2(end-(idL-idLigStart), :) = netcdf.getVar(fCdf, elevVarId, fliplr([idL-1 0]), fliplr([1 1]))';
      end

      %       lonVal2 = lon(1) + 360;

      elev = cat(2, elev1, elev2);
      %       lonVal = cat(1, lonVal1, lonVal2);
      clear elev1 elev2
   end

   if (~isempty(elev))
      if (size(elev, 1) == 2)
         if (size(elev, 2) == 2)
            o_elev(idP, 1) = elev(2, 1);
            o_elev(idP, 2) = elev(1, 1);
            o_elev(idP, 3) = elev(2, 2);
            o_elev(idP, 4) = elev(1, 2);
         else
            o_elev(idP, 1) = elev(2);
            o_elev(idP, 2) = elev(1);
         end
      else
         if (size(elev, 2) == 2)
            o_elev(idP, 1) = elev(1, 1);
            o_elev(idP, 3) = elev(1, 2);
         else
            o_elev(idP, 1) = elev;
         end
      end
   end

   clear elev
end

netcdf.close(fCdf);

clear lon lat

return

% ------------------------------------------------------------------------------
% Convert lon/lat to x/y.
%
% SYNTAX :
%  [o_x, o_y] = latLon_2_xy(a_lon, a_lat)
%
% INPUT PARAMETERS :
%   a_lon : latitudes
%   a_lat : longitudes
%
% OUTPUT PARAMETERS :
%   o_x : x
%   o_y : y
%
% EXAMPLES :
%
% SEE ALSO :
% AUTHORS  : Jean-Philippe Rannou (Altran)(jean-philippe.rannou@altran.com)
% ------------------------------------------------------------------------------
% RELEASES :
%   03/01/2022 - RNU - creation
% ------------------------------------------------------------------------------
function [o_x, o_y] = latLon_2_xy(a_lon, a_lat)

% switch in km into a local reference frame
xPosOri = a_lon(2);
yPosOri = a_lat(2);

valCos = cosd(a_lat(2))*1.852*60;
o_x = (a_lon - xPosOri) .* valCos;
o_y = (a_lat - yPosOri)*1.852*60;

return

% ------------------------------------------------------------------------------
% Create the set of locations on the search segment.
% One location per km, i.e. 2*range+1 locations on the search segment.
%
% SYNTAX :
%  [o_lon, o_lat] = get_loc_on_search_range(a_lon, a_lat, a_range, a_angle)
%
% INPUT PARAMETERS :
%   a_lon   : latitudes on the linear trajectory
%   a_lat   : longitudes on the linear trajectory
%   a_range : range dimention
%   a_angle : angle of the normal to linear trajectory
%
% OUTPUT PARAMETERS :
%   o_lon : longitudes on the search segment
%   o_lat : latitudes on the search segment
%
% EXAMPLES :
%
% SEE ALSO :
% AUTHORS  : Jean-Philippe Rannou (Altran)(jean-philippe.rannou@altran.com)
% ------------------------------------------------------------------------------
% RELEASES :
%   03/01/2022 - RNU - creation
% ------------------------------------------------------------------------------
function [o_lon, o_lat] = get_loc_on_search_range(a_lon, a_lat, a_range, a_angle)

% switch in km into a local reference frame
xPosOri = a_lon(2);
yPosOri = a_lat(2);

valCos = cosd(a_lat(2))*1.852*60;

xBis = (1:a_range).*cosd(a_angle);
yBis = (1:a_range).*sind(a_angle);

xTer = (1:a_range)*cosd(a_angle+180);
yTer = (1:a_range)*sind(a_angle+180);

lonBis = xBis./valCos + xPosOri;
latBis = yBis./(1.852*60) + yPosOri;

lonTer = xTer./valCos + xPosOri;
latTer = yTer./(1.852*60) + yPosOri;

o_lon = [fliplr(lonTer) a_lon(2) lonBis];
o_lat = [fliplr(latTer) a_lat(2) latBis];

return





% icicicicic
% ------------------------------------------------------------------------------
% Initialize global default values.
%
% SYNTAX :
%  init_global_values
%
% INPUT PARAMETERS :
%
% OUTPUT PARAMETERS :
%
% EXAMPLES :
%
% SEE ALSO :
% AUTHORS  : Jean-Philippe Rannou (Altran)(jean-philippe.rannou@altran.com)
% ------------------------------------------------------------------------------
% RELEASES :
%   01/02/2010 - RNU - creation
% ------------------------------------------------------------------------------
function init_global_values

% QC flag values (numerical)
global g_decArgo_qcInterpolated;
global g_decArgo_qcMissing;

% QC flag values (char)
global g_decArgo_qcStrGood;
global g_decArgo_qcStrProbablyGood;

% global measurement codes
global g_MC_Grounded;


% QC flag values (numerical)
g_decArgo_qcInterpolated = 8;
g_decArgo_qcMissing = 9;

% QC flag values (char)
g_decArgo_qcStrGood = '1';
g_decArgo_qcStrProbablyGood = '2';

% global measurement codes
g_MC_Grounded = 901;

return

% ------------------------------------------------------------------------------
% Get Argo attributes for a given parameter.
%
% SYNTAX :
%  [o_attributeStruct] = get_netcdf_param_attributes(a_paramName)
%
% INPUT PARAMETERS :
%   a_paramName : parameter name
%
% OUTPUT PARAMETERS :
%   o_attributeStruct : parameter associated attributes
%
% EXAMPLES :
%
% SEE ALSO :
% AUTHORS  : Jean-Philippe Rannou (Altran)(jean-philippe.rannou@altran.com)
% ------------------------------------------------------------------------------
% RELEASES :
%   02/25/2013 - RNU - creation
% ------------------------------------------------------------------------------
function [o_attributeStruct] = get_netcdf_param_attributes(a_paramName)

switch (a_paramName)

   case 'JULD'
      o_attributeStruct = struct('name', 'JULD', ...
         'longName', 'Julian day (UTC) of each measurement relative to REFERENCE_DATE_TIME', ...
         'standardName', 'time', ...
         'units', 'days since 1950-01-01 00:00:00 UTC', ...
         'conventions', 'Relative julian days with decimal part (as parts of day)', ...
         'fillValue', double(999999), ...
         'axis', 'T', ...
         'paramType', '', ...
         'paramNcType', 'NC_DOUBLE', ...
         'adjAllowed', 0);

   case 'LATITUDE'
      o_attributeStruct = struct('name', 'LATITUDE', ...
         'longName', 'Latitude of each location', ...
         'standardName', 'latitude', ...
         'units', 'degree_north', ...
         'fillValue', double(99999), ...
         'validMin', double(-90), ...
         'validMax', double(90), ...
         'axis', 'Y', ...
         'paramType', '', ...
         'paramNcType', 'NC_DOUBLE', ...
         'adjAllowed', 0);

   case 'LONGITUDE'
      o_attributeStruct = struct('name', 'LONGITUDE', ...
         'longName', 'Longitude of each location', ...
         'standardName', 'longitude', ...
         'units', 'degree_east', ...
         'fillValue', double(99999), ...
         'validMin', double(-180), ...
         'validMax', double(180), ...
         'axis', 'X', ...
         'paramType', '', ...
         'paramNcType', 'NC_DOUBLE', ...
         'adjAllowed', 0);

   case 'PRES'
      o_attributeStruct = struct('name', 'PRES', ...
         'longName', 'Sea water pressure, equals 0 at sea-level', ...
         'standardName', 'sea_water_pressure', ...
         'fillValue', single(99999), ...
         'units', 'decibar', ...
         'validMin', single(0), ...
         'validMax', single(12000), ...
         'axis', 'Z', ...
         'cFormat', '%7.1f', ...
         'fortranFormat', 'F7.1', ...
         'resolution', single(0.1), ...
         'paramType', 'c', ...
         'paramNcType', 'NC_FLOAT', ...
         'adjAllowed', 1);

   otherwise

      fprintf('ERROR: Attribute list no yet defined for parameter %s\n', a_paramName);

end

return

% ------------------------------------------------------------------------------
% Retrieve data from NetCDF file.
%
% SYNTAX :
%  [o_ncData] = get_data_from_nc_file(a_ncPathFileName, a_wantedVars)
%
% INPUT PARAMETERS :
%   a_ncPathFileName : NetCDF file name
%   a_wantedVars     : NetCDF variables to retrieve from the file
%
% OUTPUT PARAMETERS :
%   o_ncData : retrieved data
%
% EXAMPLES :
%
% SEE ALSO :
% AUTHORS  : Jean-Philippe Rannou (Altran)(jean-philippe.rannou@altran.com)
% ------------------------------------------------------------------------------
% RELEASES :
%   01/15/2014 - RNU - creation
% ------------------------------------------------------------------------------
function [o_ncData] = get_data_from_nc_file(a_ncPathFileName, a_wantedVars)

% output parameters initialization
o_ncData = [];


if (exist(a_ncPathFileName, 'file') == 2)

   % open NetCDF file
   fCdf = netcdf.open(a_ncPathFileName, 'NC_NOWRITE');
   if (isempty(fCdf))
      fprintf('ERROR: Unable to open NetCDF input file: %s\n', a_ncPathFileName);
      return
   end

   % retrieve variables from NetCDF file
   for idVar = 1:length(a_wantedVars)
      varName = a_wantedVars{idVar};

      if (var_is_present_dec_argo(fCdf, varName))
         varValue = netcdf.getVar(fCdf, netcdf.inqVarID(fCdf, varName));
         o_ncData = [o_ncData {varName} {varValue}];
      else
         %          fprintf('WARNING: Variable %s not present in file : %s\n', ...
         %             varName, a_ncPathFileName);
         o_ncData = [o_ncData {varName} {' '}];
      end

   end

   netcdf.close(fCdf);
end

return

% ------------------------------------------------------------------------------
% Check if a given variable is present in a NetCDF file.
%
% SYNTAX :
%  [o_present] = var_is_present_dec_argo(a_ncId, a_varName)
%
% INPUT PARAMETERS :
%   a_ncId    : NetCDF file Id
%   a_varName : variable name
%
% OUTPUT PARAMETERS :
%   o_present : 1 if the variable is present (0 otherwise)
%
% EXAMPLES :
%
% SEE ALSO :
% AUTHORS  : Jean-Philippe Rannou (Altran)(jean-philippe.rannou@altran.com)
% ------------------------------------------------------------------------------
% RELEASES :
%   05/27/2014 - RNU - creation
% ------------------------------------------------------------------------------
function [o_present] = var_is_present_dec_argo(a_ncId, a_varName)

o_present = 0;

[nbDims, nbVars, nbGAtts, unlimId] = netcdf.inq(a_ncId);

for idVar= 0:nbVars-1
   [varName, varType, varDims, nbAtts] = netcdf.inqVar(a_ncId, idVar);
   if (strcmp(varName, a_varName))
      o_present = 1;
      break
   end
end

return

% ------------------------------------------------------------------------------
% Interpolate between 2 dated locations.
%
% SYNTAX :
%  [o_interpLocLon, o_interpLocLat] = interpolate_between_2_locations(...
%    a_firstLocDate, a_firstLocLon, a_firstLocLat, ...
%    a_secondLocDate, a_secondLocLon, a_secondLocLat, ...
%    a_interpDate)
%
% INPUT PARAMETERS :
%   a_firstLocDate  : date of the first location
%   a_firstLocLon   : longitude of the first location
%   a_firstLocLat   : latitude of the first location
%   a_secondLocDate : date of the second location
%   a_secondLocLon  : longitude of the second location
%   a_secondLocLat  : latitude of the second location
%   a_interpDate    : date of the interpolation
%
% OUTPUT PARAMETERS :
%   o_interpLocLon : interpolated longitude
%   o_interpLocLat : interpolated latitude
%
% EXAMPLES :
%
% SEE ALSO :
% AUTHORS  : Jean-Philippe Rannou (Altran)(jean-philippe.rannou@altran.com)
% ------------------------------------------------------------------------------
% RELEASES :
%   05/18/2017 - RNU - creation
% ------------------------------------------------------------------------------
function [o_interpLocLon, o_interpLocLat] = interpolate_between_2_locations(...
   a_firstLocDate, a_firstLocLon, a_firstLocLat, ...
   a_secondLocDate, a_secondLocLon, a_secondLocLat, ...
   a_interpDate)

% output parameters initialization
o_interpLocLon = [];
o_interpLocLat = [];


% interpolate between the locations
if (((abs(a_firstLocLon) > 90) && (abs(a_secondLocLon) > 90)) && ...
      (((a_firstLocLon > 0) && (a_secondLocLon < 0)) || ((a_secondLocLon > 0) && (a_firstLocLon < 0))))
   % the float crossed the date line
   if (a_secondLocLon < 0)
      a_secondLocLon = a_secondLocLon + 360;
   else
      a_firstLocLon = a_firstLocLon + 360;
   end
   o_interpLocLon = interp1q([a_firstLocDate; a_secondLocDate], [a_firstLocLon; a_secondLocLon], a_interpDate);
   if (o_interpLocLon >= 180)
      o_interpLocLon = o_interpLocLon - 360;
   end
else
   o_interpLocLon = interp1q([a_firstLocDate; a_secondLocDate], [a_firstLocLon; a_secondLocLon], a_interpDate);
end
o_interpLocLat = interp1q([a_firstLocDate; a_secondLocDate], [a_firstLocLat; a_secondLocLat], a_interpDate);

return

% ------------------------------------------------------------------------------
% Compute geographic boundaries to plot a set of locations.
%
% SYNTAX :
%  [o_lonMin, o_lonMax, o_latMin, o_latMax] = ...
%    compute_geo_extrema(a_date, a_lon, a_lat, a_zoom)
%
% INPUT PARAMETERS :
%   a_date : date of the locations
%   a_lon  : longitude of the locations
%   a_lat  : latitude of the locations
%   a_zoom : zoom factor
%
% OUTPUT PARAMETERS :
%   o_lonMin  : min longitude of the plot
%   o_lonMax  : max longitude of the plot
%   o_latMin  : min latitude of the plot
%   o_latMax  : max latitude of the plot
%
% EXAMPLES :
%
% SEE ALSO :
% AUTHORS  : Jean-Philippe Rannou (Altran)(jean-philippe.rannou@altran.com)
% ------------------------------------------------------------------------------
% RELEASES :
%   08/01/2014 - RNU - creation
% ------------------------------------------------------------------------------
function [o_lonMin, o_lonMax, o_latMin, o_latMax] = ...
   compute_geo_extrema(a_date, a_lon, a_lat, a_zoom)

o_lonMin = [];
o_lonMax = [];
o_latMin = [];
o_latMax = [];

global g_dateDef;
global g_latDef;
global g_lonDef;

% default values initialization
init_valdef;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the geographic boundaries of the plot

if (~isempty(a_date))
   idNoData = find((a_date == g_dateDef) | (a_lon == g_lonDef) | (a_lat == g_latDef));
else
   idNoData = find((a_lon == g_lonDef) | (a_lat == g_latDef));
end
a_lon(idNoData) = [];
a_lat(idNoData) = [];

% geographic boundaries of the locations
latMin = min(a_lat);
latMax = max(a_lat);
latMarge = abs((latMax-latMin)/5);
if (latMarge == 0)
   latMarge = 1/60;
end
latMin = latMin - latMarge;
latMax = latMax + latMarge;

lonMin = min(a_lon);
lonMax = max(a_lon);

borneLonMax = 180;
if ((abs(lonMin - lonMax) > 180) && ...
      (abs(lonMin - lonMax) > abs(lonMin + lonMax)))
   id = find(a_lon < 0);
   a_lon(id) = a_lon(id) + 360;
   lonMin = min(a_lon);
   lonMax = max(a_lon);
   borneLonMax = 360;
end

lonMarge = abs((lonMax-lonMin)/5);
if (lonMarge == 0)
   lonMarge = 1/60;
end
lonMin = lonMin - lonMarge;
lonMax = lonMax + lonMarge;

% use of zoom factor
deltaLat = abs((latMax-latMin)/2);
deltaLon = abs((lonMax-lonMin)/2);
latMin = latMin - a_zoom*2*deltaLat;
latMax = latMax + a_zoom*2*deltaLat;
lonMin = lonMin - a_zoom*2*deltaLon;
lonMax = lonMax + a_zoom*2*deltaLon;

if (latMin < -90)
   latMin = -90;
end
if (latMax > 90)
   latMax = 90;
end
if (lonMin < -180)
   lonMin = -180;
end
if (lonMax > borneLonMax)
   lonMax = borneLonMax;
end

o_latMin = latMin;
o_latMax = latMax;
o_lonMin = lonMin;
o_lonMax = lonMax;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% optimization of the drawing window

m_proj('mercator', 'latitudes', [latMin latMax], 'longitudes', [lonMin lonMax]);

[xSW, ySW] = m_ll2xy(lonMin, latMin);
[xNW, yNW] = m_ll2xy(lonMin, latMax);
[xSE, ySE] = m_ll2xy(lonMax, latMin);

% use X/Y = 1.4 (with a legend, otherwise use 1.2)
coef = 1.4;
deltaX = xSE - xSW;
deltaY = yNW - ySW;
if (deltaX/deltaY > coef)
   complement = (deltaX/coef) - deltaY;
   ySW = ySW - complement/2;
   yNW = yNW + complement/2;
else
   complement = deltaY*coef - deltaX;
   xSW = xSW - complement/2;
   xSE = xSE + complement/2;
end

[lonMin, latMin] = m_xy2ll(xSW, ySW);
[bidon, latMax] = m_xy2ll(xNW, yNW);
[lonMax, bidon] = m_xy2ll(xSE, ySE);

if ((latMin >= -90) && (latMax <= 90) && (lonMin >= -180) && (lonMax <= borneLonMax))
   o_latMin = latMin;
   o_latMax = latMax;
   o_lonMin = lonMin;
   o_lonMax = lonMax;
else
   fprintf('compute_geo_extrema: cannot use the optimization of the drawing window\n');
end

return

% ------------------------------------------------------------------------------
% Retrieve the elevations of a given zone from the GEBCO 2019 file.
%
% SYNTAX :
%  [o_elev, o_lon, o_lat] = get_gebco_elev_zone( ...
%    a_lonMin, a_lonMax, a_latMin, a_latMax, a_gebcoFileName)
%
% INPUT PARAMETERS :
%   a_lonMin        : min longitude of the zone
%   a_lonMax        : max longitude of the zone
%   a_latMin        : min latitude of the zone
%   a_latMax        : max latitude of the zone
%   a_gebcoFileName : GEBCO 2019 file path name
%
% OUTPUT PARAMETERS :
%   o_elev : elevations of locations of the grid
%   o_lon  : longitudes of locations of the grid
%   o_lat  : latitudes of locations of the grid
%
% EXAMPLES :
%
% SEE ALSO :
% AUTHORS  : Jean-Philippe Rannou (Altran)(jean-philippe.rannou@altran.com)
% ------------------------------------------------------------------------------
% RELEASES :
%   04/27/2020 - RNU - creation
% ------------------------------------------------------------------------------
function [o_elev, o_lon, o_lat] = get_gebco_elev_zone( ...
   a_lonMin, a_lonMax, a_latMin, a_latMax, a_gebcoFileName)

% output parameters initialization
o_elev = [];
o_lon = [];
o_lat = [];

if (isempty(a_gebcoFileName))
   a_gebcoFileName = 'C:\Users\jprannou\_RNU\_ressources\GEBCO_2021\GEBCO_2021.nc';
end


% check inputs
if (a_latMin > a_latMax)
   fprintf('ERROR: get_gebco_elev_zone: latMin > latMax\n');
   return
else
   if (a_latMin < -90)
      fprintf('ERROR: get_gebco_elev_zone: latMin < -90\n');
      return
   elseif (a_latMax > 90)
      fprintf('ERROR: get_gebco_elev_zone: a_latMax > 90\n');
      return
   end
end
if (a_lonMin >= 180)
   a_lonMin = a_lonMin - 360;
   a_lonMax = a_lonMax - 360;
end
if (a_lonMax < a_lonMin)
   a_lonMax = a_lonMax + 360;
end

% check GEBCO file exists
if ~(exist(a_gebcoFileName, 'file') == 2)
   fprintf('ERROR: GEBCO file not found (%s)\n', a_gebcoFileName);
   return
end

% open NetCDF file
fCdf = netcdf.open(a_gebcoFileName, 'NC_NOWRITE');
if (isempty(fCdf))
   fprintf('RTQC_ERROR: Unable to open NetCDF input file: %s\n', a_gebcoFileName);
   return
end

lonVarId = netcdf.inqVarID(fCdf, 'lon');
latVarId = netcdf.inqVarID(fCdf, 'lat');
elevVarId = netcdf.inqVarID(fCdf, 'elevation');

lon = netcdf.getVar(fCdf, lonVarId);
lat = netcdf.getVar(fCdf, latVarId);
minLon = min(lon);
maxLon = max(lon);

idLigStart = find(lat <= a_latMin, 1, 'last');
idLigEnd = find(lat >= a_latMax, 1, 'first');
latVal = lat(fliplr(idLigStart:idLigEnd));

% a_lonMin is in the [-180, 180[ interval
% a_lonMax can be in the [-180, 180[ interval (case A) or [0, 360[ interval (case B)

% if ((a_lonMax - a_lonMin) > (maxLon - minLon)) we return the whole set of longitudes
% otherwise
% in case A: we should manage 3 zones
% [-180, minLon[, [minLon, maxLon] and ]maxLon, -180[, thus 5 cases
% case A1: a_lonMin and a_lonMax in [-180, minLon[
% case A2: a_lonMin in [-180, minLon[ and a_lonMax in [minLon, maxLon]
% case A3: a_lonMin in [minLon, maxLon] and a_lonMax in [minLon, maxLon]
% case A4: a_lonMin in [minLon, maxLon] and a_lonMax in ]maxLon, -180[
% case A5: a_lonMin in ]maxLon, -180[ and a_lonMax in ]maxLon, -180[
% in case B: we should manage 3 zones
% [minLon, maxLon], ]maxLon, -180[, [180, minLon+360[ and [minLon+360, maxLon+360], thus 4 cases
% case B1: a_lonMin in [minLon, maxLon] and a_lonMax in [180, minLon+360[
% case B2: a_lonMin in [minLon, maxLon] and a_lonMax in [minLon+360, maxLon+360]
% case B3: a_lonMin in ]maxLon, -180[ and a_lonMax in [180, minLon+360[
% case B4: a_lonMin in ]maxLon, -180[ and a_lonMax in [minLon+360, maxLon+360]

if ((a_lonMax - a_lonMin) <= (maxLon - minLon))
   if (a_lonMax < 180) % case A
      if ((a_lonMin >= minLon) && (a_lonMin <= maxLon) && ...
            (a_lonMax >= minLon) && (a_lonMax <= maxLon))
         % case A3
         idColStart = find(lon <= a_lonMin, 1, 'last');
         idColEnd = find(lon >= a_lonMax, 1, 'first');

         elev = nan(length(idLigStart:idLigEnd), length(idColStart:idColEnd));
         for idL = idLigStart:idLigEnd
            elev(end-(idL-idLigStart), :) = netcdf.getVar(fCdf, elevVarId, fliplr([idL-1 idColStart-1]), fliplr([1 length(idColStart:idColEnd)]))';
         end

         lonVal = lon(idColStart:idColEnd);
      elseif ((a_lonMin < minLon) && ...
            (a_lonMax >= minLon) && (a_lonMax <= maxLon))
         % case A2
         elev1 = nan(length(idLigStart:idLigEnd), 1);
         for idL = idLigStart:idLigEnd
            elev1(end-(idL-idLigStart), :) = netcdf.getVar(fCdf, elevVarId, fliplr([idL-1 length(lon)-1]), fliplr([1 1]))';
         end

         lonVal1 = lon(end);

         idColStart = 1;
         idColEnd = find(lon >= a_lonMax, 1, 'first');

         elev2 = nan(length(idLigStart:idLigEnd), length(idColStart:idColEnd));
         for idL = idLigStart:idLigEnd
            elev2(end-(idL-idLigStart), :) = netcdf.getVar(fCdf, elevVarId, fliplr([idL-1 idColStart-1]), fliplr([1 length(idColStart:idColEnd)]))';
         end

         lonVal2 = lon(idColStart:idColEnd) + 360;

         elev = cat(2, elev1, elev2);
         lonVal = cat(1, lonVal1, lonVal2);
         clear elev1 elev2 lonVal1 lonVal2
      elseif ((a_lonMin >= minLon) && (a_lonMin <= maxLon) && ...
            (a_lonMax > maxLon))
         % case A4
         idColStart = find(lon <= a_lonMin, 1, 'last');
         idColEnd = length(lon);

         elev1 = nan(length(idLigStart:idLigEnd), length(idColStart:idColEnd));
         for idL = idLigStart:idLigEnd
            elev1(end-(idL-idLigStart), :) = netcdf.getVar(fCdf, elevVarId, fliplr([idL-1 idColStart-1]), fliplr([1 length(idColStart:idColEnd)]))';
         end

         lonVal1 = lon(idColStart:idColEnd);

         elev2 = nan(length(idLigStart:idLigEnd), 1);
         for idL = idLigStart:idLigEnd
            elev2(end-(idL-idLigStart), :) = netcdf.getVar(fCdf, elevVarId, fliplr([idL-1 0]), fliplr([1 1]))';
         end

         lonVal2 = lon(1) + 360;

         elev = cat(2, elev1, elev2);
         lonVal = cat(1, lonVal1, lonVal2);
         clear elev1 elev2 lonVal1 lonVal2
      elseif ((a_lonMin < minLon) && ...
            (a_lonMax < minLon))
         % case A1
         elev1 = nan(length(idLigStart:idLigEnd), 1);
         for idL = idLigStart:idLigEnd
            elev1(end-(idL-idLigStart), :) = netcdf.getVar(fCdf, elevVarId, fliplr([idL-1 length(lon)-1]), fliplr([1 1]))';
         end

         lonVal1 = lon(end);

         elev2 = nan(length(idLigStart:idLigEnd), 1);
         for idL = idLigStart:idLigEnd
            elev2(end-(idL-idLigStart), :) = netcdf.getVar(fCdf, elevVarId, fliplr([idL-1 0]), fliplr([1 1]))';
         end

         lonVal2 = lon(1) + 360;

         elev = cat(2, elev1, elev2);
         lonVal = cat(1, lonVal1, lonVal2);
         clear elev1 elev2 lonVal1 lonVal2
      elseif ((a_lonMin > maxLon) && ...
            (a_lonMax > maxLon))
         % case A5
         elev1 = nan(length(idLigStart:idLigEnd), 1);
         for idL = idLigStart:idLigEnd
            elev1(end-(idL-idLigStart), :) = netcdf.getVar(fCdf, elevVarId, fliplr([idL-1 length(lon)-1]), fliplr([1 1]))';
         end

         lonVal1 = lon(end);

         elev2 = nan(length(idLigStart:idLigEnd), 1);
         for idL = idLigStart:idLigEnd
            elev2(end-(idL-idLigStart), :) = netcdf.getVar(fCdf, elevVarId, fliplr([idL-1 0]), fliplr([1 1]))';
         end

         lonVal2 = lon(1) + 360;

         elev = cat(2, elev1, elev2);
         lonVal = cat(1, lonVal1, lonVal2);
         clear elev1 elev2 lonVal1 lonVal2
      end
   else % case B
      if (a_lonMin <= maxLon) && (a_lonMax >= minLon + 360)
         % case B2
         idColStart = find(lon <= a_lonMin, 1, 'last');
         idColEnd = length(lon);

         elev1 = nan(length(idLigStart:idLigEnd), length(idColStart:idColEnd));
         for idL = idLigStart:idLigEnd
            elev1(end-(idL-idLigStart), :) = netcdf.getVar(fCdf, elevVarId, fliplr([idL-1 idColStart-1]), fliplr([1 length(idColStart:idColEnd)]))';
         end

         lonVal1 = lon(idColStart:idColEnd);

         idColStart = 1;
         idColEnd = find(lon >= a_lonMax - 360, 1, 'first');

         elev2 = nan(length(idLigStart:idLigEnd), length(idColStart:idColEnd));
         for idL = idLigStart:idLigEnd
            elev2(end-(idL-idLigStart), :) = netcdf.getVar(fCdf, elevVarId, fliplr([idL-1 idColStart-1]), fliplr([1 length(idColStart:idColEnd)]))';
         end

         lonVal2 = lon(idColStart:idColEnd) + 360;

         elev = cat(2, elev1, elev2);
         lonVal = cat(1, lonVal1, lonVal2);
         clear elev1 elev2 lonVal1 lonVal2
      elseif (a_lonMin <= maxLon) && (a_lonMax < minLon + 360)
         % case B1
         idColStart = find(lon <= a_lonMin, 1, 'last');
         idColEnd = length(lon);

         elev1 = nan(length(idLigStart:idLigEnd), length(idColStart:idColEnd));
         for idL = idLigStart:idLigEnd
            elev1(end-(idL-idLigStart), :) = netcdf.getVar(fCdf, elevVarId, fliplr([idL-1 idColStart-1]), fliplr([1 length(idColStart:idColEnd)]))';
         end

         lonVal1 = lon(idColStart:idColEnd);

         elev2 = nan(length(idLigStart:idLigEnd), 1);
         for idL = idLigStart:idLigEnd
            elev2(end-(idL-idLigStart), :) = netcdf.getVar(fCdf, elevVarId, fliplr([idL-1 0]), fliplr([1 1]))';
         end

         lonVal2 = lon(1) + 360;

         elev = cat(2, elev1, elev2);
         lonVal = cat(1, lonVal1, lonVal2);
         clear elev1 elev2 lonVal1 lonVal2
      elseif (a_lonMin > maxLon) && (a_lonMax >= minLon + 360)
         % case B4
         elev1 = nan(length(idLigStart:idLigEnd), 1);
         for idL = idLigStart:idLigEnd
            elev1(end-(idL-idLigStart), :) = netcdf.getVar(fCdf, elevVarId, fliplr([idL-1 length(lon)-1]), fliplr([1 1]))';
         end

         lonVal1 = lon(end);

         idColStart = 1;
         idColEnd = find(lon >= a_lonMax - 360, 1, 'first');

         elev2 = nan(length(idLigStart:idLigEnd), length(idColStart:idColEnd));
         for idL = idLigStart:idLigEnd
            elev2(end-(idL-idLigStart), :) = netcdf.getVar(fCdf, elevVarId, fliplr([idL-1 idColStart-1]), fliplr([1 length(idColStart:idColEnd)]))';
         end

         lonVal2 = lon(idColStart:idColEnd) + 360;

         elev = cat(2, elev1, elev2);
         lonVal = cat(1, lonVal1, lonVal2);
         clear elev1 elev2 lonVal1 lonVal2
      elseif (a_lonMin > maxLon) && (a_lonMax < minLon + 360)
         % case B3
         elev1 = nan(length(idLigStart:idLigEnd), 1);
         for idL = idLigStart:idLigEnd
            elev1(end-(idL-idLigStart), :) = netcdf.getVar(fCdf, elevVarId, fliplr([idL-1 length(lon)-1]), fliplr([1 1]))';
         end

         lonVal1 = lon(end);

         elev2 = nan(length(idLigStart:idLigEnd), 1);
         for idL = idLigStart:idLigEnd
            elev2(end-(idL-idLigStart), :) = netcdf.getVar(fCdf, elevVarId, fliplr([idL-1 0]), fliplr([1 1]))';
         end

         lonVal2 = lon(1) + 360;

         elev = cat(2, elev1, elev2);
         lonVal = cat(1, lonVal1, lonVal2);
         clear elev1 elev2 lonVal1 lonVal2
      end

   end
else % return the whole set of longitudes
   idColStart = 1;
   idColEnd = length(lon);

   elev = nan(length(idLigStart:idLigEnd), length(idColStart:idColEnd));
   for idL = idLigStart:idLigEnd
      elev(end-(idL-idLigStart), :) = netcdf.getVar(fCdf, elevVarId, fliplr([idL-1 idColStart-1]), fliplr([1 length(idColStart:idColEnd)]))';
   end

   lonVal = lon(idColStart:idColEnd);
end

netcdf.close(fCdf);

[longitudes, latitudes] = meshgrid(lonVal, latVal);

o_elev = elev;
o_lon = longitudes;
o_lat = latitudes;

clear lon lat elev longitudes latitudes

return



function [range,A12,A21]=distance_lpo(lat,long,argu1,argu2);
% DIST    Computes distance and bearing between points on the earth
%         using various reference spheroids.
%
%         [RANGE,AF,AR]=DIST(LAT,LONG) computes the ranges RANGE between
%         points specified in the LAT and LONG vectors (decimal degrees with
%         positive indicating north/east). Forward and reverse bearings
%         (degrees) are returned in AF, AR.
%
%         [RANGE,GLAT,GLONG]=DIST(LAT,LONG,N) computes N-point geodesics
%         between successive points. Each successive geodesic occupies
%         it's own row (N>=2)
%
%         [..]=DIST(...,'ellipsoid') uses the specified ellipsoid
%         to get distances and bearing. Available ellipsoids are:
%
%         'clarke66'  Clarke 1866
%         'iau73'     IAU 1973
%         'wgs84'     WGS 1984
%         'sphere'    Sphere of radius 6371.0 km
%
%          The default is 'wgs84'.
%
%          Ellipsoid formulas are recommended for distance d<2000 km,
%          but can be used for longer distances.

%Notes: RP (WHOI) 3/Dec/91
%         Mostly copied from BDC "dist.f" routine (copied from ....?), but
%         then wildly modified to bring it in line with Matlab vectorization.
%
%       RP (WHOI) 6/Dec/91
%         Feeping Creaturism! - added geodesic computations. This turned
%         out to be pretty hairy since there were a lot of branch problems
%         with asin, atan when computing geodesics subtending > 90 degrees
%         that were ignored in the original code!
%       RP (WHOI) 15/Jan/91
%         Fixed some bothersome special cases, like when computing geodesics
%         and N=2, or LAT=0...
%	A Newhall (WHOI) Sep 1997
%	   modified and fixed a bug found in Matlab version 5
%
%		NOTE: This routine may interfere with dist that
%			is supplied with matlab's neural net toolbox.

%C GIVEN THE LATITUDES AND LONGITUDES (IN DEG.) IT ASSUMES THE IAU SPHERO
%C DEFINED IN THE NOTES ON PAGE 523 OF THE EXPLANATORY SUPPLEMENT TO THE
%C AMERICAN EPHEMERIS.
%C
%C THIS PROGRAM COMPUTES THE DISTANCE ALONG THE NORMAL
%C SECTION (IN M.) OF A SPECIFIED REFERENCE SPHEROID GIVEN
%C THE GEODETIC LATITUDES AND LONGITUDES OF THE END POINTS
%C  *** IN DECIMAL DEGREES ***
%C
%C  IT USES ROBBIN'S FORMULA, AS GIVEN BY BOMFORD, GEODESY,
%C FOURTH EDITION, P. 122.  CORRECT TO ONE PART IN 10**8
%C AT 1600 KM.  ERRORS OF 20 M AT 5000 KM.
%C
%C   CHECK:  SMITHSONIAN METEOROLOGICAL TABLES, PP. 483 AND 484,
%C GIVES LENGTHS OF ONE DEGREE OF LATITUDE AND LONGITUDE
%C AS A FUNCTION OF LATITUDE. (SO DOES THE EPHEMERIS ABOVE)
%C
%C PETER WORCESTER, AS TOLD TO BRUCE CORNUELLE...1983 MAY 27
%C

spheroid='wgs84';
geodes=0;
if (nargin >= 3),
   if (isstr(argu1)),
      spheroid=argu1;
   else
      geodes=1;
      Ngeodes=argu1;
      if (Ngeodes <2), error('Must have at least 2 points in a goedesic!');end;
      if (nargin==4), spheroid=argu2; end;
   end;
end;

if (spheroid(1:3)=='sph'),
   A = 6371000.0;
   B = A;
   E = sqrt(A*A-B*B)/A;
   EPS= E*E/(1-E*E);
elseif (spheroid(1:3)=='cla'),
   A = 6378206.4E0;
   B = 6356583.8E0;
   E= sqrt(A*A-B*B)/A;
   EPS = E*E/(1.-E*E);
elseif(spheroid(1:3)=='iau'),
   A = 6378160.e0;
   B = 6356774.516E0;
   E = sqrt(A*A-B*B)/A;
   EPS = E*E/(1.-E*E);
elseif(spheroid(1:3)=='wgs'),

   %c on 9/11/88, Peter Worcester gave me the constants for the
   %c WGS84 spheroid, and he gave A (semi-major axis), F = (A-B)/A
   %c (flattening) (where B is the semi-minor axis), and E is the
   %c eccentricity, E = ( (A**2 - B**2)**.5 )/ A
   %c the numbers from peter are: A=6378137.; 1/F = 298.257223563
   %c E = 0.081819191
   A = 6378137.;
   E = 0.081819191;
   B = sqrt(A.^2 - (A*E).^2);
   EPS= E*E/(1.-E*E);

else
   error('dist: Unknown spheroid specified!');
end;


NN=max(size(lat));
if (NN ~= max(size(long))),
   error('dist: Lat, Long vectors of different sizes!');
end

if (NN==size(lat)), rowvec=0;  % It is easier if things are column vectors,
else                rowvec=1; end; % but we have to fix things before returning!

lat=lat(:)*pi/180;     % convert to radians
long=long(:)*pi/180;

lat(lat==0)=eps*ones(sum(lat==0),1);  % Fixes some nasty 0/0 cases in the
% geodesics stuff

PHI1=lat(1:NN-1);    % endpoints of each segment
XLAM1=long(1:NN-1);
PHI2=lat(2:NN);
XLAM2=long(2:NN);

% wiggle lines of constant lat to prevent numerical probs.
if (any(PHI1==PHI2)),
   for ii=1:NN-1,
      if (PHI1(ii)==PHI2(ii)), PHI2(ii)=PHI2(ii)+ 1e-14; end;
   end;
end;
% wiggle lines of constant long to prevent numerical probs.
if (any(XLAM1==XLAM2)),
   for ii=1:NN-1,
      if (XLAM1(ii)==XLAM2(ii)), XLAM2(ii)=XLAM2(ii)+ 1e-14; end;
   end;
end;



%C  COMPUTE THE RADIUS OF CURVATURE IN THE PRIME VERTICAL FOR
%C EACH POINT

xnu=A./sqrt(1.0-(E*sin(lat)).^2);
xnu1=xnu(1:NN-1);
xnu2=xnu(2:NN);

%C*** COMPUTE THE AZIMUTHS.  A12 (A21) IS THE AZIMUTH AT POINT 1 (2)
%C OF THE NORMAL SECTION CONTAINING THE POINT 2 (1)

TPSI2=(1.-E*E)*tan(PHI2) + E*E*xnu1.*sin(PHI1)./(xnu2.*cos(PHI2));
PSI2=atan(TPSI2);

%C*** SOME FORM OF ANGLE DIFFERENCE COMPUTED HERE??

DPHI2=PHI2-PSI2;
DLAM=XLAM2-XLAM1;
CTA12=(cos(PHI1).*TPSI2 - sin(PHI1).*cos(DLAM))./sin(DLAM);
A12=atan((1.)./CTA12);
CTA21P=(sin(PSI2).*cos(DLAM) - cos(PSI2).*tan(PHI1))./sin(DLAM);
A21P=atan((1.)./CTA21P);

%C    GET THE QUADRANT RIGHT
DLAM2=(abs(DLAM)<pi).*DLAM + (DLAM>=pi).*(-2*pi+DLAM) + ...
   (DLAM<=-pi).*(2*pi+DLAM);
A12=A12+(A12<-pi)*2*pi-(A12>=pi)*2*pi;
A12=A12+pi*sign(-A12).*( sign(A12) ~= sign(DLAM2) );
A21P=A21P+(A21P<-pi)*2*pi-(A21P>=pi)*2*pi;
A21P=A21P+pi*sign(-A21P).*( sign(A21P) ~= sign(-DLAM2) );
%%A12*180/pi
%%A21P*180/pi


SSIG=sin(DLAM).*cos(PSI2)./sin(A12);
% At this point we are OK if the angle < 90...but otherwise
% we get the wrong branch of asin!
% This fudge will correct every case on a sphere, and *almost*
% every case on an ellipsoid (wrong hnadling will be when
% angle is almost exactly 90 degrees)
dd2=[cos(long).*cos(lat) sin(long).*cos(lat) sin(lat)];
dd2=sum((diff(dd2).*diff(dd2))')';
if ( any(abs(dd2-2) < 2*((B-A)/A))^2 ),
   disp('dist: Warning...point(s) too close to 90 degrees apart');
end;
bigbrnch=dd2>2;

SIG=asin(SSIG).*(bigbrnch==0) + (pi-asin(SSIG)).*bigbrnch;

SSIGC=-sin(DLAM).*cos(PHI1)./sin(A21P);
SIGC=asin(SSIGC);
A21 = A21P - DPHI2.*sin(A21P).*tan(SIG/2.0);

%C   COMPUTE RANGE

G2=EPS*(sin(PHI1)).^2;
G=sqrt(G2);
H2=EPS*(cos(PHI1).*cos(A12)).^2;
H=sqrt(H2);
TERM1=-SIG.*SIG.*H2.*(1.0-H2)/6.0;
TERM2=(SIG.^3).*G.*H.*(1.0-2.0*H2)/8.0;
TERM3=(SIG.^4).*(H2.*(4.0-7.0*H2)-3.0*G2.*(1.0-7.0*H2))/120.0;
TERM4=-(SIG.^5).*G.*H/48.0;

range=xnu1.*SIG.*(1.0+TERM1+TERM2+TERM3+TERM4);


if (geodes),

   %c now calculate the locations along the ray path. (for extra accuracy, could
   %c do it from start to halfway, then from end for the rest, switching from A12
   %c to A21...
   %c started to use Rudoe's formula, page 117 in Bomford...(1980, fourth edition)
   %c but then went to Clarke's best formula (pg 118)

   %RP I am doing this twice because this formula doesn't work when we go
   %past 90 degrees!
   Ngd1=round(Ngeodes/2);

   % First time...away from point 1
   if (Ngd1>1),
      wns=ones(1,Ngd1);
      CP1CA12 = (cos(PHI1).*cos(A12)).^2;
      R2PRM = -EPS.*CP1CA12;
      R3PRM = 3.0*EPS.*(1.0-R2PRM).*cos(PHI1).*sin(PHI1).*cos(A12);
      C1 = R2PRM.*(1.0+R2PRM)/6.0*wns;
      C2 = R3PRM.*(1.0+3.0*R2PRM)/24.0*wns;
      R2PRM=R2PRM*wns;
      R3PRM=R3PRM*wns;

      %c  now have to loop over positions
      RLRAT = (range./xnu1)*([0:Ngd1-1]/(Ngeodes-1));

      THETA = RLRAT.*(1 - (RLRAT.^2).*(C1 - C2.*RLRAT));
      C3 = 1.0 - (R2PRM.*(THETA.^2))/2.0 - (R3PRM.*(THETA.^3))/6.0;
      DSINPSI =(sin(PHI1)*wns).*cos(THETA) + ...
         ((cos(PHI1).*cos(A12))*wns).*sin(THETA);
      %try to identify the branch...got to other branch if range> 1/4 circle
      PSI = asin(DSINPSI);

      DCOSPSI = cos(PSI);
      DSINDLA = (sin(A12)*wns).*sin(THETA)./DCOSPSI;
      DTANPHI=(1.0+EPS)*(1.0 - (E^2)*C3.*(sin(PHI1)*wns)./DSINPSI).*tan(PSI);
      %C compute output latitude (phi) and long (xla) in radians
      %c I believe these are absolute, and don't need source coords added
      PHI = atan(DTANPHI);
      %  fix branch cut stuff -
      otherbrcnh= sign(DLAM2*wns) ~= sign([sign(DLAM2) diff(DSINDLA')'] );
      XLA = XLAM1*wns + asin(DSINDLA).*(otherbrcnh==0) + ...
         (pi-asin(DSINDLA)).*(otherbrcnh);
   else
      PHI=PHI1;
      XLA=XLAM1;
   end;

   % Now we do the same thing, but in the reverse direction from the receiver!
   if (Ngeodes-Ngd1>1),
      wns=ones(1,Ngeodes-Ngd1);
      CP2CA21 = (cos(PHI2).*cos(A21)).^2;
      R2PRM = -EPS.*CP2CA21;
      R3PRM = 3.0*EPS.*(1.0-R2PRM).*cos(PHI2).*sin(PHI2).*cos(A21);
      C1 = R2PRM.*(1.0+R2PRM)/6.0*wns;
      C2 = R3PRM.*(1.0+3.0*R2PRM)/24.0*wns;
      R2PRM=R2PRM*wns;
      R3PRM=R3PRM*wns;

      %c  now have to loop over positions
      RLRAT = (range./xnu2)*([0:Ngeodes-Ngd1-1]/(Ngeodes-1));

      THETA = RLRAT.*(1 - (RLRAT.^2).*(C1 - C2.*RLRAT));
      C3 = 1.0 - (R2PRM.*(THETA.^2))/2.0 - (R3PRM.*(THETA.^3))/6.0;
      DSINPSI =(sin(PHI2)*wns).*cos(THETA) + ...
         ((cos(PHI2).*cos(A21))*wns).*sin(THETA);
      %try to identify the branch...got to other branch if range> 1/4 circle
      PSI = asin(DSINPSI);

      DCOSPSI = cos(PSI);
      DSINDLA = (sin(A21)*wns).*sin(THETA)./DCOSPSI;
      DTANPHI=(1.0+EPS)*(1.0 - (E^2)*C3.*(sin(PHI2)*wns)./DSINPSI).*tan(PSI);
      %C compute output latitude (phi) and long (xla) in radians
      %c I believe these are absolute, and don't need source coords added
      PHI = [PHI fliplr(atan(DTANPHI))];
      % fix branch cut stuff
      otherbrcnh= sign(-DLAM2*wns) ~= sign( [sign(-DLAM2) diff(DSINDLA')'] );
      XLA = [XLA fliplr(XLAM2*wns + asin(DSINDLA).*(otherbrcnh==0) + ...
         (pi-asin(DSINDLA)).*(otherbrcnh))];
   else
      PHI = [PHI PHI2];
      XLA = [XLA XLAM2];
   end;

   %c convert to degrees
   A12 = PHI*180/pi;
   A21 = XLA*180/pi;
   range=range*([0:Ngeodes-1]/(Ngeodes-1));


else

   %C*** CONVERT TO DECIMAL DEGREES
   A12=A12*180/pi;
   A21=A21*180/pi;
   if (rowvec),
      range=range';
      A12=A12';
      A21=A21';
   end;
end;

return