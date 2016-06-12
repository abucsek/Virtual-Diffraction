% Created by Ashley Bucsek, Colorado School of Mines, 2016
% Compares heXRD output files against original ff-HEDM data (max-over-all)
% Note: Does not consider detector parameters (only energy and sample-to-
% detector distance)
% Note: Does not consider intensity parameters (e.g., structure factor)
% Note: heXRD output files (grains.out and accepted_orientations.dat) must 
% be in the same directory
%
% The 'Output' variable is [x y omega], where x,y are in pixels and omega
% is in degrees. Each cell corresponds to a different grain.

clear; clc


%% Inputs
% ***** distances in mm, angles in degrees *****
%
MATERIAL_PROPS               = struct();                                   % MATERIAL PROPERTIES (define below)
MATERIAL_PROPS.latticeParams = [2.9 2.9 2.9 90 90 90];                     % [a b c alpha beta gamma]
%
DETECTOR_PROPS               = struct();                                   % DETECTOR PARAMETERS (define below)
DETECTOR_PROPS.distance      = 1012.36;                                    % Detector to sample distance
DETECTOR_PROPS.beamEnergy    = 55.618;                                     % (keV)
%
DIRECTORIES                  = struct();                                   % FILE DIRECTORIES (define below)
DIRECTORIES.heXRD = '/Users/abucsek/Documents/CHESSDec15/Analysis/50NiTiSC_0D_1/fitgrains/HEXRD/CM';             % heXRD output files (accepted_orientations and grains) directory
DIRECTORIES.maxPattern_path = '/Users/abucsek/Documents/CHESSDec15/Analysis/50NiTiSC_0D_1/fitgrains/HEXRD';   % Raw data max-over-all ge2 file path
DIRECTORIES.maxPattern_name = 'HexrdTesting_CM-max_img.ge2';               % Raw data max-over-all ge2 file
%
RING_HKLS                   = [1 0 0; 1 1 0; 1 1 1; 0 2 0; 1 2 0];         % HKLS FOR RINGS (if desired)


%% Set up
% Distance from sample to detector (m)
Distance = DETECTOR_PROPS.distance * 1e-3;

% Incoming beam wavelength (A)
beamEnergy = DETECTOR_PROPS.beamEnergy;        % (keV)
h = 6.626E-34;  c = 3.000E+08;  e = 1.602E-19;
Wavelength = h * c / ( 1000 * beamEnergy * e ) * 1e10; % (A)

% hkl families of interest
initial_hkl_list = [1 0 0;  
    1 1 0;
    1 1 1;
    0 2 0;
    1 2 0;
    1 2 1;
    2 2 0;
    2 2 1;
    0 3 0;
    1 3 0;
    1 3 1;
    2 2 2;
    2 3 0;
    2 3 1;
    0 4 0;
    2 3 2;
    1 4 0;
    1 4 1;
    3 3 0;
    3 3 0;
    2 4 0;
    2 4 1;
    3 3 2;
    2 4 2];

temp = [];
for jj = 1:size(initial_hkl_list, 1)
    vecin = cubic_symmetries( transpose( initial_hkl_list(jj,:) ) );
    temp = vertcat(temp, vecin);
end
hkl_list=unique(temp,'rows');


%% Orientations, COM's, and Strains
Orientations = importdata(fullfile(DIRECTORIES.heXRD, 'accepted_orientations.dat'));

GrainsFile = fopen(fullfile(DIRECTORIES.heXRD, 'grains.out'));
GrainsData = textscan(GrainsFile, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'headerlines', 1);
GrainIDs = GrainsData{1};

tVec_Cs(:,1) = GrainsData{7};
tVec_Cs(:,2) = GrainsData{8};
tVec_Cs(:,3) = GrainsData{9};

StrainComponents(:,1) = GrainsData{10};
StrainComponents(:,2) = GrainsData{11};
StrainComponents(:,3) = GrainsData{12};
StrainComponents(:,4) = GrainsData{13};
StrainComponents(:,5) = GrainsData{14};
StrainComponents(:,6) = GrainsData{15};


%% Virtual detector data
fprintf(['Processing ' num2str(size(Orientations,1)) ' grains']);
fprintf('\n \n')

parLights = cell(size(Orientations,1),1);
for grainNum = 1 : size(Orientations,1)
    
    Orientation = quat2rot(Orientations(grainNum, :));
    
    tVec_C = tVec_Cs(grainNum,:) * 1e-3;
    
    Strains = StrainComponents(grainNum,:);
    Vinv_S = [Strains(1) Strains(6) Strains(5);
        Strains(6) Strains(2) Strains(4);
        Strains(5) Strains(4) Strains(3)];
    
    Lights = [];
    for jj = 1:size(hkl_list,1)        
        % [omega ttheta zeta Int]
        Lights_temp = LightUp(hkl_list(jj,:), Wavelength, Distance, Orientation, MATERIAL_PROPS, tVec_C, Vinv_S);
        Lights = vertcat(Lights, Lights_temp);
    end
    
    parLights(grainNum,1) = {Lights};
end

fprintf('\n \n')
fprintf('Light Up Finished.\n\n')

Output = cell(size(Orientations,1), 1);
for ii = 1 : size(Orientations,1)
    Lights = parLights{ii,1};
    zeta_x = Lights(:,3);
    zeta_y = Lights(:,4);
    omega = Lights(:,1);
    OutputTemp = [zeta_x  zeta_y  omega*180/pi];
    ix = find(zeta_x <= 0.2048 & zeta_x >= -0.2048);
    OutputTemp = OutputTemp(ix,:);
    iy = find(OutputTemp(:,2) <= 0.2048 & OutputTemp(:,2) >= -0.2048);
    OutputTemp = OutputTemp(iy,:);
    for pp = 1 : length(OutputTemp)
        if OutputTemp(pp,3) < 0
            OutputTemp(pp,3) = OutputTemp(pp,3) + 360;
        elseif OutputTemp(pp,3) >= 360
            OutputTemp(pp,3) = OutputTemp(pp,3) - 360;
        end
    end
    [temp,ind] = sort(OutputTemp(:,3));
    OutputTemp = OutputTemp(ind,:);
    OutputTemp(:,1:2) = OutputTemp(:,1:2) / 200e-6 + 1024;
    Output{ii,1} = OutputTemp;
end


%% Figure
OrientationsStream = fopen(fullfile(DIRECTORIES.maxPattern_path,DIRECTORIES.maxPattern_name), 'r', 'n');
fseek(OrientationsStream, 8192, 'bof');
dataReal = fread(OrientationsStream, [2048 2048], '*uint16');
dataRealRot = rot90(dataReal,1);

binsize = 200e-6;       % 200 um
xbins = -0.2048+binsize/2 : binsize : 0.2048-binsize/2;  ybins = xbins;

figure;  ax=axes;
imagesc(xbins, ybins, (dataRealRot), 'Parent', ax);  hold on
for ii = 1 : size(parLights,1)
    Lights = parLights{ii,1};
    plot(Lights(:,3), Lights(:,4), '.', 'markersize', 7);  hold on
end
set(ax, 'YDir', 'normal');  axis square
colormap bone;  caxis([0 150])
xlabel('Detector x');  ylabel('Detector y')


%% Rings
a = MATERIAL_PROPS.latticeParams(1);
b = MATERIAL_PROPS.latticeParams(2);
c = MATERIAL_PROPS.latticeParams(3);
alpha = MATERIAL_PROPS.latticeParams(4) * pi/180.0;
beta = MATERIAL_PROPS.latticeParams(5) * pi/180.0;
gamma = MATERIAL_PROPS.latticeParams(6) * pi/180.0;
for i = 1 : size(RING_HKLS, 1)
    hkl = RING_HKLS(i,:);  h = hkl(1);  k = hkl(2);  l = hkl(3);
    
    dhkl = 1 / sqrt( h^2/(a^2 * sin(beta)^2) + k^2/b^2 + l^2/(c^2*sin(beta)^2) - 2 * h * l * cos(beta)/(a * c * sin(beta)^2) );
    RingTheta = asin(Wavelength / (2 * dhkl));  RingTTheta = 2 * RingTheta;
    
    RingRadius = Distance * tan(RingTTheta);
    RingAngles = 0 : 0.05 : 2*pi;
    RingX = RingRadius * cos(RingAngles);
    RingY = RingRadius * sin(RingAngles);
    figure(gcf)
    plot(RingX, RingY, 'r'); hold on
    text(RingRadius, 0.0, [num2str(h) num2str(k) num2str(l)], 'Rotation', 270, 'Color', 'r'); hold on
end