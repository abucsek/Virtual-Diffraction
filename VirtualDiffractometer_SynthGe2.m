% Created by Ashley Bucsek, Colorado School of Mines, 2016
% Produces synthetic diffraction patterns and .ge2 files
% Note: Ideal detector parameters only (i.e., only energy and sample-to-
% detector distance--no tilt, etc. )
% Note: Does not consider intensity parameters (e.g., structure factor)

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
GE2_PROPS                    = struct();                                   % GE2 FILE PROPERTIES (define below)
GE2_PROPS.deltaOmega         = 90;                                       % Omega step for GE2 files
GE2_PROPS.outPath            = '/Users/abucsek/Documents/CHESSDec15/Analysis/50NiTiSC_0D_1/fitgrains/HEXRD';
GE2_PROPS.outNameStem        = 'Synthetic';
%
% ***** SAMPLE_PROPS entries must have same number of rows *****
% y is along the loading axis (positive is up)
% z is in the beam direction (positive is toward the source)
SAMPLE_PROPS                 = struct();                                   % SAMPLE PROPERTIES (define below)
SAMPLE_PROPS.grainDimensions = [1.0 1.0 1.0;  1.0 1.0 1.0];                % Sample grain dimesnions (x,y,z)
SAMPLE_PROPS.grainCOMs       = [0.0 0.5 0.0;  0.0 -0.5 0.0];               % Sample COM locations (x,y,z)
SAMPLE_PROPS.meshSize        = 5;                                          % Sample mesh size, only increase if spots are scarcely filled
SAMPLE_PROPS.orientations    = [7.356e-1 6.616e-1 1.455e-1 -8.024e-3;  7.390e-1	6.623e-1 1.168e-1 3.913e-2]; % Sample orientations (quaternions)
SAMPLE_PROPS_strains         = {[1 0 0; 0 1 0; 0 0 1],  [1 0 0; 0 1 0; 0 0 1]};                              % Sample strains
SAMPLE_PROPS.intensities     = [100;  200];                                % Set intensities for all spots of each grain (nice way to differentiate grains by sight)


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


fprintf(['Processing ' num2str(size(SAMPLE_PROPS.grainDimensions,1)) ' grains \n \n']);

parLights = cell(size(SAMPLE_PROPS.grainDimensions,1),1);
for grain = 1 : size(SAMPLE_PROPS.grainDimensions,1)
    %% Create mosaic geometry - Parallelopiped
    SpecimenSize.x = SAMPLE_PROPS.grainDimensions(grain,1) * 1e-3;         % (m)
    SpecimenSize.y = SAMPLE_PROPS.grainDimensions(grain,1) * 1e-3;         % along vertical (loading axis)
    SpecimenSize.z = SAMPLE_PROPS.grainDimensions(grain,1) * 1e-3;         % in beam direction
    
    NumBlocks.x = SAMPLE_PROPS.meshSize;  Blocksize.x = SpecimenSize.x / NumBlocks.x;
    NumBlocks.y = NumBlocks.x;  Blocksize.y = SpecimenSize.y / NumBlocks.y;
    NumBlocks.z = NumBlocks.x;  Blocksize.z = SpecimenSize.z / NumBlocks.z;
    
    Block.x = linspace(-SpecimenSize.x/2+Blocksize.x/2,  SpecimenSize.x/2-Blocksize.x/2,  NumBlocks.x)';
    Block.y = linspace(-SpecimenSize.y/2+Blocksize.y/2,  SpecimenSize.y/2-Blocksize.y/2,  NumBlocks.y)';
    Block.z = linspace(-SpecimenSize.z/2+Blocksize.z/2,  SpecimenSize.z/2-Blocksize.z/2,  NumBlocks.z)';
    
    Volume = (SpecimenSize.x / NumBlocks.x) ^ 3;
    
    Mosaicity = zeros(NumBlocks.x*NumBlocks.y*NumBlocks.z, 4);  gg = 1;
    for ix = 1 : NumBlocks.x
        for iy = 1 : NumBlocks.y
            for iz = 1 : NumBlocks.z
                Mosaicity(gg,:) = [Block.x(ix) Block.y(iy) Block.z(iz) Volume];
                gg = gg + 1;
            end
        end
    end
    
    
    %% Calculate virtual detector data
    Orientation = quat2rot(SAMPLE_PROPS.orientations(grain,:));
    
    tVec_C = SAMPLE_PROPS.grainCOMs(grain,:) * 1e-3;
    
    Vinv_S = SAMPLE_PROPS_strains{grain};
    
    Lights = [];
    for jj = 1:size(hkl_list,1)
        % [omega ttheta zeta Int]
        Lights_temp = MosaicLightUp(hkl_list(jj,:), Wavelength, Distance, Orientation, Mosaicity, MATERIAL_PROPS, tVec_C, Vinv_S);
        Lights_temp(:,5) = SAMPLE_PROPS.intensities(grain);
        Lights = vertcat(Lights, Lights_temp);
    end
    
    parLights(grain,1) = {Lights};
end

fprintf('Light Up Finished.\n\n')

Lights = [];
Output = cell(size(SAMPLE_PROPS.grainCOMs,1), 1);
for ii = 1 : size(SAMPLE_PROPS.grainCOMs,1)
    Lights = parLights{ii,1};
    zeta_x = Lights(:,3);
    zeta_y = Lights(:,4);
    omega = Lights(:,1);
    OutputTemp = [zeta_x  zeta_y  omega*180/pi];
    ix = find(zeta_x <= 0.2048 & zeta_x >= -0.2048);
    OutputTemp = OutputTemp(ix,:);
    iy = find(OutputTemp(:,2) <= 0.2048 & OutputTemp(:,2) >= -0.2048);
    OutputTemp = OutputTemp(iy,:);
    [Temp, ia] = unique(OutputTemp(:,3));
    OutputTemp = OutputTemp(ia, :);
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

omega = [];  ttheta = [];  zeta_x = [];  zeta_y = [];  Int = [];
for grain = 1 : grain
    Lights_temp = parLights{grain};
    omega = vertcat(omega, Lights_temp(:,1));
    ttheta = vertcat(ttheta, Lights_temp(:,2));
    zeta_x = vertcat(zeta_x, Lights_temp(:,3));
    zeta_y = vertcat(zeta_y, Lights_temp(:,4));
    Int = vertcat(Int, Lights_temp(:,5));
end


%% Figure
binsize = 200e-6;
xbins = -0.2048+binsize/2 : binsize : 0.2048-binsize/2;  ybins = xbins;
[nx, idx] = histc(zeta_x, xbins);
[ny, idy] = histc(zeta_y, ybins);  idx = idx+1;  idy = idy+1;
out = zeros(length(xbins));
for jj = 1 : length(idx)
    out(idy(jj),idx(jj)) = max(out(idy(jj),idx(jj)), Int(jj)); % max
end

figure;  ax=axes;
imagesc(xbins, ybins, out, 'Parent', ax); hold on
set(ax, 'YDir', 'normal')
axis square
axis([-0.2048 0.2048 -0.2048 0.2048])
colormap bone;  caxis([0 100])


%% Sort into frames and export tiffs
clearvars -except GE2_PROPS omega zeta_x zeta_y Int Output xbins ybins
close all
omega = omega * 180/pi;
for pp = 1 : length(omega)
    if omega(pp) < 0
        omega(pp) = omega(pp) + 360;
    elseif omega(pp) >= 360
        omega(pp) = omega(pp) - 360;
    end
end

omeFrame = 0 : GE2_PROPS.deltaOmega : 360;

numFrames = 360 / GE2_PROPS.deltaOmega;
if numFrames/240 == round(numFrames/240)
    numFiles = numFrames / 240;
else
    numFiles = floor(numFrames / 240) + 1;
end

fprintf(['Writing ' num2str(numFiles) ' ge2 files \n \n'])

GE2_PROPS.outName = [GE2_PROPS.outNameStem '_01' '.ge2'];
fo = fopen(fullfile(GE2_PROPS.outPath,GE2_PROPS.outName),'w');

for k = 1 : numFiles
    
    if k == numFiles

        for i = 240*(k-1)+1 : numFrames
            deltaOmeIndexes = find(omega(:,1) >= omeFrame(i)  &  omega(:,1) < omeFrame(i+1));
            
            zeta_xFrame = zeta_x(deltaOmeIndexes);
            zeta_yFrame = zeta_y(deltaOmeIndexes);
            IntFrame = Int(deltaOmeIndexes);
            
            [nx, idx] = histc(zeta_xFrame, xbins);
            [ny, idy] = histc(zeta_yFrame, ybins);
            
            out = zeros(length(xbins));
            for ii = 1 : length(idx)
                if idx(ii)>0 && idy(ii)>0
                    out(idx(ii),idy(ii)) = max(out(idx(ii),idy(ii)), IntFrame(ii));
                end
            end
            
            if i == (k-1)*240+1
                GE2_PROPS.outName = [GE2_PROPS.outNameStem '_0' num2str(k) '.ge2'];  display(GE2_PROPS.outName)
                fo = fopen(fullfile(GE2_PROPS.outPath,GE2_PROPS.outName),'w');
                header = zeros(1,8192);
                fwrite(fo,header,'uint8');
                fseek(fo,size(header,2),'bof');
                fclose(fo);
            end
            
            fo = fopen(fullfile(GE2_PROPS.outPath,GE2_PROPS.outName),'a+');
            fwrite(fo, out, 'uint16');
            fclose(fo);
        end
        
    else
        
        for i = 240*(k-1)+1 : k*240
            deltaOmeIndexes = find(omega(:,1) >= omeFrame(i)  &  omega(:,1) < omeFrame(i+1));
            
            zeta_xFrame = zeta_x(deltaOmeIndexes);
            zeta_yFrame = zeta_y(deltaOmeIndexes);
            IntFrame = Int(deltaOmeIndexes);
            
            [nx, idx] = histc(zeta_xFrame, xbins);
            [ny, idy] = histc(zeta_yFrame, ybins);
            
            out = zeros(length(xbins));
            for ii = 1 : length(idx)
                if idx(ii)>0 && idy(ii)>0
                    out(idx(ii),idy(ii)) = max(out(idx(ii),idy(ii)), IntFrame(ii));
                end
            end
            
            if i == (k-1)*240+1
                GE2_PROPS.outName = [GE2_PROPS.outNameStem '_0' num2str(k) '.ge2'];  display(GE2_PROPS.outName)
                fo = fopen(fullfile(GE2_PROPS.outPath,GE2_PROPS.outName),'w');
                header = zeros(1,8192);
                fwrite(fo,header,'uint8');
                fseek(fo,size(header,2),'bof');
                fclose(fo);
            end
            
            fo = fopen(fullfile(GE2_PROPS.outPath,GE2_PROPS.outName),'a+');
            fwrite(fo, out, 'uint16');
            fclose(fo);
        end
        
    end
    
end

fprintf('Ge2 file creation complete \n \n')