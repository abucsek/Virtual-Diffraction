% Created by Ashley Bucsek, Colorado School of Mines, 2016
% For use with 'VirtualDiffractometer_heXRD'
% 'VirtualDiffractometer_heXRD' heXRD output files against original ff-HEDM
% data (max-over-all)
% This file calculates the virtual diffraction data for ^


function Lights = LightUp(HKL, WAVELENGTH, DISTANCE, R_C2S, MATERIAL_PROPS, tVec_C, Vinv_S)
chi = 0;

%% Initiate
omega = [];
zeta = [];
Lights = [];     % [omega ttheta zeta Int]

%% Lab frame
Lx = [1; 0; 0];
Ly = [0; 1; 0];
Lz = [0; 0; 1];

ki_L = 2*pi / WAVELENGTH * (-Lz);                                           % Incoming beam wave vector in the lab frame


%% Crystal frame
% lattice params
aMag = MATERIAL_PROPS.latticeParams(1);
bMag = MATERIAL_PROPS.latticeParams(2);
cMag = MATERIAL_PROPS.latticeParams(3);
alpha = MATERIAL_PROPS.latticeParams(4) * pi/180.0;
beta = MATERIAL_PROPS.latticeParams(5) * pi/180.0;
gamma = MATERIAL_PROPS.latticeParams(6) * pi/180.0;

% lattice vectors
a = aMag * [1; 0; 0];
b = bMag * [0; 1; 0];
c = cMag * [cos(beta); 0; sin(beta)];
latticeVol = det([a b c]);

% reciprocal lattice parameters
alphastar = acos( (cos(beta) * cos(gamma) - cos(alpha)) / (sin(beta) * sin(gamma)) );
betastar = acos( (cos(alpha) * cos(gamma) - cos(beta)) / (sin(alpha) * sin(gamma)) );
gammastar = acos( (cos(alpha) * cos(beta) - cos(beta)) / (sin(alpha) * sin(beta)) );

% Change-of-basis matrix, B, that takes vector components in the reciprocal
% lattice frame to the crystal frame
B = 1/latticeVol * [bMag*cMag*sin(alphastar)*sin(beta)*sin(gamma)  0  0;
    -bMag*cMag*sin(alphastar)*sin(beta)*cos(gamma)  aMag*cMag*sin(alphastar)*sin(beta)  0;
    -bMag*cMag*(cos(alphastar)*sin(beta)*cos(gamma)+cos(beta)*sin(gamma))  aMag*cMag*cos(alphastar)*sin(beta)  aMag*bMag*sin(gamma)];


%% hkl
Gstar = [HKL(1); HKL(2); HKL(3)];
G_C = Vinv_S * R_C2S * B * Gstar;                                          % reciprocal lattice vector in the xtal frame
G_S = Vinv_S * R_C2S * B * Gstar;                                          % reciprocal lattice vector in the sample frame


%% Calculate the set of oscillation angles
% bhat is the unit vector aligned with the beam propagation direction
bhat_L = [0; 0; -1];

% Ghat_S the reciprocal lattice vector direction
Ghat_S = unit(G_S);

% SOLVE FOR OMEGAS
theta = -asin(-WAVELENGTH/2 * norm(G_C));  ttheta = 2 * theta;
aMag = Ghat_S(3)*bhat_L(1) + sin(chi)*Ghat_S(1)*bhat_L(2) - cos(chi)*Ghat_S(1)*bhat_L(3);
bMag = Ghat_S(1)*bhat_L(1) - sin(chi)*Ghat_S(3)*bhat_L(2) + cos(chi)*Ghat_S(3)*bhat_L(3);
cMag =     -sin(theta)     - cos(chi)*Ghat_S(2)*bhat_L(2) - sin(chi)*Ghat_S(2)*bhat_L(3);
alpha = atan2(bMag,aMag);
if abs(cMag/sqrt(aMag^2+bMag^2)) <= 1
    omega(1,1) = asin(cMag/sqrt(aMag^2+bMag^2)) - alpha;
    omega(2,1) = pi - asin(cMag/sqrt(aMag^2+bMag^2)) - alpha;
end


%% Calculate the set of unit diffraction vector components in the lab frame using Eq's 29,19, and 15 for each valid omega pair
if isempty(omega) == 1
    return
else
    for ii = 1:2
        R_S2L = [cos(omega(ii))  0  sin(omega(ii));
            sin(chi)*sin(omega(ii))  cos(chi)  -sin(chi)*cos(omega(ii));
            -cos(chi)*sin(omega(ii))  sin(chi)  cos(chi)*cos(omega(ii))];
        
        G_L = R_S2L * Vinv_S * R_C2S * B * Gstar;  Ghat_L = unit(G_L);
        % dhat is the diffracted beam direction
        dhat_L = (eye(3) - 2 * dyad(Ghat_L, Ghat_L)) * bhat_L;
        
        
        %% Intensity
        Int = 100;
        
        
        %% Calculate [P4]_D from Eq. 28
        zeta0(ii,:) = [-(dhat_L(1)/dhat_L(3))*DISTANCE  -(dhat_L(2)/dhat_L(3))*DISTANCE]; % Pagan/Miller Eq. 10
        
        p = R_S2L * tVec_C';
        px = p(1);  py = p(2);  pz = p(3);
        zeta_temp = zeta0(ii,:) + [px-(dhat_L(1)/dhat_L(3))*pz   py-(dhat_L(2)/dhat_L(3))*pz];
        
        zeta = vertcat(zeta, zeta_temp);
    end

    
    Lights = [vertcat(omega(1,1), omega(2,1))  ttheta*ones(2,1)  zeta  Int*ones(2,1)];
end