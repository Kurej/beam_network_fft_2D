%{
    Fourier or Fresnel transform of a laser beam array with square or
    hexagonal lattice geometries.

    Low order aberrations such as piston, tip, tilt and defocus can be
    individually addressed to each beam.

    Assumptions :
        * Individual beams are coherent (monochromatic spectrum)
        * They share the same polarization state

    Sections of the program:
        * Defining grids in near field and far field
        * Creating near field
        * Calculating far field
        * Calculating intermediate field
%}
clear all ;


%% Defining grids in near field and far field   
grd.nbPoints = 1024*4 ;
grd.step = 40e-6 ;
grd.gridSize = grd.step*grd.nbPoints ;
[ grd.kx, grd.ky, grd.x, grd.y ] = calcGrids( grd.nbPoints, grd.gridSize ) ;


%% Définition du champ proche
if exist('CP','var'), clear CP, end
if exist('champ','var'), clear champ, end

CP.NB = 19 ; % Total number of beams
CP.lambda = 1064e-9 ; % Wavelength [m]
CP.pitch = 1.6e-3 ; % Network pitch [m]
CP.w0 = 0.6e-3 ; % 1/e² intensity beam radius [m]
CP.DLens = 1.5e-3 ; % Collimation lens diameter [m]
CP.lattice = 'hexagonal' ; % Network lattice : 'linear', 'square', or 'hexagonal'
CP.stdPower = 0/100 ; % Standard deviation of beams power
CP.power = abs(ones(CP.NB,1).*(1 + CP.stdPower*randn(CP.NB,1))) ; % Beams power

CP.pistonRange = 0 ; % Piston range for uniform law [rad]
CP.tipRange = 1e-3 ; % Tip range for uniform law [rad]
CP.tiltRange = 1e-3 ; % Tilt range for uniform law [rad]

rng('shuffle') ; % Randomizes the rng seed
CP.piston = CP.pistonRange*(-1 + 2*rand(CP.NB,1)) ; % Beams pistons
CP.tip = CP.tipRange*(-1 + 2*rand(CP.NB,1)) ; % Beams tips
CP.tilt = CP.tiltRange*(-1 + 2*rand(CP.NB,1)) ; % Beams tilts
CP.defocus = zeros(CP.NB,1) ; % Beams defocus

switch CP.lattice
    case 'square'
        [ CP.faisc, CP.posFaisc ] = champProcheMailleCarree( CP, grd ) ;
    case 'hexagonal'
        [ CP.faisc, CP.posFaisc ] = champProcheMailleHexa( CP, grd ) ;
    case 'linear'
        [ CP.faisc, CP.posFaisc ] = champProcheMailleLineaire( CP, grd ) ;
    otherwise
        error('Unknown lattice type !')
end

champ = zeros(grd.nbPoints) ;
for i=1:CP.NB
    champ(CP.faisc.ind{i}) = champ(CP.faisc.ind{i}) + CP.power(i)*CP.faisc.val{i} ;
    champ(CP.faisc.ind{i}) = champ(CP.faisc.ind{i}).*exp(1i*CP.piston(i)) ...
                             .*exp(1i*(CP.tip(i)/CP.lambda*(grd.y(CP.faisc.ind{i})-CP.posFaisc.y(i)))) ...
                             .*exp(1i*(CP.tilt(i)/CP.lambda*(grd.x(CP.faisc.ind{i})-CP.posFaisc.x(i)))) ... 
                             ;
    
end

figure(1)
    subplot(2,2,1),cla
        imagesc(grd.x(1,:),grd.y(:,1),abs(champ).^2)
        axis square,colorbar,shading flat,axis((sqrt(CP.NB)+1)*CP.pitch*[-1 1 -1 1])
        xlabel('x [m]'),ylabel('y [m]'),title('Intensité champ proche'),caxis([0 1])

    subplot(2,2,2),cla
        imagesc(grd.x(1,:),grd.y(:,1),angle(champ))
        axis square,colorbar,shading flat,axis((sqrt(CP.NB)+1)*CP.pitch*[-1 1 -1 1])
        xlabel('x [m]'),ylabel('y [m]'),title('Phase champ proche'),caxis(pi*[-1 1])
        drawnow

        
%% Fraunhofer diffraction
TFchamp = fftshift(fft2(fftshift(champ))) ;

figure(1)
    subplot(2,2,3),cla
    imagesc(grd.kx(1,:),grd.ky(:,1),abs(TFchamp).^2)
    axis square,colorbar,axis((2*pi)/CP.w0*[-1 1 -1 1])
    xlabel('kx [rad/m]'),ylabel('ky [rad/m]'),title(['Intensité CL'])
    drawnow
    
        

%% Diffraction de Fresnel du champ diffracté aux distances voulues
% dz = 2 ;
% champZ = calcBPMdirecte( champ, dz, CP, grd ) ;
% 
% figure(1)
%     subplot(2,2,4),cla
%     imagesc(grd.x(1,:),grd.y(:,1),abs(champZ).^2)
%     axis square,colorbar,shading flat,axis((sqrt(CP.NB)+1)*CP.pitch*[-1 1 -1 1])
%     xlabel('x [m]'),ylabel('y [m]'),title(['Intensité CD z=' num2str(dz(1)) 'm'])
%     drawnow


