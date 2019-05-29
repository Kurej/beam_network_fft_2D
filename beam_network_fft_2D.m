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

CP.NB = 19 ; % Nombre total de faisceaux
CP.lambda = 1064e-9 ; % Longueur d'onde du rayonnement [m]
CP.pitch = 1.6e-3 ; % Pas de la maille carrée [m] (sortie réducteur 1.6 mm)
CP.w0 = 0.6e-3 ; % Rayon à 1/e² en intensité d'un faisceau gaussien en champ proche [m] (sortie réducteur 0.6 mm) (225µm)
CP.DLens = 1.5e-3 ; % Diamètre d'une lentille de collimation du champ proche [m] (sortie réducteur 1.5 mm) (1.125 mm/2)
CP.maille = 'hexagonale' ; % Arrangement de la maille du champ proche : 'carree', 'hexagonale', ou 'lineaire'
CP.stdPower = 0/100 ;
CP.power = abs(ones(CP.NB,1).*(1 + CP.stdPower*randn(CP.NB,1))) ;

CP.pistonRange = pi ;
CP.tipRange = 1e-3 ;
CP.tiltRange = 1e-3 ;

CP.piston = CP.pistonRange*(-1 + 2*rand(CP.NB,1)) ;
CP.tip = CP.tipRange*(-1 + 2*rand(CP.NB,1)) ;
CP.tilt = CP.tiltRange*(-1 + 2*rand(CP.NB,1)) ;
CP.defocus = zeros(CP.NB,1) ;

switch CP.maille
    case 'carree'
        if mod(sqrt(CP.NB),1)~=0
            error('Veuillez entrer un nombre de faisceaux valable (4, 9, 16, 25, 36, 49, 64, 81, 100, ...) !')
        else
            [ CP.faisc, CP.posFaisc ] = champProcheMailleCarree( CP, grd ) ;
        end
    case 'hexagonale'
            [ CP.faisc, CP.posFaisc ] = champProcheMailleHexa( CP, grd ) ;
    case 'lineaire'
        [ CP.faisc, CP.posFaisc ] = champProcheMailleLineaire( CP, grd ) ;
    otherwise
        error('Type de maille inexistant !')
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
        axis square,colorbar,shading flat,colormap viridis,axis((sqrt(CP.NB)+1)*CP.pitch*[-1 1 -1 1])
        xlabel('x [m]'),ylabel('y [m]'),title('Intensité champ proche'),caxis([0 1])

    subplot(2,2,2),cla
        imagesc(grd.x(1,:),grd.y(:,1),angle(champ))
        axis square,colorbar,shading flat,colormap viridis,axis((sqrt(CP.NB)+1)*CP.pitch*[-1 1 -1 1])
        xlabel('x [m]'),ylabel('y [m]'),title('Phase champ proche'),caxis(pi*[-1 1])
        drawnow

        
%% Diffraction de Frauhofer
TFchamp = fftshift(fft2(fftshift(champ))) ;

figure(1)
    subplot(2,2,3),cla
    imagesc(grd.kx(1,:),grd.ky(:,1),abs(TFchamp).^2)
    axis square,colorbar,shading flat,axis((2*pi)/CP.w0*[-1 1 -1 1])
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


