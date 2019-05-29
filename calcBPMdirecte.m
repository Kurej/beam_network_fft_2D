function [ champPropagDirect ] = calcBPMdirecte( champInit, deltaZVoulu, CP, grid )
% calcBPMdirecte.m : calcul de diffraction via Beam Propagation Method (une seule it�ration)
%
%   Param�tres d'entr�e :
%       * champInit : champ initial � diffracter
%       * deltaZVoulu : plans de calcul choisis [m]
%       * CP : champ proche et param�tres associ�s
%       * grid.kx : grille spatiale conjugu�e apr�s meshgrid (1) [rad/m]
%       * grid.ky : grille spatiale conjugu�e apr�s meshgrid (2) [rad/m]
%       * grid.x : grille spatiale apr�s meshgrid (1) [m]
%       * grid.y : grille spatiale apr�s meshgrid (2) [m]
%
%   Param�tres de sortie :
%       * champPropagDirect : champ calcul� sur les plans choisis
    
    champPropagDirect = zeros(length(grid.x),length(grid.x),length(deltaZVoulu)) ;
    for i=1:length(deltaZVoulu)
        champPropagDirect(:,:,i) = fftshift(fft2(fftshift(champInit))) ;
        champPropagDirect(:,:,i) = champPropagDirect(:,:,i).*exp(1i*deltaZVoulu(i)*sqrt(4*pi^2/CP.lambda^2-grid.kx.^2-grid.ky.^2)) ;
        champPropagDirect(:,:,i) = ifftshift(ifft2(ifftshift(champPropagDirect(:,:,i)))) ;
    end
end

