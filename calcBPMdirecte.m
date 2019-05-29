function [ champPropagDirect ] = calcBPMdirecte( champInit, deltaZVoulu, CP, grid )
% calcBPMdirecte.m : calcul de diffraction via Beam Propagation Method (une seule itération)
%
%   Paramètres d'entrée :
%       * champInit : champ initial à diffracter
%       * deltaZVoulu : plans de calcul choisis [m]
%       * CP : champ proche et paramètres associés
%       * grid.kx : grille spatiale conjuguée après meshgrid (1) [rad/m]
%       * grid.ky : grille spatiale conjuguée après meshgrid (2) [rad/m]
%       * grid.x : grille spatiale après meshgrid (1) [m]
%       * grid.y : grille spatiale après meshgrid (2) [m]
%
%   Paramètres de sortie :
%       * champPropagDirect : champ calculé sur les plans choisis
    
    champPropagDirect = zeros(length(grid.x),length(grid.x),length(deltaZVoulu)) ;
    for i=1:length(deltaZVoulu)
        champPropagDirect(:,:,i) = fftshift(fft2(fftshift(champInit))) ;
        champPropagDirect(:,:,i) = champPropagDirect(:,:,i).*exp(1i*deltaZVoulu(i)*sqrt(4*pi^2/CP.lambda^2-grid.kx.^2-grid.ky.^2)) ;
        champPropagDirect(:,:,i) = ifftshift(ifft2(ifftshift(champPropagDirect(:,:,i)))) ;
    end
end

