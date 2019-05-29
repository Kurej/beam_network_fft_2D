function [ kx, ky, x, y ] = calcGrids( nbPoints1D, taille1D )
%   calcGrids.m : Defining grids for calculation
%
%   Input parameters :
%       * nbPoints1D : number of points in one dimension (output grid will 
%           have this dimension squared). A power of two fastens the fft
%           calculation.
%       * size1D : Spatial dimension of the grid [m]. Should be at least 5
%           to 10 times larger than the pupil to study to obtain a good far
%           field resolution.
%
%   Output parameters :
%       * x : spatial grid after meshgrid (1) [m]
%       * y : spatial grid after meshgrid (2) [m]
%       * kx : spatial frequency grid after meshgrid (1) [rad/m]
%       * ky : spatial frequency grid after meshgrid (2) [rad/m]

    %%% Spatial grid
    dX = taille1D/(nbPoints1D-1) ; % Sample size [m]
    limX = (nbPoints1D/2)*dX ; % Range limits [m]
    X = -limX:dX:limX-dX ; % Spatial vector [m]
    [x,y] = meshgrid(X) ; % Spatial grid [m]
    
    %%% Spatial frequency grid
    dNx = 1/taille1D ; % Fourier sample size [1/m]
    limNx = (nbPoints1D/2)*dNx ; % Fourier range limits [1/m]
    NX = (-limNx:dNx:limNx-dNx) ; % Spatial frequencies vector [1/m]
    KX = 2*pi*NX ; % Spatial frequencies vector [rad/m]
    [kx,ky] = meshgrid(KX) ; % Spatial frequencies grid [rad/m]
end

