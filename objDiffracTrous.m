function [ transmittance ] = objDiffracTrous( objDiffrac, CP, grid )
%   objDiffracTrousCirc.m : D�finition d'un objet diffractant en champ
%   proche
%
%   Param�tres d'entr�e :
%       * objDiffrac.diamTrous : diam�tre des trous diffractants de l'objet diffractant [m]
%       * CP.NB : nombre total de faisceaux en champ proche
%       * CP.lambda : longueur d'onde du rayonnement [m]
%       * CP.pitch : pas du r�seau [m]
%       * CP.w0 : rayon � 1/e� en intensit� d'un faisceau gaussien en champ proche [m]
%       * CP.DLens : diam�tre d'une lentille de collimation du champ proche [m]
%       * CP.maille : type d'arrangement de la maille du r�seau
%       * grid : grilles spatiales apr�s meshgrid [m, rad/m]
%
%   Param�tres de sortie :
%       * transmittance : transmittance complexe de l'objet diffractant

    %%% D�finition de l'objet diffractant
    transmittance = zeros(length(grid.x)) ; % Initialisation � un cache opaque
    transPhasePentes = zeros(length(grid.x)) ; % Initialisation � des pentes de phase plates
    courbPhase = zeros(length(grid.x)) ; % Initialisation � une courbure nulle
    ind = cell(CP.NB,1) ;
    rng('shuffle') ;
    var.x = (-1 + 2*rand(CP.NB,1))*objDiffrac.varPosTrous ;
    var.y = (-1 + 2*rand(CP.NB,1))*objDiffrac.varPosTrous ;
    var.diam = (-1 + 2*rand(CP.NB,1))*objDiffrac.varDiamTrous ;
 
    for it=1:CP.NB
        ind{it} = find( (grid.x-CP.posFaisc.x(it)-objDiffrac.transPosTrous.x-var.x(it)).^2+(grid.y-CP.posFaisc.y(it)-objDiffrac.transPosTrous.y-var.y(it)).^2 <= ((objDiffrac.diamTrous+var.diam(it))/2).^2 ) ;
        transmittance(ind{it}) = 1 ;
        if objDiffrac.penteCoeff
            transPhasePentes(ind{it}) = objDiffrac.penteCoeff ...
                                        * ( ...
                                                (CP.posFaisc.x(it)+objDiffrac.transPosTrous.x+var.x(it))*grid.x(ind{it}) ...
                                                + (CP.posFaisc.y(it)+objDiffrac.transPosTrous.y+var.y(it))*grid.y(ind{it}) ...
                                                - (CP.posFaisc.x(it)+objDiffrac.transPosTrous.x+var.x(it))^2 ...
                                                - (CP.posFaisc.y(it)+objDiffrac.transPosTrous.y+var.y(it))^2 ...
                                            ) ;
        end
        if objDiffrac.courbCoeff
            courbPhase(ind{it}) = ( (grid.x(ind{it})-CP.posFaisc.x(it)).^2/objDiffrac.courbCoeff^2 ...
                                  + (grid.y(ind{it})-CP.posFaisc.y(it)).^2/objDiffrac.courbCoeff^2 ...
                                        ) ;
        end
        transmittance(ind{it}) = transmittance(ind{it}) ...
                                    .*exp(1i*(-1+2*rand)*abs(objDiffrac.bornePhaseTrous)) ...
                                    .*exp(1i*transPhasePentes(ind{it})) ...
                                    .*exp(1i*courbPhase(ind{it})) ;
    end
end

