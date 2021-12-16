function y = skiba_production_function(G, param)

    yH = param.AH * max(G.k - param.kappa,0).^param.alpha; 
    yL = param.AL * G.k.^param.alpha;
    y  = max(yH, yL);

end