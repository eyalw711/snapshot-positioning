function trops = tropospheric_model(els)
% returns the tropospheric delay model in meters
trops = arrayfun(@(x) tropo(sin(x*pi/180),0.0,1013.0,293.0,50.0, 0.0,0.0,0.0), els);
end
