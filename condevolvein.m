function newConductance = condevolvein(previousConductance, vin,vt1,vt2,levels,gmax,gmin,mult)  
    dg = (gmax-gmin)/(levels-1);
    if vin > vt1
         newConductance = previousConductance + mult*dg;
    elseif abs(vin) > vt2
        newConductance = previousConductance - mult*dg;
    else
        newConductance = previousConductance;
    end
    if newConductance > gmax
        newConductance = previousConductance;
    elseif  newConductance < gmin
         newConductance = previousConductance;
    end
end

