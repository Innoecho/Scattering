function res=dphi(l,t)
    tStart=0;tDelta=0.04;
    taul=tStart+tDelta*(l-1);
    tau=tStart+tDelta*l;
    taur=tStart+tDelta*(l+1);
    if (taul<=t) && (t<=tau)
        res=1/tDelta;
    elseif (tau<=t) && (t<=taur)
        res=-1/tDelta;
    else
        res=0;
    end
end
    
    