function P=PIn(obs,t)
    Source=[2,0];
    strength=16.66;c=340;K=1;theta=1.3734;
    rObs=norm(obs-Source);
    P=(strength/rObs)*cos(c*t-K*rObs+theta);
end

