function P = PSca(obs,k)
    P=0;c=340;K=1;
    tStart=0;tDelta=0.04;
    eTemp=evalin('base','eTemp');
    xEqs=evalin('base','xEqs');yEqs=evalin('base','yEqs');
    tNow=(k-1)*tDelta+tStart;
    for i=1:4
        r=norm([xEqs(i),yEqs(i)]-obs);
        P=P+eTemp(i,k)*sin(c*tNow-K*r)/(4*pi*r); 
    end
end

