clear;clc;
%% 点源参数
%点源参数：密度，声速，波长，波数，源半径，源振动速度幅值，源强度A
rho=1.225;c=340;lambda=6.28;K=1;r0=0.2;amp=1;strength=rho*c*K*r0*r0*amp;
theta=atan(1/(K*r0));omega=K*c;
%点源位置
xSource=2;ySource=0;
%% 散射点参数
nPoint=360;                               %散射体一圈总点数
nSca=6;                                 %散射点个数
angle=zeros(nPoint,1);xSca=zeros(nPoint,1);ySca=zeros(nPoint,1);  %散射点方位角及坐标（以散射点编号储存）
r=zeros(nPoint,1);                %散射点到源点距离（以散射点编号储存）
tStart=0;tStep=40;tDelta=0.04;tNow=0;    %开始时间,时间总步数，时间间隔,当前时间
pIn=zeros(nPoint,tStep);            %散射点的入射声压，以散射点的编号和当前时间步数储存
xVect=zeros(nSca,1);yVect=zeros(nSca,1);
pGrad=zeros(nSca,tStep);  %散射点入射声压梯度与法矢点乘。用散射点编号和时间步数储存
%% 等效源面参数
nEqs=4;                          %等效源点个数 应小于散射点个数
kZoom=0.2;                             %等效源面缩放系数
xEqs=zeros(nEqs,1);yEqs=zeros(nEqs,1);   %等效源点坐标
eTemp=zeros(nEqs,tStep);          %SVD每步算出来的e（i）储存在eTemp(i,k),即多储存当前时间步长
%% 观察点参数
xObs=-1.2;yObs=0;rObs=0;pObs=0;        %观察点坐标，距源点距离，入射声压
eEO=zeros(nEqs,1);            %观察点到等效源点距离，用等效源编号储存
pScaObs=0;                  %观察点散射声压（等效源产生的和），以及中间计算变量每个源的散射声压
pSca=zeros(nEqs,1);
%% 散射点与等效源点间参数
rEqsSca=zeros(nEqs,nSca);           %用等效源编号与散射点编号储存，是长度，标量
rX=zeros(nEqs,nSca);rY=zeros(nEqs,nSca);rN=zeros(nEqs,nSca);   %等效源点指向散射点的单位矢量，用等效源编号和散射点编号储存
tTrans=zeros(nEqs,nSca);tau=0;      %等效源点和散射点一一对应的传播时间，等效源点的延迟时间（发生时刻）
%% 方程项
coef1=zeros(nEqs,nSca);coef2=zeros(nEqs,nSca);
A=zeros(nSca,nEqs);B=zeros(nSca,1);   %A矩阵，由等效源编号和散射点编号储存。B矩阵由散射点编号储存
e=zeros(nEqs,1);            %等效源的源强度e，用等效源点编号储存
%% 360个点的入射声压
for i=1:nPoint
    angle(i)=nPoint*i/360;
    xSca(i)=cos(angle(i)*pi/180);
    ySca(i)=sin(angle(i)*pi/180);
    r(i)=sqrt((xSca(i)-xSource)^2+(ySca(i)-ySource)^2);
    for j=1:tStep
        tNow=(j-1)*tDelta+tStart;
        pIn(i,j)=(strength/r(i))*cos(omega*tNow-K*r(i)+theta);
    end
end
%% 观察点的入射声压
obsInPressure=zeros(tStep,2);
for i=1:tStep
    tNow=(i-1)*tDelta+tStart;
    rObs=sqrt((xObs-xSource)^2+(yObs-ySource)^2);
    pObs=(strength/rObs)*cos(omega*tNow-K*rObs+theta);
    obsInPressure(i,:)=[tNow,pObs];
end
%% 6个散射点的法矢(散射点对应编号为30,90,150,210,270,330)
for i=1:nSca
    xVect(i)=-(ySca(60*i-30+1)-ySca(60*i-30-1));
    yVect(i)=xSca(60*i-30+1)-xSca(60*i-30-1);
    rTemp=sqrt(xVect(i)^2+yVect(i)^2);
    xVect(i)=-xVect(i)/rTemp;
    yVect(i)=-yVect(i)/rTemp;
end
%% 6个散射点的梯度,梯度与法矢点乘pGrad
for i=1:nSca
    for j=1:tStep
        index=(i-1)*60+30;delta=0.01;
        r1=sqrt((xSca(index)-xSource)^2+(ySca(index)+delta-ySource)^2);
        r2=sqrt((xSca(index)-xSource)^2+(ySca(index)-delta-ySource)^2);
        r3=sqrt((xSca(index)-delta-xSource)^2+(ySca(index)-ySource)^2);
        r4=sqrt((xSca(index)+delta-xSource)^2+(ySca(index)-ySource)^2);
        
        tNow=(j-1)*tDelta+tStart;
        p1=(strength/r1)*cos(omega*tNow-K*r1+theta);
        p2=(strength/r2)*cos(omega*tNow-K*r2+theta);
        p3=(strength/r3)*cos(omega*tNow-K*r3+theta);
        p4=(strength/r4)*cos(omega*tNow-K*r4+theta);

        pGradX=(p4-p3)/(2*delta);
        pGradY=(p1-p2)/(2*delta);
        pGrad(i,j)=pGradX*xVect(i)+pGradY*yVect(i);
    end
end
%% 计算方程相关系数
for i=1:nEqs
    xEqs(i)=kZoom*cos(((i-1)*360./nEqs)*pi/180);
    yEqs(i)=kZoom*sin(((i-1)*360./nEqs)*pi/180);
    for j=1:nSca
        iSca=(j-1)*60+30;             %找到6个散射点的编号   散射点坐标存了一周360个
        rEqsSca(i,j)=sqrt((xEqs(i)-xSca(iSca))^2+(yEqs(i)-ySca(iSca))^2);
        tTrans(i,j)=rEqsSca(i,j)/c;
        rX(i,j)=(xSca(iSca)-xEqs(i))/rEqsSca(i,j);
        rY(i,j)=(ySca(iSca)-yEqs(i))/rEqsSca(i,j);     
        rN(i,j)=rX(i,j)*xVect(j)+rY(i,j)*yVect(j);   %散射点法矢只存了6个
        coef1(i,j)=rN(i,j)/(c*rEqsSca(i,j));
        coef2(i,j)=rN(i,j)/(rEqsSca(i,j)^2);
    end
end
%% 观察点到等效源点距离
for i=1:nEqs
    eEO(i)=sqrt((xEqs(i)-xObs)^2+(yEqs(i)-yObs)^2);
end
%% 求矩阵A
for j=1:nSca
    for i=1:nEqs
        if tDelta-tTrans(i,j)>0
            A(j,i)=coef1(i,j)/tDelta+coef2(i,j)*(tDelta-tTrans(i,j))/tDelta;
        else
            A(j,i)=0;
        end
        A(j,i)=A(j,i)/(4*pi);
    end
end
%% 求矩阵B并计算源强度
for k=2:tStep
    tNow=(k-1)*tDelta+tStart;
    for j=1:nSca
		B(j)=0;
        for i=1:nEqs
            for l=1:k-1
                B(j)=B(j)+(coef1(i,j)*dphi(l,tNow-tTrans(i,j))+coef2(i,j)*phi(l,tNow-tTrans(i,j)))*eTemp(i,l);
            end
        end
        B(j)=pGrad(j,k)-B(j)/(4*pi);
    end
    % 调用SVD子过程求取e[nEqs]
    e=A\B;
    % 为后面求散射声压方便调用，把时间步长储存进去
    eTemp(1,k)=e(1);           
    eTemp(2,k)=e(2);
    eTemp(3,k)=e(3);
    eTemp(4,k)=e(4);
end
%% 根据求得等效源强度 用单极子源公式求散射声压
obsScaPressure=zeros(tStep,2);
for k=2:tStep
    pScaObs=0;
    for i=1:nEqs 
        tNow=(k-1)*tDelta+tStart;
        pSca(i)=eTemp(i,k)*sin(omega*tNow-K*eEO(i))/(4*pi*(eEO(i))); 
        pScaObs=pScaObs+pSca(i);      
    end
    obsScaPressure(k,:)=[tNow,pScaObs];
end
%% 云图
x=linspace(-5,5,100);
y=linspace(-5,5,100);
Pin=zeros(100,100);
Psca=zeros(100,100);
for i=1:100
    for j=1:100
        Pin(i,j)=PIn([x(i),y(j)],1.6);
        Psca(i,j)=PSca([x(i),y(j)],40);
    end
end
figure(1);set(gcf,'name','入射声压');
contourf(x,y,Pin);
axis equal
axis([-5,5,-5,5]);
figure(2);set(gcf,'name','散射声压');
contourf(x,y,Psca);
axis equal
axis([-5,5,-5,5]);
figure(3);set(gcf,'name','总声压');
contourf(x,y,Psca+Pin);
axis equal
axis([-5,5,-5,5]);