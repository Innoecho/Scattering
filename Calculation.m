clear;clc;
%% ��Դ����
%��Դ�������ܶȣ����٣�������������Դ�뾶��Դ���ٶȷ�ֵ��Դǿ��A
rho=1.225;c=340;lambda=6.28;K=1;r0=0.2;amp=1;strength=rho*c*K*r0*r0*amp;
theta=atan(1/(K*r0));omega=K*c;
%��Դλ��
xSource=2;ySource=0;
%% ɢ������
nPoint=360;                               %ɢ����һȦ�ܵ���
nSca=6;                                 %ɢ������
angle=zeros(nPoint,1);xSca=zeros(nPoint,1);ySca=zeros(nPoint,1);  %ɢ��㷽λ�Ǽ����꣨��ɢ����Ŵ��棩
r=zeros(nPoint,1);                %ɢ��㵽Դ����루��ɢ����Ŵ��棩
tStart=0;tStep=40;tDelta=0.04;tNow=0;    %��ʼʱ��,ʱ���ܲ�����ʱ����,��ǰʱ��
pIn=zeros(nPoint,tStep);            %ɢ����������ѹ����ɢ���ı�ź͵�ǰʱ�䲽������
xVect=zeros(nSca,1);yVect=zeros(nSca,1);
pGrad=zeros(nSca,tStep);  %ɢ���������ѹ�ݶ��뷨ʸ��ˡ���ɢ����ź�ʱ�䲽������
%% ��ЧԴ�����
nEqs=4;                          %��ЧԴ����� ӦС��ɢ������
kZoom=0.2;                             %��ЧԴ������ϵ��
xEqs=zeros(nEqs,1);yEqs=zeros(nEqs,1);   %��ЧԴ������
eTemp=zeros(nEqs,tStep);          %SVDÿ���������e��i��������eTemp(i,k),���ഢ�浱ǰʱ�䲽��
%% �۲�����
xObs=-1.2;yObs=0;rObs=0;pObs=0;        %�۲�����꣬��Դ����룬������ѹ
eEO=zeros(nEqs,1);            %�۲�㵽��ЧԴ����룬�õ�ЧԴ��Ŵ���
pScaObs=0;                  %�۲��ɢ����ѹ����ЧԴ�����ĺͣ����Լ��м�������ÿ��Դ��ɢ����ѹ
pSca=zeros(nEqs,1);
%% ɢ������ЧԴ������
rEqsSca=zeros(nEqs,nSca);           %�õ�ЧԴ�����ɢ����Ŵ��棬�ǳ��ȣ�����
rX=zeros(nEqs,nSca);rY=zeros(nEqs,nSca);rN=zeros(nEqs,nSca);   %��ЧԴ��ָ��ɢ���ĵ�λʸ�����õ�ЧԴ��ź�ɢ����Ŵ���
tTrans=zeros(nEqs,nSca);tau=0;      %��ЧԴ���ɢ���һһ��Ӧ�Ĵ���ʱ�䣬��ЧԴ����ӳ�ʱ�䣨����ʱ�̣�
%% ������
coef1=zeros(nEqs,nSca);coef2=zeros(nEqs,nSca);
A=zeros(nSca,nEqs);B=zeros(nSca,1);   %A�����ɵ�ЧԴ��ź�ɢ����Ŵ��档B������ɢ����Ŵ���
e=zeros(nEqs,1);            %��ЧԴ��Դǿ��e���õ�ЧԴ���Ŵ���
%% 360�����������ѹ
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
%% �۲���������ѹ
obsInPressure=zeros(tStep,2);
for i=1:tStep
    tNow=(i-1)*tDelta+tStart;
    rObs=sqrt((xObs-xSource)^2+(yObs-ySource)^2);
    pObs=(strength/rObs)*cos(omega*tNow-K*rObs+theta);
    obsInPressure(i,:)=[tNow,pObs];
end
%% 6��ɢ���ķ�ʸ(ɢ����Ӧ���Ϊ30,90,150,210,270,330)
for i=1:nSca
    xVect(i)=-(ySca(60*i-30+1)-ySca(60*i-30-1));
    yVect(i)=xSca(60*i-30+1)-xSca(60*i-30-1);
    rTemp=sqrt(xVect(i)^2+yVect(i)^2);
    xVect(i)=-xVect(i)/rTemp;
    yVect(i)=-yVect(i)/rTemp;
end
%% 6��ɢ�����ݶ�,�ݶ��뷨ʸ���pGrad
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
%% ���㷽�����ϵ��
for i=1:nEqs
    xEqs(i)=kZoom*cos(((i-1)*360./nEqs)*pi/180);
    yEqs(i)=kZoom*sin(((i-1)*360./nEqs)*pi/180);
    for j=1:nSca
        iSca=(j-1)*60+30;             %�ҵ�6��ɢ���ı��   ɢ����������һ��360��
        rEqsSca(i,j)=sqrt((xEqs(i)-xSca(iSca))^2+(yEqs(i)-ySca(iSca))^2);
        tTrans(i,j)=rEqsSca(i,j)/c;
        rX(i,j)=(xSca(iSca)-xEqs(i))/rEqsSca(i,j);
        rY(i,j)=(ySca(iSca)-yEqs(i))/rEqsSca(i,j);     
        rN(i,j)=rX(i,j)*xVect(j)+rY(i,j)*yVect(j);   %ɢ��㷨ʸֻ����6��
        coef1(i,j)=rN(i,j)/(c*rEqsSca(i,j));
        coef2(i,j)=rN(i,j)/(rEqsSca(i,j)^2);
    end
end
%% �۲�㵽��ЧԴ�����
for i=1:nEqs
    eEO(i)=sqrt((xEqs(i)-xObs)^2+(yEqs(i)-yObs)^2);
end
%% �����A
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
%% �����B������Դǿ��
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
    % ����SVD�ӹ�����ȡe[nEqs]
    e=A\B;
    % Ϊ������ɢ����ѹ������ã���ʱ�䲽�������ȥ
    eTemp(1,k)=e(1);           
    eTemp(2,k)=e(2);
    eTemp(3,k)=e(3);
    eTemp(4,k)=e(4);
end
%% ������õ�ЧԴǿ�� �õ�����Դ��ʽ��ɢ����ѹ
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
%% ��ͼ
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
figure(1);set(gcf,'name','������ѹ');
contourf(x,y,Pin);
axis equal
axis([-5,5,-5,5]);
figure(2);set(gcf,'name','ɢ����ѹ');
contourf(x,y,Psca);
axis equal
axis([-5,5,-5,5]);
figure(3);set(gcf,'name','����ѹ');
contourf(x,y,Psca+Pin);
axis equal
axis([-5,5,-5,5]);