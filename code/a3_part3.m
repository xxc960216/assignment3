% ELEC4700 - Assignment 3 part 3
% Xiaochen Xin 100989338

clearvars
clearvars -GLOBAL

global C
global X Y

    C.q_0 = 1.60217653e-19;             % electron charge
    C.hb = 1.054571596e-34;             % Dirac constant
    C.h = C.hb * 2 * pi;                % Planck constant
    C.m_0 = 9.10938215e-31;             % electron mass
    C.kb = 1.3806504e-23;               % Boltzmann constant
    C.eps_0 = 8.854187817e-12;          % vacuum permittivity
    C.mu_0 = 1.2566370614e-6;           % vacuum permeability
    C.c = 299792458;                    % speed of light
    C.g = 9.80665;                      % metres (32.1740 ft) per s²
    
mn=0.26*C.m_0; %electron mass
Temp = 300; %Given in kelvin
MTBC = 0.2e-12;
Vleft = 0.1;%voltage of left side
electronConc = 10e15;
s1 = 1;%for resistances
s2 = 0.01;

%thermal velocity
Vth = sqrt(2*C.kb*Temp/mn);

l=200*10^-9;
w=100*10^-9;
area =l*w;

size=1000;
numOfAtom=10;

X= rand(2,size);
Y= rand(2,size);

xr(1,:)= X(1,:)*l;
yr(1,:)= Y(1,:)*w;

box1left = xr>0.8e-7;
box1right = xr<1.2e-7;
box1 = box1left & box1right;
box2bottom = yr<0.4e-7;
box1bottom = box2bottom & box1;
box2top = yr>0.6e-7;
boxtop = box2top & box1;

checkboxes = boxtop | box1bottom;
while(sum(checkboxes)>0)
    
    xr(checkboxes) = rand*l;
    yr(checkboxes) = rand*w;
    
    box1left = xr>0.8e-7;
    box1right = xr<1.2e-7;
    box1 = box1left & box1right;
    box2bottom = yr<0.4e-7;
    box1bottom = box2bottom & box1;
    box2top = yr>0.6e-7;
    boxtop = box2top & box1;

    checkboxes = boxtop | box1bottom;
end
colour = rand(1,numOfAtom);

spacStep = 0.01*w;
dt = spacStep/Vth;
steps = 1000;

Vthn = Vth/sqrt(2);
xv = Vthn*randn(1,size);
yv = Vthn*randn(1,size);

xv(1,:) = xv(1,:)*dt;
yv(1,:) = yv(1,:)*dt;

Pscat=1-exp(-(dt/MTBC));

MFPcount = zeros(1,size);

Efield = Vleft/l;
force = Efield*C.q_0;
acceleration = force/mn;
accelVelocity = acceleration*(dt^2);

xbox = [0.8e-7 1.2e-7];
yboxbottom = [0 0.4e-7];
yboxtop = [0.6e-7 w];

squares = 100;
xResolution = l/squares;
yResolution = w/squares;

[MX,MY] = meshgrid(0:xResolution:l,0:yResolution:w);

xBoxlogic = MX>=xbox(1) & MX<=xbox(2);
yBoxlogic = MY>=yboxtop(1) | MY<=yboxbottom(2);
Boxlogic = xBoxlogic & yBoxlogic;

Smap = zeros(squares+1);
Smap(Boxlogic) = s2;
Smap(~Boxlogic) = s1;

voltage = 0.8;
G = sparse(squares+1);
B = zeros(squares+1,1);
for i =1:1:squares
    for j =1:1:squares
        n = j+(i-1)*squares;
        nxm = j+(i-2)*squares;
        nxp = j+i*squares;
        nyp = j+1+ (i-1)*squares;
        nym = j-1+ (i-1)*squares;
        
        if(i==1)
            G(n,:) = 0;
            G(n,n) = Smap(i,j);
            B(n) = voltage;
        elseif(i==squares)
            G(n,:) = 0;
            G(n,n) = Smap(i,j);
            B(n) = 0;
        elseif(j==1)
            G(n,:) = 0;
            G(n,nxm) = (Smap(i-1,j)+Smap(i,j))/2;
            G(n,nxp) = (Smap(i+1,j)+Smap(i,j))/2;
            G(n,nyp) = (Smap(i,j+1)+Smap(i,j))/2;
            G(n,n) = -(G(n,nxm)+G(n,nxp)+G(n,nyp));
        elseif(j==squares)
            G(n,:) = 0;
            G(n,nxm) = (Smap(i-1,j)+Smap(i,j))/2;
            G(n,nxp) = (Smap(i+1,j)+Smap(i,j))/2;
            G(n,nym) = (Smap(i,j-1)+Smap(i,j))/2;
            G(n,n) = -(G(n,nxm)+G(n,nxp)+G(n,nym));
        else          
            G(n,:) = 0;
            G(n,nxm) = (Smap(i-1,j)+Smap(i,j))/2;
            G(n,nxp) = (Smap(i+1,j)+Smap(i,j))/2;
            G(n,nyp) = (Smap(i,j+1)+Smap(i,j))/2;
            G(n,nym) = (Smap(i,j-1)+Smap(i,j))/2;
            G(n,n) = -(G(n,nxm)+G(n,nxp)+G(n,nyp)+G(n,nym));   
        end
    end
end

V = G\B;

%map
Vmap = zeros(squares);
for i =1:1:squares
    for j =1:1:squares
        n=i+(j-1)*squares;
        Vmap(i,j) =V(n);
    end
end

[Ex,Ey] = gradient(Vmap*10^6);

forcex = -Ex*C.q_0;
forcey = -Ey*C.q_0;
accelerationX = forcex/mn;
accelerationY = forcey/mn;
accelVelocityX = accelerationX*(dt^2);
accelVelocityY = accelerationY*(dt^2);


figure(7)
boxplotX = [0.8e-7 0.8e-7 1.2e-7 1.2e-7];
boxplotY = [0 0.4e-7 0.4e-7 0];
plot(boxplotX,boxplotY,'color',[0 0 0]);
hold on
boxplotY = [1e-7 0.6e-7 0.6e-7 1e-7];
plot(boxplotX,boxplotY,'color',[0 0 0]);

for i = 1:1:steps
    
    %determine which accelerations to use
    for L = 1:1:squares
        for W = 1:1:squares
            axlogic = xr<L*xResolution & xr>(L-1)*xResolution;
            aylogic = yr<W*yResolution & yr>(W-1)*yResolution;
            
            xv(axlogic) = xv(axlogic)+ accelVelocityX(L,W);
            yv(aylogic) = yv(aylogic)+ accelVelocityY(L,W);
        end
    end

    scattered=rand(1,size);
    scatterCheck = scattered<=Pscat;
    velocity = Vthn*randn(1,size);
    xv(scatterCheck) = velocity(scatterCheck)*dt;
    velocity = Vthn*randn(1,size);
    yv(scatterCheck) = velocity(scatterCheck)*dt;
    tvelocity = sqrt((xv/dt).^2 +(yv/dt).^2);
    MFPcount(~scatterCheck) = MFPcount(~scatterCheck)+spacStep;
    
    box1left1ref = (xr + xv)>(xbox(1)-spacStep);
    box1right1ref= (xr + xv)<(xbox(2)+spacStep);
    box2bottom1ref  = (yr + yv)>yboxbottom(1) &(yr + yv)<yboxbottom(2);
    bottombox1 = box1left1ref & box1right1ref & box2bottom1ref;
    xv(bottombox1) = xv(bottombox1).*(-1);
    
    box1left2ref = (xr + xv)>xbox(1);
    box1right2ref= (xr + xv)<xbox(2);
    box2bottom2ref  = (yr + yv)>yboxbottom(1) &(yr + yv)<(yboxbottom(2)+spacStep);
    box2bottom = box1left2ref & box1right2ref & box2bottom2ref;
    yv(box2bottom) = yv(box2bottom).*(-1);
   
    box2top1ref  = (yr + yv)>yboxtop(1) &(yr + yv)<yboxtop(2);
    topbox1 = box1left1ref & box1right1ref & box2top1ref;
    xv(topbox1) = xv(topbox1).*(-1);
    
    box2top2ref  = (yr + yv)>(yboxtop(1)-spacStep) &(yr + yv)<yboxtop(2);
    topbox2 = box1left2ref & box1right2ref & box2top2ref;
    yv(topbox2) = yv(topbox2).*(-1);
    
 
    checkXright = xr +xv>2e-7;
    xr(checkXright) = xr(checkXright)+xv(checkXright)-l;
    checkXleft = xr +xv<0;
    xr(checkXleft) = xr(checkXleft) +xv(checkXleft)+l;
    
    leftover = ~(checkXright | checkXleft);
    
    xr(leftover) = xr(leftover) +xv(leftover);
    
    checkY = (yr+yv>1e-7 | yr+yv<0);
    yv(checkY) = yv(checkY).*(-1);
    yr(1,:) = yr(1,:)+yv(1,:);
    
    
    prevX(i,:) =xr(1,:);
    prevY(i,:) =yr(1,:);

end

for j = 1:1:numOfAtom
    plot(prevX(:,j),prevY(:,j),'color',[colour(1,j) 0 j/numOfAtom])

    xlim([0 l])
    ylim([0 w])
    hold on
    drawnow
end
title('Plot of trajectories'),xlabel('X'),ylabel('Y')

figure(8)
hist3([xr',yr'],[50,50]);
view(34,45)
title('Electron Density Map')

disp('For next step, the resolution of G matrix can be raised. ')







