clear all;
height=1; % max height of the net 
g=9.81;
Cd=0.65;
m=0.058;
ro=1.2;
d=0.065;
f=40;
v0=60;
alfa=pi/2+pi/42
beta=pi/10
% v0=13;
% alfa=pi/6+pi/7
% beta=pi/10
x(1)=-1;
y(1)=2.5;
z(1)=2.5;
time(1)=0;
Vx(1)=v0*sin(alfa);
Vy(1)=v0*sin(beta)*sin(alfa);
Vz(1)=v0*cos(alfa);

i=1;
c=10.e-3;
h=0.001; %step

%net function (parabolic)
M=[10.97^2 10.97 1; (10.97/2)^2 10.97/2 1; 0 0 1];
N=[height; 0.914; height];
coefficients=linsolve(M,N);
net_height=@(y) coefficients(1)*y.^2+coefficients(2)*y+coefficients(3);

%drift 
Fdx=@(Vx,Vy,Vz) Cd*ro*pi*(d^2/8)*sqrt(Vx^2+Vy^2+Vz^2)*(-Vx);
Fdy=@(Vx,Vy,Vz) Cd*ro*pi*(d^2/8)*sqrt(Vx^2+Vy^2+Vz^2)*(-Vy);
Fdz=@(Vx,Vy,Vz) Cd*ro*pi*(d^2/8)*sqrt(Vx^2+Vy^2+Vz^2)*(-Vz);

%magnus effect
Ix=1; 
Iy=1;
Iz=1;
Cm=@(Vx,Vy,Vz) 1/(2+0.98*sqrt(Vx^2+Vy^2+Vz^2)/(pi*d*f));
Fmx=@(Vx,Vy,Vz) Cm(Vx,Vy,Vz)*ro*pi*(d^2/8)*((Vx^2+Vy^2+Vz^2)*(Iy*Vz-Iz*Vy))/sqrt((Iy*Vz-Iz*Vy)^2+(-Ix*Vz+Iz*Vx)^2+(Ix*Vy-Iy*Vx)^2);
Fmy=@(Vx,Vy,Vz) Cm(Vx,Vy,Vz)*ro*pi*(d^2/8)*((Vx^2+Vy^2+Vz^2)*(-Ix*Vz+Iz*Vx))/sqrt((Iy*Vz-Iz*Vy)^2+(-Ix*Vz+Iz*Vx)^2+(Ix*Vy-Iy*Vx)^2);
Fmz=@(Vx,Vy,Vz) Cm(Vx,Vy,Vz)*ro*pi*(d^2/8)*((Vx^2+Vy^2+Vz^2)*(Ix*Vy-Iy*Vx))/sqrt((Iy*Vz-Iz*Vy)^2+(-Ix*Vz+Iz*Vx)^2+(Ix*Vy-Iy*Vx)^2);


%equations with all forces applied 
f1=@(time,x,y,z,Vx,Vy,Vz)  Vx; % define x
f2=@(time,x,y,z,Vx,Vy, Vz) Vy; % define y
f3=@(time,x,y,z,Vx,Vy, Vz) Vz; % define z
L1=@(time,x,y,z,Vx,Vy, Vz) Fdx(Vx, Vy, Vz)/m+Fmx(Vx, Vy, Vz)/m; % define Vx
L2=@(time,x,y,z,Vx,Vy, Vz) Fdy(Vx, Vy, Vz)/m+Fmy(Vx, Vy, Vz)/m; % define Vy
L3=@(time,x,y,z,Vx,Vy, Vz) Fdz(Vx, Vy, Vz)/m+Fmz(Vx, Vy, Vz)/m-g; % define Vz

while x(i)<=23.78 && y(i)>0 && y(i)<10.97 && z(i)>=0 
    index=i;
    a=i;
    b=i+1;
    quit=0;
    quitRK = 0;
    for i=a:b % RK method 4
        if x(i)>=23.78 || y(i)<0 || y(i)>10.97 || z(i)<0 
            quitRK = 1;
            break
        end
        time(i+1)=time(i)+h;
        
        k11=x(i);
        k12=y(i);
        k13=z(i);
        k14=Vx(i);
        k15=Vy(i);
        k16=Vz(i);
        
        k21=x(i)+(h/2)*f1(time(i),k11,k12,k13,k14,k15,k16);
        k22=y(i)+(h/2)*f2(time(i),k11,k12,k13,k14,k15,k16);
        k23=z(i)+(h/2)*f3(time(i),k11,k12,k13,k14,k15,k16);
        k24=Vx(i)+(h/2)*L1(time(i),k11,k12,k13,k14,k15,k16);
        k25=Vy(i)+(h/2)*L2(time(i),k11,k12,k13,k14,k15,k16);
        k26=Vz(i)+(h/2)*L3(time(i),k11,k12,k13,k14,k15,k16);
        
        k31=x(i)+(h/2)*f1(time(i)+(h/2),k21,k22,k23,k24,k25,k26);
        k32=y(i)+(h/2)*f2(time(i)+(h/2),k21,k22,k23,k24,k25,k26);
        k33=z(i)+(h/2)*f3(time(i)+(h/2),k21,k22,k23,k24,k25,k26);
        k34=Vx(i)+(h/2)*L1(time(i)+(h/2),k21,k22,k23,k24,k25,k26);
        k35=Vy(i)+(h/2)*L2(time(i)+(h/2),k21,k22,k23,k24,k25,k26);
        k36=Vz(i)+(h/2)*L3(time(i)+(h/2),k21,k22,k23,k24,k25,k26);
        
        k41=x(i)+h*f1(time(i)+(h/2),k31,k32,k33,k34,k35,k36);
        k42=y(i)+h*f2(time(i)+(h/2),k31,k32,k33,k34,k35,k36);
        k43=z(i)+h*f3(time(i)+(h/2),k31,k32,k33,k34,k35,k36);
        k44=Vx(i)+h*L1(time(i)+(h/2),k31,k32,k33,k34,k35,k36);
        k45=Vy(i)+h*L2(time(i)+(h/2),k31,k32,k33,k34,k35,k36);
        k46=Vz(i)+h*L3(time(i)+(h/2),k31,k32,k33,k34,k35,k36);
        
        x(i+1)=x(i)+(h/6)*(f1(time(i),k11,k12,k13,k14,k15,k16)+2*f1(time(i)+(h/2),k21,k22,k23,k24,k25,k26)+2*f1(time(i)+(h/2),k31,k32,k33,k34,k35,k36)+f1(time(i)+h,k41,k42,k43,k44,k45,k46));
        y(i+1)=y(i)+(h/6)*(f2(time(i),k11,k12,k13,k14,k15,k16)+2*f2(time(i)+(h/2),k21,k22,k23,k24,k25,k26)+2*f2(time(i)+(h/2),k31,k32,k33,k34,k35,k36)+f2(time(i)+h,k41,k42,k43,k44,k45,k46));
        z(i+1)=z(i)+(h/6)*(f3(time(i),k11,k12,k13,k14,k15,k16)+2*f3(time(i)+(h/2),k21,k22,k23,k24,k25,k26)+2*f3(time(i)+(h/2),k31,k32,k33,k34,k35,k36)+f3(time(i)+h,k41,k42,k43,k44,k45,k46));
        Vx(i+1)=Vx(i)+(h/6)*(L1(time(i),k11,k12,k13,k14,k15,k16)+2*L1(time(i)+(h/2),k21,k22,k23,k24,k25,k26)+2*L1(time(i)+(h/2),k31,k32,k33,k34,k35,k36)+L1(time(i)+h,k41,k42,k43,k44,k45,k46));
        Vy(i+1)=Vy(i)+(h/6)*(L2(time(i),k11,k12,k13,k14,k15,k16)+2*L2(time(i)+(h/2),k21,k22,k23,k24,k25,k26)+2*L2(time(i)+(h/2),k31,k32,k33,k34,k35,k36)+L2(time(i)+h,k41,k42,k43,k44,k45,k46));
        Vz(i+1)=Vz(i)+(h/6)*(L3(time(i),k11,k12,k13,k14,k15,k16)+2*L3(time(i)+(h/2),k21,k22,k23,k24,k25,k26)+2*L3(time(i)+(h/2),k31,k32,k33,k34,k35,k36)+L3(time(i)+h,k41,k42,k43,k44,k45,k46));
    
    end
    
    i=index+2;
    
    while z(i)>=0 
        if x(i)>=23.78 || y(i)<0 || y(i)>10.97
            quit = 1;
            break
        end
        
        if x(i)>= 11.89 && x(i-1)<11.89 % checking if the ball can fly above the net
            x_net=11.89;
            y(i)=((x_net-x(i-1))*(y(i)-y(i-1))/(x(i)-x(i-1)))+y(i-1);
            z(i)=((x_net-x(i-1))*(z(i)-z(i-1))/(x(i)-x(i-1)))+z(i-1);
            Vx(i)=((x_net-x(i-1))*(Vx(i)-Vx(i-1))/(x(i)-x(i-1)))+Vx(i-1);
            Vy(i)=((x_net-x(i-1))*(Vy(i)-Vy(i-1))/(x(i)-x(i-1)))+Vy(i-1);
            Vz(i)=(((x_net-x(i-1))*(Vz(i)-Vz(i-1))/(x(i)-x(i-1)))+Vz(i-1));
            x(i)=x_net;
            
            if z(i)<= net_height(y(i))
                quit=1;
                break
            end
        end
        time(i+1)=time(i)+h;
        
        x_star=x(i)+h*f1(time(i),x(i),y(i),z(i),Vx(i),Vy(i),Vz(i));
        y_star=y(i)+h*f2(time(i),x(i),y(i),z(i),Vx(i),Vy(i),Vz(i));
        z_star=z(i)+h*f3(time(i),x(i),y(i),z(i),Vx(i),Vy(i),Vz(i));
        Vx_star=Vx(i)+h*L1(time(i),x(i),y(i),z(i),Vx(i),Vy(i),Vz(i));
        Vy_star=Vy(i)+h*L2(time(i),x(i),y(i),z(i),Vx(i),Vy(i),Vz(i));
        Vz_star=Vz(i)+h*L3(time(i),x(i),y(i),z(i),Vx(i),Vy(i),Vz(i));
        
        x(i+1)=x(i)+(h/24)*(9*f1(time(i+1),x_star,y_star,z_star,Vx_star,Vy_star,Vz_star)+19*f1(time(i),x(i),y(i),z(i),Vx(i),Vy(i),Vz(i))-5*f1(time(i-1),x(i-1),y(i-1),z(i-1),Vx(i-1),Vy(i-1),Vz(i-1))+f1(time(i-2),x(i-2),y(i-2),z(i-2),Vx(i-2),Vy(i-2),Vz(i-2)));
        y(i+1)=y(i)+(h/24)*(9*f2(time(i+1),x_star,y_star,z_star,Vx_star,Vy_star,Vz_star)+19*f2(time(i),x(i),y(i),z(i),Vx(i),Vy(i),Vz(i))-5*f2(time(i-1),x(i-1),y(i-1),z(i-1),Vx(i-1),Vy(i-1),Vz(i-1))+f2(time(i-2),x(i-2),y(i-2),z(i-2),Vx(i-2),Vy(i-2),Vz(i-2)));
        z(i+1)=z(i)+(h/24)*(9*f3(time(i+1),x_star,y_star,z_star,Vx_star,Vy_star,Vz_star)+19*f3(time(i),x(i),y(i),z(i),Vx(i),Vy(i),Vz(i))-5*f3(time(i-1),x(i-1),y(i-1),z(i-1),Vx(i-1),Vy(i-1),Vz(i-1))+f3(time(i-2),x(i-2),y(i-2),z(i-2),Vx(i-2),Vy(i-2),Vz(i-2)));
        Vx(i+1)=Vx(i)+(h/24)*(9*L1(time(i+1),x_star,y_star,z_star,Vx_star,Vy_star,Vz_star)+19*L1(time(i),x(i),y(i),z(i),Vx(i),Vy(i),Vz(i))-5*L1(time(i-1),x(i-1),y(i-1),z(i-1),Vx(i-1),Vy(i-1),Vz(i-1))+L1(time(i-2),x(i-2),y(i-2),z(i-2),Vx(i-2),Vy(i-2),Vz(i-2)));
        Vy(i+1)=Vy(i)+(h/24)*(9*L2(time(i+1),x_star,y_star,z_star,Vx_star,Vy_star,Vz_star)+19*L2(time(i),x(i),y(i),z(i),Vx(i),Vy(i),Vz(i))-5*L2(time(i-1),x(i-1),y(i-1),z(i-1),Vx(i-1),Vy(i-1),Vz(i-1))+L2(time(i-2),x(i-2),y(i-2),z(i-2),Vx(i-2),Vy(i-2),Vz(i-2)));
        Vz(i+1)=Vz(i)+(h/24)*(9*L3(time(i+1),x_star,y_star,z_star,Vx_star,Vy_star,Vz_star)+19*L3(time(i),x(i),y(i),z(i),Vx(i),Vy(i),Vz(i))-5*L3(time(i-1),x(i-1),y(i-1),z(i-1),Vx(i-1),Vy(i-1),Vz(i-1))+L3(time(i-2),x(i-2),y(i-2),z(i-2),Vx(i-2),Vy(i-2),Vz(i-2)));
        
        expression1=abs((x(i+1)-x_star)/x_star);
        expression2=abs((y(i+1)-y_star)/y_star);
        expression3=abs((z(i+1)-z_star)/z_star);
        expression4=abs((Vx(i+1)-Vx_star)/Vx_star);
        expression5=abs((Vy(i+1)-Vy_star)/Vy_star);
        expression6=abs((Vz(i+1)-Vz_star)/Vz_star);
        
        while (expression1>c) || (expression2>c) || (expression3>c) || (expression4>c) || (expression5>c) || (expression6>c)
            x_star=x(i+1);
            y_star=y(i+1);
            z_star=z(i+1);
            Vx_star=Vx(i+1);
            Vy_star=Vy(i+1);
            Vz_star=Vz(i+1);
            
            x(i+1)=x(i)+(h/24)*(9*f1(time(i+1),x_star,y_star,z_star,Vx_star,Vy_star,Vz_star)+19*f1(time(i),x(i),y(i),z(i),Vx(i),Vy(i),Vz(i))-5*f1(time(i-1),x(i-1),y(i-1),z(i-1),Vx(i-1),Vy(i-1),Vz(i-1))+f1(time(i-2),x(i-2),y(i-2),z(i-2),Vx(i-2),Vy(i-2),Vz(i-2)));
            y(i+1)=y(i)+(h/24)*(9*f2(time(i+1),x_star,y_star,z_star,Vx_star,Vy_star,Vz_star)+19*f2(time(i),x(i),y(i),z(i),Vx(i),Vy(i),Vz(i))-5*f2(time(i-1),x(i-1),y(i-1),z(i-1),Vx(i-1),Vy(i-1),Vz(i-1))+f2(time(i-2),x(i-2),y(i-2),z(i-2),Vx(i-2),Vy(i-2),Vz(i-2)));
            z(i+1)=z(i)+(h/24)*(9*f3(time(i+1),x_star,y_star,z_star,Vx_star,Vy_star,Vz_star)+19*f3(time(i),x(i),y(i),z(i),Vx(i),Vy(i),Vz(i))-5*f3(time(i-1),x(i-1),y(i-1),z(i-1),Vx(i-1),Vy(i-1),Vz(i-1))+f3(time(i-2),x(i-2),y(i-2),z(i-2),Vx(i-2),Vy(i-2),Vz(i-2)));
            Vx(i+1)=Vx(i)+(h/24)*(9*L1(time(i+1),x_star,y_star,z_star,Vx_star,Vy_star,Vz_star)+19*L1(time(i),x(i),y(i),z(i),Vx(i),Vy(i),Vz(i))-5*L1(time(i-1),x(i-1),y(i-1),z(i-1),Vx(i-1),Vy(i-1),Vz(i-1))+L1(time(i-2),x(i-2),y(i-2),z(i-2),Vx(i-2),Vy(i-2),Vz(i-2)));
            Vy(i+1)=Vy(i)+(h/24)*(9*L2(time(i+1),x_star,y_star,z_star,Vx_star,Vy_star,Vz_star)+19*L2(time(i),x(i),y(i),z(i),Vx(i),Vy(i),Vz(i))-5*L2(time(i-1),x(i-1),y(i-1),z(i-1),Vx(i-1),Vy(i-1),Vz(i-1))+L2(time(i-2),x(i-2),y(i-2),z(i-2),Vx(i-2),Vy(i-2),Vz(i-2)));
            Vz(i+1)=Vz(i)+(h/24)*(9*L3(time(i+1),x_star,y_star,z_star,Vx_star,Vy_star,Vz_star)+19*L3(time(i),x(i),y(i),z(i),Vx(i),Vy(i),Vz(i))-5*L3(time(i-1),x(i-1),y(i-1),z(i-1),Vx(i-1),Vy(i-1),Vz(i-1))+L3(time(i-2),x(i-2),y(i-2),z(i-2),Vx(i-2),Vy(i-2),Vz(i-2)));
            
            expression1=abs((x(i+1)-x_star)/x_star);
            expression2=abs((y(i+1)-y_star)/y_star);
            expression3=abs((z(i+1)-z_star)/z_star);
            expression4=abs((Vx(i+1)-Vx_star)/Vx_star);
            expression5=abs((Vy(i+1)-Vy_star)/Vy_star);
            expression6=abs((Vz(i+1)-Vz_star)/Vz_star);
        end
        i=i+1;
    end
    
    if quit == 1
        break
    end
    z0=0;
    x(i)=((z0-z(i-1))*(x(i)-x(i-1))/(z(i)-z(i-1)))+x(i-1);
    y(i)=((z0-z(i-1))*(y(i)-y(i-1))/(z(i)-z(i-1)))+y(i-1);
    Vx(i)=((z0-z(i-1))*(Vx(i)-Vx(i-1))/(z(i)-z(i-1)))+Vx(i-1);
    Vy(i)=((z0-z(i-1))*(Vy(i)-Vy(i-1))/(z(i)-z(i-1)))+Vy(i-1);
    Vz(i)=-(((z0-z(i-1))*(Vz(i)-Vz(i-1))/(z(i)-z(i-1)))+Vz(i-1));
    z(i)=z0;
    if quit == 1
        break
    end
end

A1=[0 0 0];
A2=[0 10.97 0];
pts=[A1;A2];
plot3(pts(:,1), pts(:,2), pts(:,3),'b');
hold on;
B1=[0 0 0];
B2=[23.78 0 0];
pts=[B1;B2];
plot3(pts(:,1), pts(:,2), pts(:,3),'b');
hold on;
A11=[23.78 0 0];
A22=[23.78 10.97 0];
pts=[A11;A22];
plot3(pts(:,1), pts(:,2), pts(:,3),'b');
hold on;
B11=[0 10.97 0];
B22=[23.78 10.97 0];
pts=[B11;B22];
plot3(pts(:,1), pts(:,2), pts(:,3),'b');
hold on;
C1=[0 1.37 0];
C2=[23.78 1.37 0];
pts=[C1;C2];
plot3(pts(:,1), pts(:,2), pts(:,3),'b');
hold on;
C11=[0 9.6 0];
C22=[23.78 9.6 0];
pts=[C11;C22];
plot3(pts(:,1), pts(:,2), pts(:,3),'b');
D1=[5.49 1.37 0];
D2=[5.49 9.6 0];
pts=[D1;D2];
plot3(pts(:,1), pts(:,2), pts(:,3),'b');
hold on;
D11=[18.29 1.37 0];
D22=[18.29 9.6 0];
pts=[D11;D22];
plot3(pts(:,1), pts(:,2), pts(:,3),'b');
hold on;
E1=[5.49 5.485 0];
E2=[18.29 5.485 0];
pts=[E1;E2];
plot3(pts(:,1), pts(:,2), pts(:,3),'b');
hold on;
E11=[11.89 0 0];
E22=[11.89 10.97 0];
pts=[E11;E22];
plot3(pts(:,1), pts(:,2), pts(:,3),'r');
hold on;
y_net=0:h:10.97;
x_net=11.89.*ones(length(y_net),1);
plot3(x_net, y_net, net_height(y_net),'r');
hold on;
G1=[11.89 0 0];
G2=[11.89 0 1];
pts=[G1;G2];
plot3(pts(:,1), pts(:,2), pts(:,3),'r');
hold on;
G11=[11.89 10.97 0];
G22=[11.89 10.97 1];
pts=[G11;G22];
plot3(pts(:,1), pts(:,2), pts(:,3),'r');
hold on;


plot3(x,y,z, 'm');
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
grid on;


