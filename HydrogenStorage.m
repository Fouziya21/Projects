m_in = 0.00025;  %mass inflow rate
m_out = 0.0002; %mass outflow rate
m_o=0.5;        %initial mass
T_i=300;        %infinity temperature
T_f=300;        %fluid temperature
T_0=300;        %initial temperature
a_f=24e-3;      %Heat Transfer Coeffcients
A_s=0.175929;   %Surface Area
c_v=10.316;     %Specific heat at constant volume
c_p=14.615;     %Specific heat at constant pressure
y=c_p/c_v;      %ratio of specific heats
R=8.314;        %Gas constant
V=0.1508;       %Volume of Tank
M_h2=2.0159*10^-3;  %Mass of Hydrogen 


m=zeros(8000,1);
Time=zeros(8000,1);
tou=zeros(8000,1);
alpha=zeros(8000,1);

T_char=zeros(8000,1);
T=zeros(8000,1);
P=zeros(8000,1);
H=zeros(8000,1);
Q=zeros(8000,1);
U=zeros(8000,1);

T_char_vinout=zeros(8000,1);
T_vinout=zeros(8000,1);
P_vinout=zeros(8000,1);
H_vinout=zeros(8000,1);
Q_vinout=zeros(8000,1);
U_vinout=zeros(8000,1);

T_char_cinvout=zeros(8000,1);
T_cinvout=zeros(8000,1);
P_cinvout=zeros(8000,1);
H_cinvout=zeros(8000,1);
Q_cinvout=zeros(8000,1);
U_cinvout=zeros(8000,1);

Filling_time=[40,190,370];
A_2=[1.204,1.195,1.173];
C_2=[1.631,1.000,1.026];
A_3=[-35.17,-9.504,-10.51];
B_3=[36.16,10.50,11.50];
C_3=[0.04595,0.05526,0.04841];
T_2=zeros(21,3);
T_3=zeros(21,3);
M=zeros(21,1);

for i=1:8000
    Time(i)=i;
    if i<=2000
        m(i)=m_o + m_in*i;
        alpha(i)=(a_f*A_s)/(c_v*m_in);
        tou(i)=(i*m_in)/m_o;

        T_char(i)=((y*T_i)+(alpha(i)*T_f))/(1+alpha(i));
        T(i)=T_char(i)-((T_char(i)-T_0)/((1+tou(i))^(1+alpha(i))));
        P(i)=(m(i)*R*T(i))/(M_h2*V);
        H(i)=m_in*c_p*T_i;
        Q(i)=a_f*A_s*(T_f-T(i));
        U(i)=H(i)+Q(i);

        T_char_vinout(i)=(alpha(i)*T_f)/(1+alpha(i)-y);
        T_vinout(i)=T_char_vinout(i)-((T_char_vinout(i)-T_0)/((1+tou(i))^(1+alpha(i)-y)));
        P_vinout(i)=(m(i)*R*T_vinout(i))/(M_h2*V);
        H_vinout(i)=m_in*c_p*T_vinout(i);
        Q_vinout(i)=a_f*A_s*(T_f-T_vinout(i));
        U_vinout(i)=H_vinout(i)+Q_vinout(i);

        T_char_cinvout(i)=((y*T_i)+(alpha(i)*T_f))/(1+alpha(i));
        T_cinvout(i)=T_char_cinvout(i)-((T_char_cinvout(i)-T_0)/((1+tou(i))^(1+alpha(i))));
        P_cinvout(i)=(m(i)*R*T_cinvout(i))/(M_h2*V);
        H_cinvout(i)=m_in*c_p*T_i;
        Q_cinvout(i)=a_f*A_s*(T_f-T_cinvout(i));
        U_cinvout(i)=H_cinvout(i)+Q_cinvout(i);


    elseif   i>2000 && i<=4000
        m(i)=m(i-1);
        tou(i)=((i-2000)*a_f*A_s)/(m(2000)*c_v);

         T_char(i)=T_f;
         T(i)=T_char(i) - ((T_char(i)-T(2000))*(exp(-tou(i))));
         P(i)=(m(i)*R*T(i))/(M_h2*V);
         H(i)=0;
         Q(i)=a_f*A_s*(T_f-T(i));
         U(i)=H(i)+Q(i);

         T_char_vinout(i)=T_f;
         T_vinout(i)=T_char_vinout(i) - ((T_char_vinout(i)-T_vinout(2000))*(exp(-tou(i))));
         P_vinout(i)=(m(i)*R*T_vinout(i))/(M_h2*V);
         H_vinout(i)=0;
         Q_vinout(i)=a_f*A_s*(T_f-T_vinout(i));
         U_vinout(i)=H_vinout(i)+Q_vinout(i);

         T_char_cinvout(i)=T_f;
         T_cinvout(i)=T_char_cinvout(i) - ((T_char_cinvout(i)-T_cinvout(2000))*(exp(-tou(i))));
         P_cinvout(i)=(m(i)*R*T_cinvout(i))/(M_h2*V);
         H_cinvout(i)=0;
         Q_cinvout(i)=a_f*A_s*(T_f-T_cinvout(i));
         U_cinvout(i)=H_cinvout(i)+Q_cinvout(i);


    elseif i>4000 && i<=6000
        m(i)=1- m_out*(i-4000);
        alpha(i)=(a_f*A_s)/(c_v*m_out*(-1));
        tou(i)=((i-4000)*m_out*(-1))/m(4000);

        T_char(i)=((y*T_i)+(alpha(i)*T_f))/(1+alpha(i));
        T(i)=T_char(i) - ((T_char(i)-T(4000))/((1+tou(i))^(1+alpha(i))));
        P(i)=(m(i)*R*T(i))/(M_h2*V);
        H(i)=-(m_out*T_i*c_p);
        Q(i)=a_f*A_s*(T_f-T(i));
        U(i)=H(i)+Q(i);

        T_char_vinout(i)=(alpha(i)*T_f)/(1+alpha(i)-y);
        T_vinout(i)=T_char_vinout(i)-((T_char_vinout(i)-T_vinout(4000))/((1+tou(i))^(1+alpha(i)-y)));
        P_vinout(i)=(m(i)*R*T_vinout(i))/(M_h2*V);
        H_vinout(i)=(-1)*m_out*c_p*T_vinout(i);
        Q_vinout(i)=a_f*A_s*(T_f-T_vinout(i));
        U_vinout(i)=H_vinout(i)+Q_vinout(i);

        T_char_cinvout(i)=((alpha(i)*T_f))/(1+alpha(i)-y);
        T_cinvout(i)=T_char_cinvout(i)-((T_char_cinvout(i)-T_cinvout(4000))/((1+tou(i))^(1+alpha(i)-y)));
        P_cinvout(i)=(m(i)*R*T_cinvout(i))/(M_h2*V);
        H_cinvout(i)=(-1)*m_out*c_p*T_cinvout(i);
        Q_cinvout(i)=a_f*A_s*(T_f-T_cinvout(i));
        U_cinvout(i)=H_cinvout(i)+Q_cinvout(i);

    else
        m(i)=m(i-1);
        tou(i)=((i-6000)*(A_s*a_f))/(m(6000)*c_v);

        T_char(i)=T_f;
        T(i)=T_char(i)-((T_char(i)-T(6000))*exp(-tou(i)));
        P(i)=(m(i)*R*T(i))/(M_h2*V);
        H(i)=0;
        Q(i)=a_f*A_s*(T_f-T(i));
        U(i)=H(i)+Q(i);

        T_char_vinout(i)=T_f;
        T_vinout(i)=T_char_vinout(i) - ((T_char_vinout(i)-T_vinout(6000))*(exp(-tou(i))));
        P_vinout(i)=(m(i)*R*T_vinout(i))/(M_h2*V);
        H_vinout(i)=0;
        Q_vinout(i)=a_f*A_s*(T_f-T_vinout(i));
        U_vinout(i)=H_vinout(i)+Q_vinout(i);

        T_char_cinvout(i)=T_f;
        T_cinvout(i)=T_char_cinvout(i) - ((T_char_cinvout(i)-T_cinvout(6000))*(exp(-tou(i))));
        P_cinvout(i)=(m(i)*R*T_cinvout(i))/(M_h2*V);
        H_cinvout(i)=0;
        Q_cinvout(i)=a_f*A_s*(T_f-T_cinvout(i));
        U_cinvout(i)=H_cinvout(i)+Q_cinvout(i);
    end
end

for j=1:3
    for i=1:21 
        M(i) = 0.75+ (0.25*i);
        T_2(i,j)=A_2(j) - (A_2(j)-1)*(M(i)^(-C_2(j)));
        T_3(i,j)=(A_3(j) +(B_3(j)*(M(i)^0.5)))^C_3(j);
    end
end

for i=1:21
    T_2_1(i) = T_2(i,1);
    T_2_2(i) = T_2(i,2);
    T_2_3(i) = T_2(i,3);
    T_3_1(i) = T_3(i,1);
    T_3_2(i) = T_3(i,2);
    T_3_3(i) = T_3(i,3);
end


hold on ;
plot(M,T_2_1,'Marker','+','Color','r','LineWidth',1.5);
plot(M,T_2_2,'Marker','diamond','Color','r','LineWidth',1.5);
plot(M,T_2_3,'Marker','hexagram','Color','r','LineWidth',1.5);
plot(M,T_3_1,'Marker','+','Color','b','LineWidth',1.5);
plot(M,T_3_2,'Marker','diamond','Color','b','LineWidth',1.5);
plot(M,T_3_3,'Marker','hexagram','Color','b','LineWidth',1.5);
grid on;
xlabel('Final/Inital Mass m/mo');
ylabel('Final/Initial Temperature T/To');
legend({'40s Fit','190s Fit','370s Fit','40s Data','190s Data','370sData'},'Location','southeast');
hold off;


figure;
plot(Time,U_cinvout);
xlabel('Time(secs)');
ylabel('Internal Energy(KJ)');
title('Constant inflow and Variable outflow temperatures');
grid on;

figure;
plot(Time,Q_cinvout);
xlabel('Time(secs)');
ylabel('Heat Flow(KJ)');
title('Constant inflow and Variable outflow temperatures');
grid on;

figure;
plot(Time,H_cinvout);
xlabel('Time(secs)');
ylabel('Enthalpy(KJ)');
title('Constant inflow and Variable outflow temperatures');
grid on;

figure;
plot(Time,P_cinvout);
xlabel('Time(secs)');
ylabel('Pressure(Pa)');
title('Constant inflow and Variable outflow temperatures');
grid on;


figure;
plot(Time,T_cinvout);
xlabel('Time(secs)');
ylabel('Temperature(K)');
title('Constant inflow and Variable outflow temperatures');
grid on;

figure;
plot(Time,m);
xlabel('Time(secs)');
ylabel('Mass(Kg)');
title('Constant inflow and Variable outflow temperatures');
grid on;

figure;
plot(Time,U_vinout);
xlabel('Time(secs)');
ylabel('Internal Energy(KJ)');
title('Variable inflow and outflow temperatures');
grid on;

figure;
plot(Time,Q_vinout);
xlabel('Time(secs)');
ylabel('Heat Flow(KJ)');
title('Variable inflow and outflow temperatures');
grid on;

figure;
plot(Time,H_vinout);
xlabel('Time(secs)');
ylabel('Enthalpy(KJ)');
title('Variable inflow and outflow temperatures');
grid on;


figure;
plot(Time,P_vinout);
xlabel('Time(secs)');
ylabel('Pressure(Pa)');
title('Variable inflow and outflow temperatures');
grid on;


figure;
plot(Time,T_vinout);
xlabel('Time(secs)');
ylabel('Temperature(K)');
title('Variable inflow and outflow temperatures');
grid on;

figure;
plot(Time,m);
xlabel('Time(secs)');
ylabel('Mass(Kg)');
title('Variable inflow and outflow temperatures');
grid on;


figure;
plot(Time,U);
xlabel('Time(secs)');
ylabel('Internal Energy(KJ)');
title('Constant inflow and outflow temperatures');
grid on;

figure;
plot(Time,Q);
xlabel('Time(secs)');
ylabel('Heat Flow(KJ)');
title('Constant inflow and outflow temperatures');
grid on;

figure;
plot(Time,H);
xlabel('Time(secs)');
ylabel('Enthalpy(KJ)');
title('Constant inflow and outflow temperatures');
grid on;

figure;
plot(Time,P);
xlabel('Time(secs)');
ylabel('Pressure(Pa)');
title('Constant inflow and outflow temperatures');
grid on;

figure;
plot(Time,T);
xlabel('Time(secs)');
ylabel('Temperature(K)');
title('Constant inflow and outflow temperatures');
grid on;


figure;
plot(Time,m);
xlabel('Time(secs)');
ylabel('Mass(Kg)');
title('Constant inflow and outflow temperatures');
grid on;
