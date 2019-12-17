%% Part 1: Finding Angle of Twist Based on Target Tip Speed Ratio and AOA
close all
clear all
clc

%----Given Conditions----%
U = 10; % m/s (Wind Speed)
R = (30/200); % m
rho = 1.23; %kg/m^3
kin_visc=1.516*10^(-5); %m^2/s
R_hub = (1.25/2)*(2.54/100); %m (Hub Radius)
mu_hub = R_hub/R;
w_hub = (1.5)*(2.54/100); %m (Hub Width)
c_min = (15/1000); %m (Minimum Cord at Tip)
c_max = (45/1000); %m (Max Cord)
res = (5/1000)*(2.54/100); %m (3D Printing Resolution)

%----Controllable Variables----%
n = 100; % Number of Blade Element Sections
dr = R/n; % m (Size of each radius step)
B = 3; % Number of Blades
c = [linspace(c_max,c_max,1*n/10),linspace(c_max,c_min,9*n/10)]; % m (Cord over the span)

%----Target Variables----%
alpha = [linspace(0,0,n/10),linspace(0,5,4*n/10),linspace(5,5,5*n/10)]; % degrees (Target Angle of Attack)
lamda = 5; % Target Tip Speed Ratio
omega = (lamda*U)/R; % rad/s (Angular Velocity)

%----Airfoil Data----%
% Columns: AOA , Cl , Cd
Data = csvread('xf-e214-il-100000.csv',11,0);
thick_cord_ratio = 0.111; % 11.1% Thickness at 33% Cord

%For use in post-stall Cl and Cd calculations (Viterna)
c_mean=mean(c); %Mean Cord
AR = R/c_mean; %Aspect Ratio
Cl_stall = 1; %Cl at Stall
Cd_stall = 0.3; %Cd at Stall
alpha_stall = 17; %degrees
Cd_max = 1.11+0.018*AR; %Max Cd at 90 degrees AOA
B1 = Cd_max;
B2 = (Cd_stall-Cd_max*(sind(alpha_stall)^2))/cosd(alpha_stall);
A1 = B1/2;
A2 = (Cl_stall-Cd_max*sind(alpha_stall)*cosd(alpha_stall))*sind(alpha_stall)/(cosd(alpha_stall)^2);


%----BEM Theory Calculation to Find Twist----%
%Creating Arrays for Variables that Vary with Radius
r = linspace(0,R*0.99,n); %radius array (0.99 of R else F-->1)
mu = r./R; %r/R
dCp = linspace(0,0,n); %Differential Power Coefficient
dT = linspace(0,0,n); %Differential Thrust
dM = linspace(0,0,n); %Differential Moment
dQ = linspace(0,0,n); %Differential Torque
dFt = linspace(0,0,n); %Differential Tangential Force
dP = linspace(0,0,n); %Differential Power
beta = linspace(0,0,n); %Angle of Twist
phi = linspace(0,0,n);
Axial_Ind_Factor = linspace(0,0,n); %Axial Induction Factor
Tangential_Ind_Factor = linspace(0,0,n); %Tangential Induction Factor
F = linspace(0,0,n); %Prantl's Tip Loss Factor
sigma = linspace(0,0,n); %Solidity
C_T = linspace(0,0,n);
Cl = linspace(0,0,n);
Cd = linspace(0,0,n);
Re = linspace(0,0,n);

%Looping from r=0 to r=R
for i = 1:n

    %Evaluating only after hub radius
    if r(i) < R_hub
    %Calculating Differential Values at r(i)
    dT(i) = 0;
    dM(i) = 0;
    dQ(i) = 0;
    dP(i) = 0;
    dCp(i) = 0;
    Axial_Ind_Factor(i) = 0;
    Tangential_Ind_Factor(i) = 0;

    else
        %Initialize a and a' and set tolerance
        a = 0;
        a_prime = 0;
        a_old = 0;
        a_prime_old = 0;
        target_tol = 0.000005;
        tolerance_a = 1; %To ensure first iteration occurs
        tolerance_a_prime = 1; %To ensure first iteration occurs

        %Calculate Local Solidity
        sigma(i) = (B*c(i))/(2*pi*r(i));

        %Finding a and a' loop ending when target tolerance is met
        while (tolerance_a > target_tol && tolerance_a_prime > target_tol)

            %Compute phi
            phi(i) = atand((U*(1-a))/(omega*r(i)*(1+a_prime)));

            %Ensuring Twist Angle is never negative
            if phi(i) < alpha(i)
                alpha(i)=phi(i);
            end

            %Compute beta
            beta(i) = phi(i)-alpha(i);

            %Coeff of Lift and Drag
            k = 1;
            %Checking whether or not alpha is within airfoil data range
            %If larger, use post stall theory (Viterna)
            if alpha(i) > Data(end,1)
                Cl(i) = A1*sind(2*alpha(i))+A2*(cosd(alpha(i))^2)/sind(alpha(i));
                Cd(i) = B1*(sind(alpha(i))^2)+B2*cosd(alpha(i));

            %If smaller, use lowest AOA data available
            elseif alpha(i) < Data(1,1)
                Cl(i) = Data(1,2);
                Cd(i) = Data(1,3);

            %Otherwise use the data in the airfoil csv file
            else
                while alpha(i) >= Data(k,1)
                    k=k+1;
                end
                %Interpolate
                Cl(i) = ((Data(k,2)-Data(k-1,2))/(Data(k,1)-Data(k-1,1)))*(alpha(i)-Data(k-1,1))+Data(k-1,2);
                Cd(i) = ((Data(k,3)-Data(k-1,3))/(Data(k,1)-Data(k-1,1)))*(alpha(i)-Data(k-1,1))+Data(k-1,3);
            end

%             %Test Cl and Cd Values
%             Cl(i) =0.8;
%             Cd(i) =0.02;

            %Find Cx and Cy (Force Coefficients)
            Cx = Cl(i)*cosd(phi(i))+Cd(i)*sind(phi(i));
            Cy = Cl(i)*sind(phi(i))-Cd(i)*cosd(phi(i));

            %Prantls Tip Loss Factor
            %With Losses
            f=(B/2)*(R-r(i))/(r(i)*sind(phi(i)));
            F(i)=(2/pi)*acos(exp(-f));

            %Check Local Thrust Coefficient
            C_T(i) = sigma(i)*(1-a)^2*Cx/(sind(phi(i))^2);

            %Save old and find new a and a_prime
            a_old = a;
            a_prime_old = a_prime;

%             %Accounting for Turbulent Wake State if C_T > 0.96
%             %Using Glauerts Empirical Relations
%             if C_T(i) > 0.96
%                 a = (1/F(i))*(0.143+(0.0203-0.6427*(0.889-C_T(i)))^0.5);
%
%             %Otherwise regular BEM equations used
%             else
%                 a = (sigma(i)*Cx)/(4*F(i)*sind(phi(i))^2+(sigma(i)*Cx));
%             end

            a = (sigma(i)*Cx)/(4*F(i)*sind(phi(i))^2+(sigma(i)*Cx));
            a_prime = (sigma(i)*Cy)/((4*F(i)*sind(phi(i))*cosd(phi(i)))-(sigma(i)*Cy));

            %Calculating Tolerance
            tolerance_a = a-a_old;
            tolerance_a_prime = a_prime-a_prime_old;
        end

    %Calculating Differential Values at r(i)
    dT(i) = 2*F(i)*rho*(U^2)*a*(1-a)*2*pi*r(i)*dr;
    dM(i) = dT(i)*r(i);
    dQ(i) = 2*F(i)*a_prime*(1-a)*rho*U*omega*(r(i)^2)*2*pi*r(i)*dr;
    dFt(i) = dQ(i)/r(i);
    dP(i) = omega*dQ(i);
    dCp(i) = dP(i)/((rho/2)*(U^3)*pi*(R^2));
    Axial_Ind_Factor(i) = a;
    Tangential_Ind_Factor(i) = a_prime;
    W(i) = ((U^2)*((1-a)^2)+(omega^2)*(r(i)^2)*(1+a_prime)^2)^0.5;
    Re(i) = (W(i)*c(i))/kin_visc;
    end

end


%----Calculating Totals----%
Cp_target = sum(dCp)
Thrust = sum(dT)
Torque = sum(dQ)
Power = sum(dP);

%----Plots----%
%Differential Power Coefficient along span
figure
plot(mu,dCp)
xlabel('r/R')
ylabel('Differential Power Coefficient')
ylim([0 inf])

%Induction Factors along span
figure
plot(mu,Axial_Ind_Factor,'LineWidth',2)
hold on
plot(mu,Tangential_Ind_Factor,'LineWidth',2)
legend('Axial Ind Factor','Tangential Ind Factor','Location','northwest')
xlabel('r/R')
ylabel('Induction Factor')
ylim([0 inf])
hold off

%Torque and Moment Plot
figure
plot(mu,dM,'LineWidth',2)
hold on
plot(mu,dQ,'LineWidth',2)
legend('Differential Moment','Differential Torque','Location','northwest')
xlabel('r/R')
ylabel('Torque [Nm]')
ylim([0 inf])
hold off

dM_dQ = dM./dQ;
%Moment to Torque Ratio
figure
plot(mu,dM_dQ,'LineWidth',2)
title('dM_dQ versus r/R')
hold off

%Angle of Attack and Twist Plot
figure
plot(mu,alpha,'LineWidth',2)
hold on
plot(mu,beta,'LineWidth',2)
legend('Angle of Attack','Angle of Twist','Location','northeast')
xlabel('r/R')
ylabel('Angle [Degrees]')
hold off

%Prantl's Correction Factor Plot
figure
plot(mu,F,'LineWidth',2)
xlabel('r/R')
ylabel('Prantl Tip Loss Factor')

%Phi
figure
plot(mu,phi,'LineWidth',2)
xlabel('r/R')
title('Phi versus r/R')
ylabel('Angle [Degrees]')

%Local Thrust Coefficient
figure
plot(mu,C_T,'LineWidth',2)
xlabel('r/R')
ylabel('Thrust Coefficient')

%Reynolds
figure
plot(mu,Re,'LineWidth',2)
xlabel('r/R')
ylabel('Reynolds Number')

%Reynolds
Ft_Fn = dFt./dT;
figure
plot(mu,Ft_Fn,'LineWidth',2)
xlabel('r/R')
ylabel('Ft/Fn')

%% Part 2: Use Determined Twist to Evaluate Performance at Various Tip Speed Ratios
%----Varying Tip Speed Ratios----%

Data = csvread('xf-e214-il-100000.csv',11,0);

U=10; %(Unhide to Analyze Performance starting - ie slower wind speeds)
lamda_start=1;
lamda_end=7;
n2 = lamda_end-lamda_start+1;
lamda = linspace(lamda_start,lamda_end,n2); % (Tip Speed Ratios)
%lamda = linspace(lamda_start,lamda_end,20); % (Tip Speed Ratios)
omega = (lamda.*U)./R; % rad/s (Angular Velocity)

%---BEM Theory Calculation for Varying Tip Speed Ratios----%
%Variables that Vary with Radius and Tip Speed Ratio
dCp = zeros(length(lamda),n); % Differential Power Coefficient
dT = zeros(length(lamda),n); % Differential Thrust
dM = zeros(length(lamda),n); % Differential Moment at Center
dM_root = zeros(length(lamda),n); % Differential Moment at root
dM_1_2span = zeros(length(lamda),n); % Differential Moment at 1/2 span
dM_3_4span = zeros(length(lamda),n); % Differential Moment at 3/4 span
dFt = zeros(length(lamda),n); % Differential Tangential Force
dQ = zeros(length(lamda),n); % Differential Torque
dP = zeros(length(lamda),n); % Differential Power
Axial_Ind_Factor = zeros(length(lamda),n); % Axial Induction Factor
Tangential_Ind_Factor = zeros(length(lamda),n); % Tangential Induction Factor
F = zeros(length(lamda),n); % Prantl's Tip Loss Factor
sigma = zeros(length(lamda),n); % Solidity
alpha = zeros(length(lamda),n); % degrees (Angle of Attack)
phi = zeros(length(lamda),n); % degrees (Angle of Attack)
Re = zeros(length(lamda),n); %Reynolds Number
W = zeros(length(lamda),n); %Relative Velocity
Cl = zeros(length(lamda),n); % Lift Coeff
Cd = zeros(length(lamda),n); % Drag Coeff

%Variables that Vary with Tip Speed Ratio
Cp_total = zeros(1,length(lamda));
Thrust = zeros(1,length(lamda));
Torque = zeros(1,length(lamda));
Power = zeros(1,length(lamda));
Moment = zeros(1,length(lamda)); %At the center

%Looping through the tip speed ratios
for j = 1:length(lamda)

    %Looping from r=0 to r=R
    for i = 1:n

        %Evaluating only after hub radius
        if r(i) < R_hub
            i_hub = i+1;
        else
            %Initialize a and a' and set tolerance
            a = 0;
            a_prime = 0;
            a_old = 0;
            a_prime_old = 0;
            target_tol = 0.000005;
            tolerance_a = 1; %To ensure first iteration occurs
            tolerance_a_prime = 1; %To ensure first iteration occurs

            %Calculate Solidity
            sigma(j,i) = (B*c(i))/(2*pi*r(i));

            %Finding a and a' loop ending when target tolerance is met
            while (tolerance_a > target_tol && tolerance_a_prime > target_tol)

                %Compute phi
                phi(j,i) = atand((U*(1-a))/(omega(j)*r(i)*(1+a_prime)));

                %Compute Angle of Attack
                alpha(j,i) = phi(j,i) - beta(i);

                %Coeff of Lift and Drag
                k = 1;
                %Checking whether or not it is within data range
                %If larger, use post stall theory
                if alpha(j,i) > Data(end,1)
                    Cl(j,i) = A1*sind(2*alpha(j,i))+A2*(cosd(alpha(j,i))^2)/sind(alpha(j,i));
                    Cd(j,i) = B1*(sind(alpha(j,i))^2)+B2*cosd(alpha(j,i));

                %If smaller, use lowest AOA data available
                elseif alpha(j,i) < Data(1,1)
                    Cl(j,i) = Data(1,2);
                    Cd(j,i) = Data(1,3);

                %Otherwise use the data in the csv files
                else
                    while alpha(j,i) >= Data(k,1)
                        k=k+1;
                    end

                    %Interpolate
                    Cl(j,i) = ((Data(k,2)-Data(k-1,2))/(Data(k,1)-Data(k-1,1)))*(alpha(j,i)-Data(k-1,1))+Data(k-1,2);
                    Cd(j,i) = ((Data(k,3)-Data(k-1,3))/(Data(k,1)-Data(k-1,1)))*(alpha(j,i)-Data(k-1,1))+Data(k-1,3);
                end

%                 %Test Cl and Cd Values
%                 Cl(j,i) =0.8;
%                 Cd(j,i) =0.02;

                %Find Cx and Cy
                Cx = Cl(j,i)*cosd(phi(j,i))+Cd(j,i)*sind(phi(j,i));
                Cy = Cl(j,i)*sind(phi(j,i))-Cd(j,i)*cosd(phi(j,i));

                %Prantls Tip Loss Factor
                %With Losses
                f=(B/2)*(R-r(i))/(r(i)*sind(phi(j,i)));
                F(j,i)=(2/pi)*acos(exp(-f));

                %Check Local Thrust Coefficient
                C_T = sigma(j,i)*((1-a)^2)*Cx/(sind(phi(j,i))^2);

                %Save old and find new a and a_prime
                a_old = a;
                a_prime_old = a_prime;

%                 %Accounting for Turbulent Wake State if C_T > 0.96
%                 %Using Glauerts Empirical Relations
%                 if C_T > 0.96
%                     a = (1/F(j,i))*(0.143+((0.0203-0.6427*(0.889-C_T))^0.5));
%                     a_prime = (sigma(j,i)*Cy)/((4*F(j,i)*sind(phi(j,i))*cosd(phi(j,i)))-(sigma(j,i)*Cy));
%
%                 %Otherwise regular BEM equations used
%                 else
%                     a = (sigma(j,i)*Cx)/((4*F(j,i)*(sind(phi(j,i))^2))+(sigma(j,i)*Cx));
%                     a_prime = (sigma(j,i)*Cy)/((4*F(j,i)*sind(phi(j,i))*cosd(phi(j,i)))-(sigma(j,i)*Cy));
%                 end

                a = (sigma(j,i)*Cx)/((4*F(j,i)*(sind(phi(j,i))^2))+(sigma(j,i)*Cx));
                a_prime = (sigma(j,i)*Cy)/((4*F(j,i)*sind(phi(j,i))*cosd(phi(j,i)))-(sigma(j,i)*Cy));

                %Calculating Tolerance
                tolerance_a = a-a_old;
                tolerance_a_prime = a_prime-a_prime_old;
            end

        %Calculating Differential Values at r(i)
        dT(j,i) = 2*F(j,i)*rho*(U^2)*a*(1-a)*2*pi*r(i)*dr;
        dM(j,i) = dT(j,i)*r(i);
        dQ(j,i) = 2*F(j,i)*a_prime*(1-a)*rho*U*omega(j)*(r(i)^2)*2*pi*r(i)*dr;
        dFt(j,i) = dQ(j,i)/r(i);
        dP(j,i) = omega(j)*dQ(j,i);
        dCp(j,i) = dP(j,i)/((rho/2)*(U^3)*pi*(R^2));
        Axial_Ind_Factor(j,i) = a;
        Tangential_Ind_Factor(j,i) = a_prime;
        W(j,i) = ((U^2)*((1-a)^2)+(omega(j)^2)*(r(i)^2)*(1+a_prime)^2)^0.5;
        Re(j,i) = (W(j,i)*c(i))/kin_visc;
        end

    end

    %----Calculating Totals----%
    Cp_total(j) = sum(dCp(j,:));
    Thrust(j) = sum(dT(j,:));
    Torque(j) = sum(dQ(j,:));
    Power(j) = sum(dP(j,:));
    Moment(j) = sum(dM(j,:));

end

Mn = zeros(length(lamda),n);
Mt = zeros(length(lamda),n);
M_total = zeros(length(lamda),n);
stress = zeros(length(lamda),n);
I_XX = linspace(0,0,n);
I_YY = linspace(0,0,n);
t = c.*thick_cord_ratio;
E = 3.1 * 10^9 ; % youngs modulus for abs plastic


for j = 1:length(lamda)
  for i = i_hub:n

      %Setting Values Back to Zero
    dMn = 0;
    Mn_temp = 0;
    dMt = 0;
    Mt_temp = 0;

    for k = i:n
        dMn = dT(j,k)*(r(k)-r(i));
        dMt = dFt(j,k)*(r(k)-r(i));
        Mn_temp = Mn_temp+dMn;
        Mt_temp = Mt_temp+dMt;
    end

    I_YY(i) = (1/4)*pi*((c(i)/2)^3)*t(i);
    I_XX(i) = (1/4)*pi*((t(i)/2)^3)*c(i);
    Mn(j,i) = Mn_temp;
    Mt(j,i) = Mt_temp;
    M_total(j,i) = sqrt((Mn(j,i)^2)+(Mt(j,i)^2));
    SF = 1; %Safety Factor
    stress(j,i) = (Mn(j,i)*c(i)/(2*I_XX(i))+Mt(j,i)*t(i)/(2*I_XX(i)))*SF;
  end
end

F_total = zeros(length(lamda),n);
Y_A = zeros(length(lamda),n);
% deflection loading scenario (Parabolic Loading)
for j = 1:length(lamda)
  for i = i_hub:n
      F_total(j,i) = M_total(j,i)/ r(i);
      Y_A(j,i) = ((F_total(j,i)*(R^4)/(12*E*I_XX(i)))*((-1*(mu(i).^6))/(18*(1+(r(i)-1)*mu(i)).^6)));
      Y_A2(j,i) = (F_total(j,i)*(R^4)/(12*E*I_XX(i)))*((((1/(r(i)-1)^5)*((3*(r(i)^4))-(27*(r(i)^3))-(47*(r(i)^2))+(13*r(i))-2))/(r(i)^3))+((60*log(r(i)))/(r(i)-1)));
  end
end


%F_tot = (sqrt((Fn(n)^2)+(Ft(n)^2)));
%Area = pi*w_chosen*c_chosen; % area of ellipse

%----Plots for Varying Tip Speed Ratio----%
%Differential Power Coefficient along span
figure
plot(mu,dCp,'LineWidth',2)
legend('TSR 1','TSR 2','TSR 3','TSR 4','TSR 5','TSR 6','TSR 7','TSR 8','TSR 9','TSR 10','Location','Northwest')
xlabel('r/R')
ylabel('Differential Power Coefficient')
ylim([-0.015 inf])
hold off

%Differential Power Coefficient along span
figure
plot(mu,dCp,'LineWidth',2)
legend('TSR 1','TSR 2','TSR 3','TSR 4','TSR 5','TSR 6','TSR 7','TSR 8','TSR 9','TSR 10','Location','Northwest')
xlabel('r/R')
ylabel('Differential Power Coefficient')
ylim([-0.015 inf])
hold off

%Axial Induction Factor along span
figure
%subplot(1,2,1)
plot(mu,M_total,'LineWidth',2)
legend('TSR 1','TSR 2','TSR 3','TSR 4','TSR 5','TSR 6','TSR 7','TSR 8','TSR 9','TSR 10','Location','Northwest')
xlabel('r/R')
ylabel('M_total')
ylim([0 1])

figure
%subplot(1,2,1)
plot(mu,F_total,'LineWidth',2)
legend('TSR 1','TSR 2','TSR 3','TSR 4','TSR 5','TSR 6','TSR 7','TSR 8','TSR 9','TSR 10','Location','Northwest')
xlabel('r/R')
ylabel('Ftotal')

figure
%subplot(1,2,1)
plot(mu,Y_A,'LineWidth',2)
legend('TSR 1','TSR 2','TSR 3','TSR 4','TSR 5','TSR 6','TSR 7','TSR 8','TSR 9','TSR 10','Location','Northwest')
xlabel('r/R')
ylabel('Y_A')

%Tangential Induction Factor along span
figure
%subplot(1,2,2)
plot(mu,Tangential_Ind_Factor,'LineWidth',2)
legend('TSR 1','TSR 2','TSR 3','TSR 4','TSR 5','TSR 6','TSR 7','TSR 8','TSR 9','TSR 10','Location','Northeast')
xlabel('r/R')
ylabel('Tangential Induction Factor')
ylim([0 0.6])
hold off

%Angle of Attack
figure
plot(mu,alpha,'LineWidth',2)
legend('TSR 1','TSR 2','TSR 3','TSR 4','TSR 5','TSR 6','TSR 7','TSR 8','TSR 9','TSR 10','Location','Northwest')
xlabel('r/R')
ylabel('Angle of Attack [Degrees]')
ylim([-15 inf])
hold off

%Prantl's Tip Loss Factor Plot
figure
plot(mu,F,'LineWidth',2)
legend('TSR 1','TSR 2','TSR 3','TSR 4','TSR 5','TSR 6','TSR 7','TSR 8','TSR 9','TSR 10','Location','Northwest')
xlabel('r/R')
ylabel('Prantls Tip Loss Factor')
hold off

%Power Coefficient for Various Tip Speed Ratios
figure
plot(lamda,Cp_total,'LineWidth',2)
xlabel('Tip Speed Ratio')
ylabel('Power Coefficient')
ylim([0 0.6])
hold off


dphi = zeros(1,length(phi));
for i=1:length(phi)
    dphi(i)=(phi(1,i)-phi(5,i));
end

%Change in Phi from TSR of 1 to 5
figure
plot(mu,dphi,'LineWidth',2)
xlabel('r/R')
ylabel('Change in Phi [degrees]')
hold off

%Renolds Number For Various TSR
figure
plot(mu,Re,'LineWidth',2)
xlabel('r/R')
ylabel('Reynolds Number')
ylim([0 100000])
legend('TSR 1','TSR 2','TSR 3','TSR 4','TSR 5','TSR 6','TSR 7','TSR 8','TSR 9','TSR 10','Location','Northwest')
hold off


%Stress
figure
plot(mu,stress,'LineWidth',2)
xlabel('r/R')
ylabel('Stress [Pa]')
ylim([0 inf])
legend('TSR 1','TSR 2','TSR 3','TSR 4','TSR 5','TSR 6','TSR 7','TSR 8','TSR 9','TSR 10','Location','Northeast')
hold off

%Prantl's Tip Loss Factor Plot
figure
plot(mu,dFt,'LineWidth',2)
legend('TSR 1','TSR 2','TSR 3','TSR 4','TSR 5','TSR 6','TSR 7','TSR 8','TSR 9','TSR 10','Location','Northwest')
xlabel('r/R')
ylabel('Differential Tangential Force [N]')
hold off

%Cord Length Variation
figure
plot(mu,c*1000,'LineWidth',2)
xlabel('r/R')
ylabel('Cord Length [mm]')

TABLE = [lamda',Cp_total']