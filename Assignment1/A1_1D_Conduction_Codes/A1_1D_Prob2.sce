clear
printf("\n*********************2D Heat Conduction***********\n");
//user input 
rho = 7750.0; cp = 500.0; k = 16.2;
Lx=1.0; imax=12; jmax = imax;
T0 = 30; Twb = 100; Tsb = 200; Teb = 300; Tnb = 400;
Q_vol_gen=0; epsilon_st=0.0001; 
alpha = k/(rho*cp);
DTc = Tnb - Twb; //maximum temperature difference 
Dx = Lx/(imax-1); Dy = Dx; Dt =  0.99*0.25*Dx*Dx/(alpha); //grid paramters defined 
Da = Dx*Dy //area of the square plate
Q_gen = Q_vol_gen*Da;
//ICs and BCs defined here 
T(1, 2:imax-1) = Twb;
T(jmax, 2:imax-1) = Teb;
for j=2:jmax-2
   T(j, 1) = Tsb; 
   T(j, imax) = Tnb;
end
for i =2:imax-1
    for j = 2:imax-1
        T(j,i) = T0;
    end
end
for j=1:jmax-1,
        X(j,1:imax-1) = linspace(0, Lx, imax-1);
    for i=1:imax-1,
        Y(j,i) = (j-1)*Dy;
    end
end
Tx(1,1) = X(1,1); Tx(1, imax) = X(1, imax-1); Tx(jmax, 1) = X(jmax-1, 1); Tx(jmax, imax) = X(jmax-1, imax-1);
Ty(1,1) = Y(1,1); Ty(1, imax) = Y(1, imax-1); Ty(jmax, 1) = Y(jmax-1, 1); Ty(jmax, imax) = Y(jmax-1, imax-1);
for j=2:jmax-1,
    for i=2:imax-1,
        Tx(j, i) = (X(j,i)+X(j,i-1))/2;
        Ty(j, i) = (Y(j,i)+Y(j-1,i))/2;
    end
end
for i=2:imax-1
    Tx(1, i) = (X(1,i)+X(1,i-1))/2;
    Ty(1, i) = Y(1, 1);
    Tx(jmax, i) = (X(jmax-1,i)+X(jmax-1,i-1))/2;
    Ty(jmax, i) = Y(jmax-1, 1);
end
for j=2:jmax-1
    Ty(j, 1) = (Y(j,1)+Y(j-1,1))/2;
    Tx(j, 1) = X(1, 1);
    Ty(j, imax) = (Y(j,imax-1)+Y(j-1,imax-1))/2;
    Tx(j, imax) = X(1, imax-1);
end
figure(1)
plot(X,Y,'m-s');
xlabel('X'); ylabel('Y');
plot(Tx,Ty,'k-o');
unsteadiness_nd = 1;
n = 0 ;
T_old1 = T
while (unsteadiness_nd >= epsilon_st)//note that the y coordinate is always written first in a matrix notation
T_old = T
n = n + 1 ;
for j = 2:imax-1
    for i = 1:imax-1
        if (i == 1) | (i == imax-1)
            qx(j,i) = -k*(T_old(j,i+1) - T_old(j,i))/(Dx/2.0); 
        else 
            qx(j,i) = -k*(T_old(j,i+1) - T_old(j,i))/(Dx);
        end
    end     
end
for j = 1:imax-1
    for i = 2:imax-1
        if (j == 1) | (j == imax-1)
            qy(j,i) = -k*(T_old(j+1,i) - T_old(j,i))/(Dx/2.0); 
        else 
            qy(j,i) = -k*(T_old(j+1,i) - T_old(j,i))/(Dx);
        end
    end     
end
for j = 2:imax-1
    for i = 2:imax-1   
        Q_cond(j,i) = (qx(j,i-1) - qx(j,i))*Dy + (qy(j-1,i) - qy(j,i))*Dx 
        T(j,i) = T_old(j,i) + (Dt/(rho*cp*Da))*(Q_cond(j,i)+Q_gen);
    end
end
unsteadiness=max(abs(T-T_old))/Dt;
unsteadiness_nd=unsteadiness*Lx*Lx/(alpha*DTc); //STEP 5: Steady state converegnce criterion
printf("Time step no. %5d, Unsteadiness) = %8.4e\n", n , unsteadiness);
end
//
//figure(2)
//title('TEMPERATURE PROFILE PLOT','color','black','fontsize',3);
//Sgrayplot(linspace(0,Lx,12),linspace(0,Lx,12),T, strf='090')
xset("colormap", jetcolormap(64)),colorbar(0,150),Sgrayplot(linspace(0,Lx,12),linspace(0,Lx,12),T, strf='100')


