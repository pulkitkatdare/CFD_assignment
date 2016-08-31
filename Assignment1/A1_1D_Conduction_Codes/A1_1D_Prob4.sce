clear
L=1.0; imax=12;xi=linspace(0,1,imax-1);
Beta=1.2;
Beta_p1=Beta+1;Beta_m1=Beta-1;Beta_p1_div_m1=(Beta_p1/Beta_m1)^(2*xi-1);
num=(Beta_p1*Beta_p1_div_m1)-Beta_m1;den=2*(1+Beta_p1_div_m1);
x=L*num./den;
y=x;//by symmetry of the situation  
xc(2:imax-1)=(x(2:imax-1)+x(1:imax-2))/2;
xc(1)=x(1);xc(imax)=x(imax-1);yc=xc
//grid parameters
Dx(2:imax-1)=x(2:imax-1)-x(1:imax-2);
dx(1:imax-1)=xc(2:imax)-xc(1:imax-1);
Dy = Dx; 
dy = dx;
printf("\n************** ONE-DIMENSIONAL HEAT CONDUCTION ***************");
//step 1: User inputs
rho = 7750.0; cp = 500.0; k = 16.2;
alpha=k/(rho*cp);
T0=30; T_wb=100.0;T_inf=30.0;h=100.0;qw = 10000.0
Q_vol_gen=0; epsilon_st=0.0001; Dt = 1000 //10e-7* (Dx(5)*Dx(5)/alpha);;
DTc = T_wb - T_inf
Q_gen=Q_vol_gen.*Dy'*Dx; // Total Heat Generation
aP0 = rho*cp*Dy'*Dx/Dt ;
aEW = k*Dy'*(1./dx)';// East and the west boundary
aNS = k*(1./dy)*Dx;// North and south boundary
T(2:imax-1,1) = T_wb;
T(2:imax-1,2:imax-1) = T0;
for j=2:imax-1
    for i =2:imax-1
        aP(j,i) = aP0(j,i) + aEW(j,i-1)+aEW(j,i) + aNS(j-1,i) + aNS(j,i);
    end
end

//add an imax term if necessary
unsteadiness_nd=1; n=0; alpha=k/(rho*cp); DTc=T_wb-T_inf;
while unsteadiness_nd>=epsilon_st
    n=n+1;
    T(imax,2:imax-1)=(k*T(imax-1,2:imax-1))+(h*dx(imax-1)*T_inf);
    T(imax,2:imax-1)=T(imax,2:imax-1)/(k+h*dx(imax-1));
    T(2:imax-1,imax)= T(2:imax-1,imax-1)  - (qw*dx(imax-1)/k);
    T(1,2:imax-1)= T(2,2:imax-1) //  - (qw*dx(imax-1)/k);
    T_old = T;
    for j=2:imax-1
        for i = 2:imax-1
            b(j,i) = aP0(j,i)*T_old(j,i) + Q_gen(j,i);
        end
    end
    epsilon = 0.0001;
    N = 0;
    Error = 1;
    while(Error >= epsilon)
        T_old_iter = T;
        N = N+1;
        for j = 2:imax-1
            for i =2:imax-1
            T(j,i) = aEW(j,i)*T(j,i+1)+ aEW(j,i-1)*T(j,i-1) + aNS(j,i)*T(j+1,i)+ aNS(j-1,i)*T(j-1,i)+ b(j,i)
            T(j,i)=T(j,i)/aP(j,i);
            end
        end
        Error=max(abs(T-T_old_iter));            
    end
    unsteadiness=max(abs(T-T_old))/Dt;
    unsteadiness_nd=unsteadiness*L*L/(alpha*DTc);
    printf("Time step no. %5d, Unsteadiness_nd = %8.4e\n", n , unsteadiness_nd);    
end
xset("colormap", jetcolormap(64)),colorbar(0,150),Sgrayplot(linspace(0,L,12),linspace(0,L,12),T, strf='100')
