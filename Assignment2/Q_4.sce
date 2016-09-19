clear 
//step 1: User inputs
rho = 7750; mu = 775;
L = 6; H = 1;
u = 1 ;
v = 0;
cp = 500;
k = cp*mu //Pr = 1;
imax = 62; jmax = 22
T_i = 50;
DTc = 100;
epsilon_st = 0.000001;
dx = L/(imax-2);
dy = L/(jmax-2);
Dt = 1e-3;
alpha=k/(rho*cp);
Beta=5;
xi=linspace(0,L,imax-1);
yi=linspace(0,H,jmax-1);
Beta_p1=Beta+1;Beta_m1=Beta-1;Beta_p1_div_m1=(Beta_p1/Beta_m1)^(2*xi-1);
num=(Beta_p1*Beta_p1_div_m1)-Beta_m1;den=2*(1+Beta_p1_div_m1);
x=L*num./den;
Beta_p1=Beta+1;Beta_m1=Beta-1;Beta_p1_div_m1=(Beta_p1/Beta_m1)^(2*yi-1);
num=(Beta_p1*Beta_p1_div_m1)-Beta_m1;den=2*(1+Beta_p1_div_m1);
y=H*num./den;
//x = xi;
//y = yi;
xc(2:imax-1)=(x(2:imax-1)+x(1:imax-2))/2;
xc(1)=x(1);xc(imax)=x(imax-1);
yc(2:jmax-1)=(y(2:jmax-1)+y(1:jmax-2))/2;
yc(1)=x(1);yc(jmax)=y(jmax-1);
//grid parameters
Dx(2:imax-1)=x(2:imax-1)-x(1:imax-2);
dx(1:imax-1)=xc(2:imax)-xc(1:imax-1);
Dy(2:jmax-1)=y(2:jmax-1)-y(1:jmax-2);
dy(1:jmax-1)=yc(2:jmax)-yc(1:jmax-1);
//Similarly to be written for y and x axis repectively
W_2(2:imax-1) = (2*Dx(1:imax-2)+Dx(2:imax-1))./(Dx(1:imax-2)+Dx(2:imax-1));
W_2(1) = 1
W_1(2:imax-1) = -Dx(2:imax-1)/(Dx(1:imax-2)+Dx(2:imax-1))
W_1(1) = 0
//step 2: Initialising the parameters
T = 0*ones(imax,jmax);
T(:,1) = 0;
T(:,jmax) = 0;
T(1,:) = 100;
unsteadiness_nd = 1;
N = 0;
while (unsteadiness_nd >= epsilon_st)
    N = N+1
    T_old = T;
    for i = 2:imax-1
        for j = 2:jmax-1
            if (i~=2 & j~=2)
            b(i,j) = rho*u*((W_2(i-1)-1)*T(i-1,j)+W_1(i-1)*T(i-2,j)-((W_2(i)-1)*T(i,j)+W_1(i)*T(i-1,j)))*Dy(j)
            end
        end
    end
    
    epsilon = 0.01;
    n = 0;
    Error = 1;
    while(Error >= epsilon)
        T_old_iter = T
        n = n+1
        for i = 2:imax-1
            for j = 2:jmax-1
                T(i,j) = T_old(i,j) + (Dt/(rho*cp*Dx(i)*Dy(j)))*(rho*cp*u*(T(i-1,j)-T(i,j))*Dy(j)+b(i,j)+k*((Dy(j)/dx(i-1))*(T(i-1,j)-T(i,j))-(Dy(j)/dx(i))*(T(i,j)-T(i+1,j))+(Dx(i)/dy(j-1))*(T(i,j-1)-T(i,j))-(Dx(i)/dy(j))*(T(i,j)-T(i,j+1))))
            end
        end
        T(imax,:) = T(imax-1,:)
        //T(:,jmax) = T(:,jmax-1)
        Error=max(abs(T-T_old_iter));
    end
    unsteadiness=max(abs(T-T_old))/DTc;
    unsteadiness_nd = unsteadiness*(L/(u*Dt));    
    printf("Time step no. %5d, Unsteadiness_nd = %8.4e\n", N , unsteadiness_nd);  
end
xset("colormap", jetcolormap(64)),colorbar(0,100),Sgrayplot(linspace(0,L,imax),linspace(0,H,jmax),T, strf='100')
