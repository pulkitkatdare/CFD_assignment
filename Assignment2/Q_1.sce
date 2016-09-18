clear 
//step 1: User inputs
rho = 1000; 
L = 1; H = 1;
u = 1 ;
v = 1;
cp = 4180;
imax = 32; jmax = 32
T_i = 50;
DTc = 100;
epsilon_st = 0.000001;
dx = L/(imax-2);
dy = L/(jmax-2);
Dt = 0.001;
//step 2: Initialising the parameters
grid = (L/(imax-1))*ones(imax-1,imax-1);//distance between ith and (i+1)th grid
grid(1,:) = 0.5*(L/(imax-1));//in grid(i,j) 'i' corresponds to x axis, and vice versa
grid(imax-1,:) = 0.5*(L/(imax-1));
grid(:,1) = 0.5*(L/(imax-1));
grid(:,jmax-1) = 0.5*(L/(imax-1));
T = 50*ones(imax,jmax);
T(:,1) = 0;
T(1,:) = 100;
unsteadiness_nd = 1;
N = 0;
while(unsteadiness_nd >= epsilon_st)
   N = N+1;
   T_old = T;
   for i = 1:imax
       for j = 1:jmax
       h_x(i,j) = T(i,j)*cp*rho*u*dy;
       h_y(i,j) = T(i,j)*cp*rho*v*dx;
       end
   end
   for i = 2:imax
       for j = 2:imax
           if (i == 2 & j ~=2 & j ~= imax)    
             T(i,j) = T_old(i,j) + (Dt/(rho*cp*dx*dy))*((0.375*h_x(i,j)+0.75*h_x(i-1,j)-0.125*h_x(i-1,j))-(0.375*h_x(i+1,j)+0.75*h_x(i,j)-0.125*h_x(i-1,j)) + (0.375*h_y(i,j)+0.75*h_y(i,j-1)-0.125*h_y(i,j-2)) - (0.375*h_y(i,j+1)+0.75*h_y(i,j)-0.125*h_y(i,j-1)))  
           //T(i,j) = T_old(i,j) + (Dt/(rho*cp*dx*dy))*(h_x(i-1,j)-h_x(i,j) + h_y(i,j-1) - h_y(i,j))//FOU
       elseif (i == 2 & j ==2)
           T(i,j) = T_old(i,j) + (Dt/(rho*cp*dx*dy))*((0.375*h_x(i,j)+0.75*h_x(i-1,j)-0.125*h_x(i-1,j))-(0.375*h_x(i+1,j)+0.75*h_x(i,j)-0.125*h_x(i-1,j)) + (0.375*h_y(i,j)+0.75*h_y(i,j-1)-0.125*h_y(i,j-1)) - (0.375*h_y(i,j+1)+0.75*h_y(i,j)-0.125*h_y(i,j-1)))
       elseif (i==2 & j==jmax)
           T(i,j) = T_old(i,j) + (Dt/(rho*cp*dx*dy))*((0.375*h_x(i,j)+0.75*h_x(i-1,j)-0.125*h_x(i-1,j))-(0.375*h_x(i+1,j)+0.75*h_x(i,j)-0.125*h_x(i-1,j)) + (0.375*h_y(i,j)+0.75*h_y(i,j-1)-0.125*h_y(i,j-2)) - (0.375*h_y(i,j)+0.75*h_y(i,j)-0.125*h_y(i,j-1)))
       elseif (i==imax & j~=2 & j~=jmax)
           T(i,j) = T_old(i,j) + (Dt/(rho*cp*dx*dy))*((0.375*h_x(i,j)+0.75*h_x(i-1,j)-0.125*h_x(i-2,j))-(0.375*h_x(i,j)+0.75*h_x(i,j)-0.125*h_x(i-1,j)) + (0.375*h_y(i,j)+0.75*h_y(i,j-1)-0.125*h_y(i,j-2)) - (0.375*h_y(i,j+1)+0.75*h_y(i,j)-0.125*h_y(i,j-1)))
       elseif (j==2 & i~= imax & i~= 2) 
           T(i,j) = T_old(i,j) + (Dt/(rho*cp*dx*dy))*((0.375*h_x(i,j)+0.75*h_x(i-1,j)-0.125*h_x(i-2,j))-(0.375*h_x(i+1,j)+0.75*h_x(i,j)-0.125*h_x(i-1,j)) + (0.375*h_y(i,j)+0.75*h_y(i,j-1)-0.125*h_y(i,j-1)) - (0.375*h_y(i,j+1)+0.75*h_y(i,j)-0.125*h_y(i,j-1)))
       elseif (j==2 & i==imax)
           T(i,j) = T_old(i,j) + (Dt/(rho*cp*dx*dy))*((0.375*h_x(i,j)+0.75*h_x(i-1,j)-0.125*h_x(i-2,j))-(0.375*h_x(i,j)+0.75*h_x(i,j)-0.125*h_x(i-1,j)) + (0.375*h_y(i,j)+0.75*h_y(i,j-1)-0.125*h_y(i,j-1)) - (0.375*h_y(i,j+1)+0.75*h_y(i,j)-0.125*h_y(i,j-1)))
       elseif (j==jmax & i~=2 & i~=imax)
           T(i,j) = T_old(i,j) + (Dt/(rho*cp*dx*dy))*((0.375*h_x(i,j)+0.75*h_x(i-1,j)-0.125*h_x(i-2,j))-(0.375*h_x(i+1,j)+0.75*h_x(i,j)-0.125*h_x(i-1,j)) + (0.375*h_y(i,j)+0.75*h_y(i,j-1)-0.125*h_y(i,j-2)) - (0.375*h_y(i,j)+0.75*h_y(i,j)-0.125*h_y(i,j-1)))
       elseif (i==imax & j==jmax)
           T(i,j) = T_old(i,j) + (Dt/(rho*cp*dx*dy))*((0.375*h_x(i,j)+0.75*h_x(i-1,j)-0.125*h_x(i-2,j))-(0.375*h_x(i,j)+0.75*h_x(i,j)-0.125*h_x(i-1,j)) + (0.375*h_y(i,j)+0.75*h_y(i,j-1)-0.125*h_y(i,j-2)) - (0.375*h_y(i,j)+0.75*h_y(i,j)-0.125*h_y(i,j-1)))
       else
           T(i,j) = T_old(i,j) + (Dt/(rho*cp*dx*dy))*((0.375*h_x(i,j)+0.75*h_x(i-1,j)-0.125*h_x(i-2,j))-(0.375*h_x(i+1,j)+0.75*h_x(i,j)-0.125*h_x(i-1,j)) + (0.375*h_y(i,j)+0.75*h_y(i,j-1)-0.125*h_y(i,j-2)) - (0.375*h_y(i,j+1)+0.75*h_y(i,j)-0.125*h_y(i,j-1)))
           end
       end
   end
   unsteadiness=max(abs(T-T_old))/DTc;
   unsteadiness_nd = unsteadiness*(L/(u*Dt));    
   printf("Time step no. %5d, Unsteadiness_nd = %8.4e\n", N , unsteadiness_nd); 
end
xset("colormap", jetcolormap(64)),colorbar(0,100),Sgrayplot(linspace(0,L,imax),linspace(0,L,jmax),T, strf='100')
// SOU
//for i = 2:imax
//       for j = 2:imax
//           if (i == 2 & j ~=2)    
//             T(i,j) = T_old(i,j) + (Dt/(rho*cp*dx*dy))*((1.5*h_x(i-1,j)-0.5*h_x(i-1,j))-(1.5*h_x(i,j)-0.5*h_x(i-1,j)) + (1.5*h_y(i,j-1)-0.5*h_y(i,j-2)) - (1.5*h_y(i,j)-0.5*h_y(i,j-1)))  
//           //T(i,j) = T_old(i,j) + (Dt/(rho*cp*dx*dy))*(h_x(i-1,j)-h_x(i,j) + h_y(i,j-1) - h_y(i,j))//FOU
//       elseif (j == 2 & i ~=2)
//           T(i,j) = T_old(i,j) + (Dt/(rho*cp*dx*dy))*((1.5*h_x(i-1,j)-0.5*h_x(i-2,j))-(1.5*h_x(i,j)-0.5*h_x(i-1,j)) + (1.5*h_y(i,j-1)-0.5*h_y(i,j-1)) - (1.5*h_y(i,j)-0.5*h_y(i,j-1)))
//       elseif (i==2 & j==2)
//           T(i,j) = T_old(i,j) + (Dt/(rho*cp*dx*dy))*((1.5*h_x(i-1,j)-0.5*h_x(i-1,j))-(1.5*h_x(i,j)-0.5*h_x(i-1,j)) + (1.5*h_y(i,j-1)-0.5*h_y(i,j-1)) - (1.5*h_y(i,j)-0.5*h_y(i,j-1)))
//       else 
//           T(i,j) = T_old(i,j) + (Dt/(rho*cp*dx*dy))*((1.5*h_x(i-1,j)-0.5*h_x(i-2,j))-(1.5*h_x(i,j)-0.5*h_x(i-1,j)) + (1.5*h_y(i,j-1)-0.5*h_y(i,j-2)) - (1.5*h_y(i,j)-0.5*h_y(i,j-1)))
//           end
//       end
//   end
//   unsteadiness=max(abs(T-T_old))/DTc;
//   unsteadiness_nd = unsteadiness*(L/(u*Dt));    
//   printf("Time step no. %5d, Unsteadiness_nd = %8.4e\n", N , unsteadiness_nd); 
//end

