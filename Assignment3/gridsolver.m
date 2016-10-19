clear
clc
% User-Input
rho = 1; L=1; Re=400; mu=1/Re;
imax=42; jmax=42;
epsilon_st= 10e-8; 

% Geometrical Parameters
Dx = L/(imax-2); Dy = L/(jmax-2); Dt=0.001;

% Coefficients in pressure correction
aE=(Dt*rho*Dy)/Dx;  aW=(Dt*rho*Dy)/Dx;
aN=(Dt*rho*Dx)/Dy;  aS=(Dt*rho*Dx)/Dy;
aP=aE+aW+aN+aS;

% IC and BCs
u=zeros(jmax,imax-1);  v=zeros(jmax-1,imax);  p=zeros(jmax,imax);  P_correction=zeros(jmax,imax); u_star=zeros(jmax,imax-1);  v_star=zeros(jmax-1,imax);
u_star(:,1)=0; v_star(:,1)=0; P_correction(:,1)= P_correction(:,2);
u_star(:,imax-1)=0; v_star(:,imax)=0; P_correction(:,imax)= P_correction(:,imax-1);
u_star(1,:)=0; v_star(1,:)=0; P_correction(1,:)= P_correction(2,:);
u_star(jmax,:)=1; v_star(jmax-1,:)=0; P_correction(jmax,:)= P_correction(jmax-1,:);


x=linspace(0,L,imax-1);    % Coordinates of face centers
y=linspace(0,L,jmax-1); 

% Coordinates of Cell Centers; NEEDED FOR PLOTTING.
xc(1)=x(1); xc(imax)=x(imax-1); 
for i=2:imax-1
xc(i)=(x(i)+x(i-1))/2;
end
yc(1)=y(1); yc(jmax)=y(jmax-1); 
for i=2:jmax-1
yc(i)=(y(i)+y(i-1))/2;
end 

% Grid points ploting
figure;
for i=1:imax-1
    plot(x(i),yc(2:jmax-1),'g-s'); hold on;
end
for j=1:jmax-1
    plot(xc(2:imax-1),y(j),'r-v'); hold on;
end
for i=1:imax
    plot(xc(i),yc,'y-o'); hold on;
end
hold off;
 
unsteadiness_nd=1; 
% Time-Marching to Explicit Unsteady State LAEs
while unsteadiness_nd>=0.0001
    u_old=u_star;
    v_old=v_star;
    
    for j=2:jmax-1                                % Fluxes in X direction to find u velocity
       for i=1:imax-2 
         mx1(j,i)=((u_old(j,i)+u_old(j,i+1))*rho)/2;
         ax1(j,i)=(max(mx1(j,i),0)*u_old(j,i))+(min(mx1(j,i),0)*u_old(j,i+1));
         dx1(j,i)=((u_old(j,i+1)-u_old(j,i))/Dx)*mu;
       end 
    end
    
    for j=1:jmax-1                               % Fluxes in Y direction to find u velocity
       for i=2:imax-2 
         my1(j,i)=((v_old(j,i)+v_old(j,i+1))*rho)/2;
         ay1(j,i)=(max(my1(j,i),0)*u_old(j,i))+(min(my1(j,i),0)*u_old(j+1,i));
         dy1(j,i)=((u_old(j+1,i)-u_old(j,i))/Dy)*mu;
       end 
    end
    
    for j=2:jmax-1                                % Temporal evolution of u velocity
       for i=2:imax-2 
         A(j,i)=((ax1(j,i)-ax1(j,i-1))*Dy)+((ay1(j,i)-ay1(j-1,i))*Dx);
         D(j,i)=((dx1(j,i)-dx1(j,i-1))*Dy)+((dy1(j,i)-dy1(j-1,i))*Dx);
         S(j,i)=(p(j,i)-p(j,i+1))*Dy;
         u_star(j,i)=u_old(j,i)+((Dt/(rho*Dx*Dy))*(D(j,i)-A(j,i)+S(j,i))) ;
       end 
    end
    
    for j=2:jmax-2                                % Fluxes in X direction to find v velocity
       for i=1:imax-1 
         mx2(j,i)=((u_old(j,i)+u_old(j+1,i))*rho)/2;
         ax2(j,i)=(max(mx2(j,i),0)*v_old(j,i))+(min(mx2(j,i),0)*v_old(j,i+1));
         dx2(j,i)=((v_old(j,i+1)-v_old(j,i))/Dx)*mu;
       end 
    end
    
    for j=1:jmax-2                               % Fluxes in Y direction to find v velocity
       for i=2:imax-1 
         my2(j,i)=((v_old(j,i)+v_old(j+1,i))*rho)/2;
         ay2(j,i)=(max(my2(j,i),0)*v_old(j,i))+(min(my2(j,i),0)*v_old(j+1,i));
         dy2(j,i)=((v_old(j+1,i)-v_old(j,i))/Dy)*mu;
       end 
    end
    
    for j=2:jmax-2                               % Temporal evolution of  v velocity
       for i=2:imax-1 
         A1(j,i)=((ax2(j,i)-ax2(j,i-1))*Dy)+((ay2(j,i)-ay2(j-1,i))*Dx);
         D1(j,i)=((dx2(j,i)-dx2(j,i-1))*Dy)+((dy2(j,i)-dy2(j-1,i))*Dx);
         S1(j,i)=(p(j,i)-p(j+1,i))*Dx;
         v_star(j,i)=v_old(j,i)+((Dt/(rho*Dx*Dy))*(D1(j,i)-A1(j,i)+S1(j,i))) ;
       end 
    end
         
 while (1)      % Pressure correction loop
  for j=2:jmax-1
    for i=2:imax-1
        Div(j,i)=((u_star(j,i)-u_star(j,i-1))*Dy)+((v_star(j,i)-v_star(j-1,i))*Dx);
    end    
  end  
  if (max(max(abs(Div(j,i)))) <= epsilon_st)     
      break;

  else
      Error=1;
    while (Error>=epsilon_st)
    P_old=P_correction; 
     for j=2:jmax-1
       for i=2:imax-1
        P_correction(j,i)=((aE*P_correction(j,i+1))+(aW*P_correction(j,i-1))+(aN*P_correction(j+1,i))+(aS*P_correction(j-1,i))-(Div(j,i)))/aP;
       end    
     end 
        P_correction(:,imax)= P_correction(:,imax-1); 
        P_correction(jmax,:)= P_correction(jmax-1,:);
        P_correction(:,1)= P_correction(:,2); 
        P_correction(1,:)= P_correction(2,:); 
        Error=max(max(abs(P_correction-P_old)));
    end
    for j=2:jmax-1                              
       for i=2:imax-2 
          u_star(j,i)=u_star(j,i)+ ((Dt/(rho*Dx))*(P_correction(j,i)-P_correction(j,i+1))) ;
       end 
     end
    
    for j=2:jmax-2                              
       for i=2:imax-1 
          v_star(j,i)=v_star(j,i) + ((Dt/(rho*Dy))*(P_correction(j,i)-P_correction(j+1,i))) ;
       end 
    end  
  end   
      
 end 
for j=1:jmax
    for i=1:imax
        p(j,i) = p(j,i) + P_correction(j,i);
    end
end 
u_star(:,imax-1)=0;  
v_star(jmax-1,:)=0; 
unsteadiness_nd=max((max(max(abs(u_star-u_old)))),(max(max(abs(v_star-v_old)))));  % Steady state convergence criterion.
unsteadiness_nd 
end

v=linspace(-1,1,51); v1=linspace(-1,1,26);
figure; 
[C,h]=contourf(x,yc,u_star,v1); clabel(C,h); xlabel('X'); ylabel('Y'); title('Steady State U Velocity Contour');
figure;
[C,h]=contourf(xc,yc,p,v); clabel(C,h); xlabel('X'); ylabel('Y'); title('Steady State V Velocity Contour');

ypA = [1 0.9766 0.9688 0.9609 0.9531 0.8516 0.7344 0.6172 0.5 0.4531 0.2813 0.1719 0.1016 0.0703 0.0625 0.0547 0];           % Benchmark Results
ucA100 = [1 0.84123 0.78871 0.73722 0.68717 0.23151 0.00332 -0.13641 -0.20581 -0.2109 -0.15662 -0.1015 -0.06434 -0.04775 -0.04192 -0.03717 0];
ucA400 = [1 0.75837 0.68439 0.61756 0.55892 0.29093 0.16256 0.02135 -0.11477 -0.17119 -0.32726 -0.24299 -0.14612 -0.10338 -0.09266 -0.08186 0];
vcA100 = [0 -0.05906 -0.07391 -0.08864 -0.10313 -0.16914 -0.22445 -0.24533 0.05454 0.17527 0.17507 0.16077 0.12317 0.1089 0.100091 0.09233 0];
vcA400 = [0 -0.12146 -0.15663 -0.19254 -0.22847 -0.23827 -0.44993 -0.38598 0.05186 0.30174 0.30203 0.28124 0.22965 0.2092 0.19713 0.1836 0];
vcA1000 = [0 -0.21388 -0.27669 -0.33714 -0.39188 -0.5155 -0.42665 -0.31966 0.02526 0.32235 0.33075 0.37095 0.32627 0.30353 0.29012 0.27485 0];
xpA = [1 0.9688 0.9609 0.9531 0.9453 0.9063 0.8594 0.8047 0.5 0.2344 0.2266 0.1563 0.0938 0.0781 0.0703 0.0625 0];

figure;
plot((u_star(:,imax/2)+u_star(:,(imax/2)+1))/2,yc,'r-s'); xlabel('U Velocity','FontSize',13,'FontWeight','bold'); ylabel('Y','FontSize',13,'FontWeight','bold'); title('U'); title('Variation of U-Velocity along the Verticle Centerline'); grid on; hold on;
if (Re==100)
    plot(ucA100,ypA,'k-o'); hold off;
end

if (Re==400) 
    plot(ucA400,ypA,'k-o'); hold off;
end
legend('By Code','Benchmark Result','Location','southeast') ;
figure;
plot(xc,(v_star(jmax/2,:)+v_star((jmax/2)+1,:))/2,'g-*'); xlabel('X','FontSize',13,'FontWeight','bold'); ylabel('V Velolcity','FontSize',13,'FontWeight','bold'); title('U'); title('Variation of V-Velocity along the Horizontal Centerline'); grid on; hold on;
if (Re==100) 
    plot(xpA,vcA100,'k-o'); hold off;
end

if (Re==400)
    plot(xpA,vcA400,'k-o'); hold off;
end
legend('By Code','Benchmark Result') ;

figure;
u_star(:,imax)=0; v_star(jmax,:)=0;
quiver(xc,yc,u_star,v_star,1.25); xlabel('X','FontSize',13,'FontWeight','bold'); ylabel('Y','FontSize',13,'FontWeight','bold'); title('U'); title('Velocity Vector'); 
