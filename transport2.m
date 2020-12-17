clc
phi=0.3;
uw=0.001;
%define t and x
ts=3000;
xs=100;
Cinj=1;
D=0.01;
tn=200;
tot=xs+1;
dx=0.1;
dt=tn/ts;
for i=1:xs
    x(i)=i*dx;
end
for i=1:tn
    t(i)=i*dt;
end

C(1:(tot),1:ts)=0;

Sw=0.4;

for m=1:1
const1=D*dt/(2*dx*dx);
const2=uw*dt/(2*phi*Sw*dx);
vectorofone=ones(xs,1);
A(1:xs,1:xs)=0; 
B(1:xs,1)=0;
A=diag((2*const1+1+const2)*vectorofone,0)+diag(-const1*vectorofone(1:xs-1),1)+diag((-const1-const2)*vectorofone(1:xs-1),-1);



for j=2:ts
for i=1:xs
     
if i==1
    B(i,1)=C(i,j-1)+const1*(C(i+1,j-1)+Cinj-2*C(i,j-1))-const2*(C(i,j-1)-Cinj)+Cinj*(const1+const2);
elseif i==xs
     B(i,1)=C(i,j-1)+const1*(C(i-1,j-1)-2*C(i,j-1))-const2*(C(i,j-1)-C(i-1,j-1));
else
  B(i,1)=C(i,j-1)+const1*(C(i+1,j-1)+C(i-1,j-1)-2*C(i,j-1))-const2*(C(i,j-1)-C(i-1,j-1));
end  
end

New=ThomasAlgo(A,B);
C(1:xs,j)=New; 

end

end
V=C(1:xs,ts);
W=[Cinj; V];

plot(x,W(1:xs,1));

hold on;

for i=1:xs
C2(i,1)=Cinj/2*(erfc((x(i)-uw*tn)/(2*sqrt(D*tn)))+exp(x(i)*uw/D)*(erfc((x(i)+uw*tn)/(2*sqrt(D*tn)))));
end
 
C3=[Cinj;C2];
plot(x,C3(1:xs,1));
title('COMPARISON BETWEEN NUMERICAL AND ANALYTICAL VALUES','Color','m');
xlabel('Grid number');
ylabel('Concentration');
LEG=legend('Numerical solution','Analytical Solution');
LEG.FontSize = 14;
r = corrcoef(C3(1:xs), W(1:xs));
disp(r(1,2));

str=['       R= ',num2str(r(1,2))]
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'left');
hold off;
