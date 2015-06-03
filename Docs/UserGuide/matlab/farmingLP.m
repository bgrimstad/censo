
close all

A = [-4000 -3000;
       -60   -80];
b = [-100000 -2000]';
lb = [0 0]';
c = [-7000 -6000];
[x0,x1]=meshgrid(-1:0.1:26,-1:0.1:26);
f = c(1).*x0 + c(2).*x1;
plotregion(A,b,lb,[],[1 1 1],0.7)
set(gca,'Color',[0.9 0.9 0.9]);
hold on
contour(x0,x1,f,20);
xlabel('x_0 (Apples [tonnes])');
ylabel('x_1 (Bananas [tonnes])');
plot(14.2857,14.2857,'bo')
plot(0,0,'bo')
text(1,1,'z_0')
text(15,15,'x^*=(14.2857,14.2857)')
legend('Feasible area')
grid
axis ([-1 26 -1 26])