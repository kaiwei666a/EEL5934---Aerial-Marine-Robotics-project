figure(1)
xd = out.xd;
yd = out.yd;
plot(xd,yd,'red','Linewidth',2)
hold on
x = out.x;
y = out.y;
plot(x,y,'blue')
legend('desired', 'actual');
xlabel('x(m)');
ylabel('y(m)');

figure(2)
t = out.t;
x_e = out.x_e;
y_e = out.y_e;
plot(t,x_e,'blue')
hold on
plot(t,y_e,'red')
legend('x_e(m)', 'y_e(m)');
xlabel('t(s)')

figure(3)
u = out.u;
v = out.v;
r = out.r;
plot(t,u,'blue')
hold on
plot(t,v,'red')
hold on
plot(t,r,'magenta')
legend('u(m/s)', 'v(m/s)','r(m/s)');
xlabel('t(s)')

figure(4)
e_u = out.e_u;
e_v = out.e_v;
plot(t,e_u,'blue')
hold on
plot(t,e_v,'red')
legend('e_u(m/s)', 'e_v(m/s)');
xlabel('t(s)')

figure(5)
subplot(2,1,1)
t_new = linspace(0,20,2000);
tau_u = out.tau_u;
plot(t_new,tau_u,'blue')
xlabel('t(s)')
ylabel('tau_u(N)')
subplot(2,1,2)
tau_r = out.tau_r;
plot(t_new,tau_r,'blue')
xlabel('t(s)')
ylabel('tau_r(N*M)')
