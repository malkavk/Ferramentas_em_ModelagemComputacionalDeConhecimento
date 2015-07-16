dt = 0.1;
n = 100;

measnoise = 10;
accelnoise = 0.2;
a = [1 dt; 0 1];
b = [dt^2/2; dt];
c=[1 1];
x = [0;0];
xhat = x;
Sz = measnoise^2;
Sw = accelnoise^2*[dt^4 dt^3/2; dt^3/2 dt^2];
P = [50 50; 50 50];

Sk = [];
qk = [];
Innt =[];
pos = [];
poshat =[];
posmeas = [];
vel = [];
velhat = [];
err = [];

for t=0:dt:n,
    u =1
    
    rand("normal")
    ProcessNoise = accelnoise * [(dt^2/2)*rand(); dt*rand()];
    x = a*x+b*u+ProcessNoise;
    MeasNoise = measnoise * rand();
    y = c*x+ MeasNoise;
    xhat = a*xhat + b*u;
    Inn = y -c*xhat;
    s = c*P*c' +Sz;
    K = a*P*c'*inv(s);
    xhat = xhat + K*Inn;
    P = a*P*a' -a*P*c'*inv(s)*c*P*a'+ Sw;
        Sk = [Sk, c*P*c'+ Sz];
    Innt = [Innt, Inn];
    qk = [qk, Inn*inv(c*P*c'+ Sz)*Inn]
    pos =[pos; x(1)];
    posmeas = [posmeas; y(1)];
    poshat = [poshat; xhat(1)];
    vel = [vel; x(2)];
    velhat = [velhat; xhat(2)];
    err =[err; abs(y(1)-xhat(1))]

end

//clf()
//t = [0:dt:n];
//plot(t,pos,t, posmeas, t, poshat);
//histplot(10, err)



//Teste de hipótese
aa =[];
for j = 1:1:n,
    aa = [aa,((poshat(j)-posmeas(j)).^2)/posmeas(j)];
end
xx = sum(aa);

//1º Gráfico
tt =[1:1:100];
plot(t,aa);

//2º Gráfico
tt =[1:1:100];
plot(tt, aa, t, -2*sqrt(Sk),t, 2*sqrt(Sk));

//3º Gráfico
plot(t, Innt^2);
