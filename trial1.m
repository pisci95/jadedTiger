sigma = 1;
offset = 3;
up = 5;
total = 10;
symbol = [zeros(1,offset), ones(1,up), zeros(1,total-(up+offset))];
token = conv([1], symbol);
token = conv([1], token);

msg = [ 1 0 1 0  ];
pilot = [0, 0, 1];
msg = [pilot, msg];
symbolSize = 10;
fin = [];
xk = [];
goTo = 5;

for n = 1:length(msg)
 if (msg(n) == 1)
     fin = [fin, ones(1, symbolSize)];
 else 
     fin = [fin, zeros(1, symbolSize)];
 end
end
for i = 1:length(fin)
 if (fin(i) == 1)
     xk = [xk, 1];
 else 
     xk = [xk, -1];
 end
end
c = conv(symbol, xk);
transmitsignalZ = c;
     
%plot(c)
a = randn(size(transmitsignalZ));
b = randn(size(transmitsignalZ));
noiseZ= sigma/sqrt(2)*(a + j*b);

rsZ = transmitsignalZ + noiseZ;

rsZ = conv(rsZ, symbol);

%syncing
s = -1;
for i = 1:(length(rsZ)-31)
    if(s == -1)
        bool = up*length(pilot);
        for j = 1:length(pilot)
            ex = pilot(j);
            for k = 1:up
                t = rsZ(i + k +(j-1)*total+offset);
                if ((ex == 1) & (t > .1))
                    bool = bool -1;
                elseif ((ex == 0) & (t < -.1))
                    bool = bool -1;
                    
                end
            end
        end
        if(bool == 0)
            s = i + up*length(pilot);
            break;
        end
        i
    end
end
s
%rsZ = rsZ(s:end);


figure(1)
clf
subplot(2,1,1)
plot(real(transmitsignalZ),'b')
hold on
plot(imag(transmitsignalZ),'r')
legend('real','imag')
ylabel('xI(t)  and  xQ(t)')
xlabel('Time in samples')
subplot(2,1,2)
plot(real(rsZ),'b')
hold on
plot(imag(rsZ),'r')
zoom xon
legend('real','imag')
ylabel('yI(t)  and  yQ(t)')
xlabel('Time in samples')
hold off