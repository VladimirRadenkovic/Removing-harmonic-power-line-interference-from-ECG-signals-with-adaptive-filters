clear all;
close all;
clc

Fs = 1000;
mi1 = 0.01;
mi2 = 0.01;
mi3 = 0.01;
mi4 = 0.01;

load('Z.mat');
load('D.mat');


rng(1);
t = 0:1/Fs:(6000 - 1)/Fs;
rng(1);
A = 3*rand(1);
x_f50_sum = cos(2*pi*50*t) +  A ;
w1 = zeros(1,length(t)); 
w2 = zeros(1,length(t)); 
w3 = zeros(1,length(t)); 
w4 = zeros(1,length(t)); 
x_f50 = zeros(1, length(t));
x_f100 = zeros(1, length(t));
x_f150 = zeros(1, length(t));
x_f200 = zeros(1, length(t));

for i = 6:length(t)
    x_f50(i) = x_f50_sum(i) - w1(i);
    w1(i + 1) = w1(i) + mi1*x_f50(i);
%     xac = x_f50_sum(i:-1:i-4);
%     x_f(i) = 2*20/(pi*(xac*xac' + 0.000001))/xac(1);
    x_f100(i) = 2*x_f50(i)^2 - w2(i);
    w2(i + 1) = w2(i) + mi2*x_f100(i);
    x_f150(i) = 4*x_f50(i)^3 - 3*x_f50(i) - w3(i);
    w3(i + 1) = w3(i) + mi3*x_f150(i);
    x_f200(i) = 8*x_f50(i)^4 - 4*x_f100(i) - w4(i);
    w4(i + 1) = w4(i) + mi4*x_f200(i);
end

figure
hold on;
plot(w1);
plot(w2);
plot(w3);
plot(w4);
legend('w1','w2','w3','w4'); title('Težine za predobradu reference tokom LMS algoritma');xlabel('t[ms]');
hold off;

%dodaj
for i = 6:length(t)
    x_f50(i) = x_f50_sum(i) - w1(end);
    x_f100(i) = 2*x_f50(i)^2 - w2(end);
    x_f150(i) = 4*x_f50(i)^3 - 3*x_f50(i) - w3(end);
    x_f200(i) = 8*x_f50(i)^4 - 4*x_f100(i) - w4(end);
end

N = 2^nextpow2(length(x_f50));
f = linspace(0,1,N/2+1)*Fs/2;
Xf50 = fft(x_f50,N)/length(x_f50)*2;
Xf100 = fft(x_f100,N)/length(x_f100)*2;
Xf150 = fft(x_f150,N)/length(x_f150)*2;
Xf200 = fft(x_f200,N)/length(x_f200)*2;

X = Xf50+Xf100+Xf150+Xf200;


figure; hold all;
plot(f,abs(X(1:N/2+1)).^2);
title('Spektar snage reference nakon filtriranja'); xlabel('f[Hz]'); ylabel('Snaga')

%%
mi = 0.01;
for j = 1:10
z = Z(j,:);
d = D(j,:);
e = zeros(1, 6000);
W = zeros(8, 6000);

for i = 6:6000
  
    x_q50 =  x_f50(i - round(Fs./(4*50)));
    x_q100 =  x_f100(i - round(Fs./(4*100)));
    x_q150 =  x_f150(i - round(Fs./(4*150)));
    x_q200 =  x_f200(i - round(Fs./(4*200)));
    xn  = [x_f50(i); x_q50; x_f100(i); x_q100; x_f150(i); x_q150; x_f200(i); x_q200];
    
    e(i) = d(i) - W(:,i)'*xn;
    W(:,i + 1) = W(:, i) + mi*e(i)*xn;
    
end

for i = 6:6000
  
    x_q50 =  x_f50(i - round(Fs./(4*50)));
    x_q100 =  x_f100(i - round(Fs./(4*100)));
    x_q150 =  x_f150(i - round(Fs./(4*150)));
    x_q200 =  x_f200(i - round(Fs./(4*200)));
    xn  = [x_f50(i); x_q50; x_f100(i); x_q100; x_f150(i); x_q150; x_f200(i); x_q200];
    
    e(i) = d(i) - W(:,end-1)'*xn;
    
end

figure
plot(z); title('Željeni signal z(n)');
xlabel('t[ms]');ylabel('Amplituda');
%filename = ['Prez\Željeni signal z(n)_', num2str(j), '.jpg'];saveas(gcf, filename, 'jpg');
figure
plot(d); title('Signal sa smentnjama d(n)');
xlabel('t[ms]');ylabel('Amplituda');
%filename = ['Prez\Signal sa smentnjama d(n)_', num2str(j), '.jpg'];saveas(gcf, filename, 'jpg');
figure
plot(e); title('Signal greške e(n)');
xlabel('t[ms]');ylabel('Amplituda');
%filename = ['Prez\Signal greške e(n)_', num2str(j), '.jpg'];saveas(gcf, filename, 'jpg');

N = 2^nextpow2(length(z));
f = linspace(0,1,N/2+1)*Fs/2;

D1 = fft(d,N)/length(d)*2;
Z1 = fft(z,N)/length(z)*2;
E = fft(e,N)/length(e)*2;

figure, plot(f,abs(Z1(1:N/2+1)).^2);
title('Spektar snage željenog signala z(n)'); xlabel('f[Hz]'); ylabel('Snaga');
%filename = ['Prez\Spektar snage željenog signala z(n)_', num2str(j), '.jpg'];saveas(gcf, filename, 'jpg');
figure, plot(f,abs(D1(1:N/2+1)).^2);
title('Spektar snage signala sa smetnjama d(n)'); xlabel('f[Hz]'); ylabel('Snaga');
%filename = ['Prez\Spektar snage signala sa smetnjama d(n)_', num2str(j), '.jpg'];saveas(gcf, filename, 'jpg');
figure, plot(f,abs(E(1:N/2+1)).^2);
title('Spektar snage signala greške e(n)'); xlabel('f[Hz]'); ylabel('Snaga');
%filename = ['Prez\Spektar snage signala greške e(n)_', num2str(j), '.jpg'];saveas(gcf, filename, 'jpg');

figure
hold all;
plot(W(1,:));
plot(W(2,:));
plot(W(3,:));
plot(W(4,:));
plot(W(5,:));
plot(W(6,:));
plot(W(7,:));
plot(W(8,:));
legend('w1','w2','w3','w4','w5','w6','w7','w8'); title('Težine LMS algoritma za uklanjanje smetnji iz EKG signala \mu = 0.01');xlabel('t[ms]');
%filename = ['Prez\Težine LMS algoritma za uklanjanje smetnji iz EKG signala_', num2str(j), '.jpg'];saveas(gcf, filename, 'jpg');
end




%%
% idealni slucaj
x_f50 = cos(2*pi*50*t);
x_f100 = 2*x_f50.^2 - 1;
x_f150 = 4*x_f50.^3 - 3*x_f50;
x_f200 = 8*x_f50.^4 - 4*x_f100 - 3;

e = zeros(1, 6000);
W = zeros(8, 6000);

for i = 6:6000
  
    x_q50 =  x_f50(i - round(Fs./(4*50)));
    x_q100 =  x_f100(i - round(Fs./(4*100)));
    x_q150 =  x_f150(i - round(Fs./(4*150)));
    x_q200 =  x_f200(i - round(Fs./(4*200)));
    xn  = [x_f50(i); x_q50; x_f100(i); x_q100; x_f150(i); x_q150; x_f200(i); x_q200];
    
    e(i) = d(i) - W(:,i)'*xn;
    W(:,i + 1) = W(:, i) + mi*e(i)*xn;
    
end

for i = 6:6000
  
    x_q50 =  x_f50(i - round(Fs./(4*50)));
    x_q100 =  x_f100(i - round(Fs./(4*100)));
    x_q150 =  x_f150(i - round(Fs./(4*150)));
    x_q200 =  x_f200(i - round(Fs./(4*200)));
    xn  = [x_f50(i); x_q50; x_f100(i); x_q100; x_f150(i); x_q150; x_f200(i); x_q200];
    
    e(i) = d(i) - W(:,end)'*xn;
    
end

figure
plot(z); title('Željeni signal (n)');
xlabel('t[ms]');ylabel('Amplituda');
figure
plot(d); title('Signal sa smentnjama d(n)');
xlabel('t[ms]');ylabel('Amplituda');
figure
plot(e); title('Signal greške e(n)');
xlabel('t[ms]');ylabel('Amplituda');

N = 2^nextpow2(length(z));
f = linspace(0,1,N/2+1)*Fs/2;

D = fft(d,N)/length(d)*2;
Z = fft(z,N)/length(z)*2;
E = fft(e,N)/length(e)*2;

figure, plot(f,20*log10(abs((Z(1:N/2+1))).^2));
title('Spektar snage željenog signala z(n)'); xlabel('f[Hz]'); ylabel('Snaga[dB]')
figure, plot(f,20*log10(abs((D(1:N/2+1))).^2));
title('Spektar snage signala sa smetnjama d(n)'); xlabel('f[Hz]'); ylabel('Snaga[dB]')
figure, plot(f,20*log10(abs((Z(1:N/2+1))).^2));
title('Spektar snage signala greške e(n)'); xlabel('f[Hz]'); ylabel('Snaga[dB]')

figure
hold all;
plot(W(1,:));
plot(W(2,:));
plot(W(3,:));
plot(W(4,:));
plot(W(5,:));
plot(W(6,:));
plot(W(7,:));
plot(W(8,:));
legend('w1','w2','w3','w4','w5','w6','w7','w8'); title('Težine LMS algoritma za uklanjanje smetnji iz EKG signala');xlabel('t[ms]');
hold off;
