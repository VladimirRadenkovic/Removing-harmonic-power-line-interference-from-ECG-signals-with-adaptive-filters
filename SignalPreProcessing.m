close all; clear; clc;
folder_path = 'Training_PTB';  
file_pattern = '*.mat';  
file_list = dir(fullfile(folder_path, file_pattern));
n = 20000;
N = 2^nextpow2(n);
Fs = 1000;
f = linspace(0,1,N/2+1)*Fs/2;
f_hp = 0.5; 
w_hp = f_hp/(Fs/2);
[b_hp, a_hp] = butter(1, w_hp, 'high'); 
f_lp1 = 40 ; 
w_lp1 = f_lp1/(Fs/2);
[b_lp1, a_lp] = butter(5, w_lp1, 'low'); 
    
D = zeros(length(file_list),6000);
Z = zeros(length(file_list),6000);
for i = 1:10
    file_name = fullfile(folder_path, file_list(i).name);  
    z = load(file_name).val(1,:);  
    N = 2^nextpow2(length(z));
    f = linspace(0,1,N/2+1)*Fs/2;
    Z1 = fft(z,N)/length(z)*2;
    y = filter(b_lp1, a_lp, filter(b_hp, a_hp, z));
    figure, plot(z)
    title('Primer željenog signala z(n) pre BP filtriranja'); xlabel('t[ms]'); ylabel('Amplituda');
    filename = ['Primer željenog signala z(n) pre BP filtriranja_', num2str(i), '.jpg'];saveas(gcf, filename, 'jpg');
    figure, plot(y)
    title('Primer željenog signala z(n) nakon BP filtriranja'); xlabel('t[ms]'); ylabel('Amplituda')
    filename = ['Primer željenog signala z(n) nakon BP filtriranja_', num2str(i), '.jpg'];saveas(gcf, filename, 'jpg');
    
    A = max(abs(y(1001:7000)));
    y = y(1001:7000)/A*100;
    t = 0:1/Fs:(length(y)-1)/Fs;
    x1 = 100*sin(2*pi*50*t + rand(1)*pi);
    x2 = 100*sin(2*pi*100*t + rand(1)*pi);
    x3 = 100*sin(2*pi*150*t + rand(1)*pi);
    x4 = 100*sin(2*pi*200*t + rand(1)*pi);
    y1 = y + x1 + x2 + x3 + x4;
    D(i,:) = y1;
    Z(i,:) = y;
    Y = fft(y,N)/length(y)*2;
    Y1 = fft(y1,N)/length(y1)*2;
    

   
   figure, plot(f,20*log10(abs((Y(1:N/2+1))).^2));
   title('Spektar snage željenog signala z(n) nakon BP filtriranja'); xlabel('f[Hz]'); ylabel('Snaga[dB]')
   filename = ['Spektar snage željenog signala z(n)_', num2str(i), '.jpg'];saveas(gcf, filename, 'jpg');
   
   figure, plot(f,20*log10(abs((Y1(1:N/2+1))).^2));   
   title('Spektar snage signala d(n) sa dodatim smetnjama'); xlabel('f[Hz]'); ylabel('Snaga[dB]')
   filename = ['Spektar snage signala d(n) sa dodatim smetnjama_', num2str(i), '.jpg'];saveas(gcf, filename, 'jpg');
   
   figure, plot(y)
   title('Željeni signal z(n)'); xlabel('t[ms]'); ylabel('Amplituda')
   filename = ['Željeni signal z(n)_', num2str(i), '.jpg'];saveas(gcf, filename, 'jpg');
   
   figure, plot(y1)
   title('Signal sa smetnjama d(n)'); xlabel('t[ms]'); ylabel('Amplituda')
   filename = ['Signal sa smetnjama d(n)_', num2str(i), '.jpg'];saveas(gcf, filename, 'jpg');
end

save('D.mat', 'D');
save('Z.mat', 'Z');
