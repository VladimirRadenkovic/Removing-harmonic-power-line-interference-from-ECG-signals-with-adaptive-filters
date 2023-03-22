close all; clear; clc;
folder_path = 'Training_PTB'; 
file_pattern = '*.mat';  
file_list = dir(fullfile(folder_path, file_pattern));  
n = 20000;
N = 2^nextpow2(n);
Fs = 1000;
f = linspace(0,1,N/2+1)*Fs/2;
f_hp = 1; 
w_hp = f_hp/(Fs/2); 
[b_hp1, a_hp1] = butter(1, w_hp, 'high'); 
[b_hp2, a_hp2] = butter(5, w_hp, 'high'); 
f_lp = 40 ; 
w_lp = f_lp/(Fs/2); 
[b_lp, a_lp] = butter(2, w_lp, 'low'); 
for i = 1:100
    file_name = fullfile(folder_path, file_list(i).name);
    z = load(file_name).val;% 
    x = z(1,:);  
    y1 = filter(b_lp, a_lp, filter(b_hp1, a_hp1, x));
    %figure, plot(y1);
end
for i = 1:length(file_list)
    file_name = fullfile(folder_path, file_list(i).name);
    z = load(file_name).val;
    for j = 1:12
        x = z(j,11001:11000+n); 
%         figure, plot(x)
        y1 = filter(b_lp, a_lp, filter(b_hp, a_hp1, x));
        figure, plot(y1);
    end
    close all
    
end



    

for i = 1:length(file_list)
    file_name = fullfile(folder_path, file_list(i).name);  
    z = load(file_name).val(10001:10000+n);  
    figure, plot(z)
    Z = fft(z,N)/length(z)*2;
    figure, plot(f,abs(Z(1:N/2+1)));
    f_lp = 40;
    w_lp = f_lp/(Fs/2); 
    [b_lp, a_lp] = butter(6, w_lp, 'low'); 
    y1 = filter(b_lp, a_lp, z);
    figure, plot(y1);
    Y = fft(y1,N)/length(y1)*2;
    figure,plot(f,abs(Y(1:N/2+1)));
    end
