clear all
close all
clc

jay = sqrt(-1);

data_cur = csvread("39bus_step001_v33_noitr_loadchange_new_3.csv", 1,0);
data_cur = csvread("39bus_step005_v33_itr_loadchange_new_compensateY_eachstep.csv", 1,0);
data_cur = csvread("39bus_step005_v33_noitr_loadchange_new_compensateY_eachstep.csv", 1,0);
vmag8_cur=data_cur(:,16);
vmag508_cur=data_cur(:,19);
vang8_cur=data_cur(:,22);
vang508_cur=data_cur(:,25);
vt8_cur = vmag8_cur.*(cos(vang8_cur)+jay*sin(vang8_cur));
vt508_cur = vmag508_cur.*(cos(vang508_cur)+jay*sin(vang508_cur));
cur508_cur = (vt8_cur-vt508_cur)/(jay*0.003);
s508_cur = vt508_cur.*conj(cur508_cur);

% data_pq = csvread("39bus_step001_v33_itr_loadchange_new_2.csv", 1,0);
% vmag8_pq=data_pq(:,16);
% vmag508_pq=data_pq(:,19);
% vang8_pq=data_pq(:,22);
% vang508_pq=data_pq(:,25);
% vt8_pq = vmag8_pq.*(cos(vang8_pq)+jay*sin(vang8_pq));
% vt508_pq = vmag508_pq.*(cos(vang508_pq)+jay*sin(vang508_pq));
% cur508_pq = (vt8_pq-vt508_pq)/(jay*0.003);
% s508_pq = vt508_pq.*conj(cur508_pq);


figure();
plot(real(s508_cur), 'r'); hold on;
% plot(real(s508_pq), 'b'); hold off;
% title ('real power at bus 508 compare');
% legend('Constant Current', 'Constant Power');

figure();
plot(imag(s508_cur), 'r'); hold on;

figure();
plot(abs(vt508_cur), 'r');hold on;
% plot(abs(vt508_pq), 'b');hold off;
% title ('voltage magnitude at bus 508 compare');
% legend('Constant Current', 'Constant Power');


% % figure();
% % plot(real(cur508_cur))
% % 
% % figure();
% % plot(imag(cur508_cur))
% % 
% % figure();
% % plot(real(s508_cur)./abs(vt508_cur))
% % 
% % figure();
% % plot(imag(s508_cur)./abs(vt508_cur))
