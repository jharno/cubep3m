%SetDefault
close all
!scp bnu_ztj_haoran@172.16.21.103:/vol-th/home/bnu_ztj_haoran/haoran/tides/LOS5/log_*.txt .
a=load('log_dm.txt');
b=load('log_nu.txt');
z_max=20;

% nts, z, time/hour, den_buf

va = diff([0;a(:,3)]) * 3600; % speed^-1: sec/timestep
vb = diff([0;b(:,3)]) * 3600; % speed^-1: sec/timestep
v = [va;vb];

nhour_max = max(a(:,3)) + max(b(:,3));
nts_max=max(b(:,1));

nt=[a(:,1);b(:,1)];
z=[a(:,2);b(:,2)];
t=[a(:,3);max(a(:,3))+b(:,3)];
db=[a(:,4)*9;b(:,4)];

subplot(2,1,1)

figure(1);hold
    [AX,H1,H2]=plotyy(nt,z,nt,db,'plot','plot');
    set(get(AX(1),'Ylabel'),'String','redshift') 
    set(get(AX(2),'Ylabel'),'String','density buffer')
    %set(H1,'LineStyle','--','Color',[0,0,0])
    %set(H2,'LineStyle','-','Color',[0,0,0])

time_int = 1; % /hour    
    
for nhour = time_int : time_int : nhour_max
    nts_temp = interp1(t,nt,nhour);
    figure(1)
    plot([nts_temp,nts_temp],[0,z_max],'r--')
    text(nts_temp+0.02*nts_max,0.97*z_max,[num2str(nhour),'h'],'Color','r')
end    
    
subplot(2,1,2)
    plot(nt,v)
    xlabel('timestep')
    ylabel('sec/timestep')
