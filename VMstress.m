%Code to determine the yield stress of the Ag junction from in situ
%experiments
clc
close all 
load Raw_data.mat
DFf=[8.9 8 27 40 16 16]; %nN: uncertainty in friction force
DFn=[7.6 7.4 37 55 21 20]; %nN: uncertainty in normal force
Dd=.9; %nm: uncertainty in position
figure 
hold on
for i=1:6
clearvars -except i DFf DFn Dd rawdata stressdata
dist=rawdata(i).dist; dist=dist-min(dist);Ff=rawdata(i).Ff; diam=rawdata(i).diam; Fn=rawdata(i).Fn;  
[Y,I]=sort(dist); dist=Y;Ff=Ff(I);Fn=Fn(I);diam=diam(I);
shear= Ff./(pi*(diam.^2));
Dshear=sqrt((DFf(i)./((pi*diam).^2/4)).^2+(Dd*2*Ff./(pi*diam.^3/4)).^2);
norm=(Fn./(pi*(diam/2).^2));
Dnorm=sqrt((DFn(i)./(pi*diam.^2/4)).^2+(Dd*2*Fn./(pi*diam.^3/4)).^2);
width=12;height=12; 
VM=sqrt(norm.^2+3*shear.^2);
DVM=sqrt((.5*VM.^-1.*(2*norm.*Dnorm)).^2+(.5*VM.^-1.*(6*shear.*Dshear)).^2);

stressdata(i).dist=dist;
stressdata(i).VM=VM;
stressdata(i).norm=norm;
stressdata(i).shear=shear;

scatter(stressdata(i).dist,stressdata(i).VM)
end
clearvars -except stressdata rawdata
legend('1 experiment','2','3','4','5','6','Location','north')
xlabel('distance [nm]')
ylabel('Von Mises stress [GPa]')