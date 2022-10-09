close all 
clear all


%%
tempIn = 20; 
SalIn = 0.003; 
pressIn = 1; 
PO4In = 0;
SIn = 0;
TAlkIn = 1.5: 1e-2 : 3;
DICin = (1.5 : 1e-2 : 3)';  % important:  note orientation of vector.  
Ca = 0;
Mg = 0; 
SO4 = 0; 

out = co2syslite_v0_2_test(tempIn,SalIn,pressIn,PO4In,SIn,TAlkIn,DICin,Ca,Mg,SO4);
%%

% - So your co2syslite_v0_2 now produces an X x Y x 4 array, where X is the
% l3ngth of the DIC input and Y is the length of the ALK input.  so you can
% run this function on vectors of DICs and ALKs, and e.g. produce controur
% plots in a single hit:  

figure(99)

subplot(2,2,1)
out1 = out(:,:,1);
contourf(TAlkIn, DICin, out1);
xlabel('TA')
ylabel('DIC')
c = colorbar;
c.Label.String = 'pH';
title('Subplot 1: pH')

subplot(2,2,2)
out2 = out(:,:,2);
contourf(TAlkIn, DICin, out2);
xlabel('TA')
ylabel('DIC')
c = colorbar;
c.Label.String = 'CO2';
title('Subplot 2: CO2')

subplot(2,2,3)
out3 = out(:,:,3);
contourf(TAlkIn, DICin, out3);
xlabel('TA')
ylabel('DIC')
c = colorbar;
c.Label.String = 'Omega';
title('Subplot 3: Omega')

subplot(2,2,4)
out4 = out(:,:,4);
contourf(TAlkIn, DICin, out4);
xlabel('TA')
ylabel('DIC')
c = colorbar;
c.Label.String = 'CO3';
title('Subplot 4: CO3')





