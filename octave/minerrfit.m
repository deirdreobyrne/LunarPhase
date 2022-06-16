pkg load optim;
format long;

ang=csvread("i.csv");
t = (ang([1:524288],1)-ang(1,1))*86400;
dt = 45.0 + 50.0*t/86400/36525;
t = t + dt / 86400;
phaseangle=ang([1:524288],2);
dt=t(2)-t(1);
n = rows(t);
freq=linspace(0,1/dt,n);
degjcy = freq*360*36525*86400;

phase = (1 - cos(phaseangle))/2;

args = [
   4.455915442782913e+00
   3.789022549880648e+05
   4.847299115959568e+00
   4.060747530046794e+05
   3.711660706912899e+00
   3.741949447565621e+05
   6.244955226928489e+00
   5.022680478319048e+06
   1.094058217228406e-01
  -3.645478117891136e-02
   2.217866961747020e-02
   1.439518894443362e-02
   3.727205257421440e-03
   1.952669837057546e-03
   1.106825135693690e-03
   1.008726072053177e-03
  -1.240627205414290e-03
   9.787354909447901e-04
   8.494457711525009e-04
   7.699865140427416e-04
  -7.125138902358125e-04
  -5.179824666526639e-04
  -9.960396371095429e-05
   2.747671633223807e-04
   8.579460694360643e-05
   1.732183870670095e-04
   1.484419451667800e-04
   8.376357047334948e-07
  -4.644009089094595e-06
   9.738762955525674e-07
   6.866396058689803e-05
  -4.767805350497088e-06
  -3.832300081587489e-06
  -2.770725440201805e-06
  -8.371312035697797e-05
  -2.442889793212899e-04
   9.412768821066482e-05
  -1.405875007764292e-04
   6.507285578148511e-05
   1.904091431941523e-03
   6.294098181004421e-05
  -4.980318753955695e-05
   9.120882679826181e-04
   ];

fitfunc = @(p) max(abs((1-cos( ...
    mod(polyval([1/p(4),p(3)],t),2*pi) ...
  + p(9)*sin(polyval([1/p(2),p(1)],t)) ... # mp
  + p(10)*sin(polyval([1/p(8),p(7)],t)) ... # M
  + p(11)*sin(2*polyval([1/p(4),p(3)],t)-polyval([1/p(2),p(1)],t)) ... # 2D-mp
  + p(12)*sin(2*polyval([1/p(4),p(3)],t)) ... # 2D
  + p(13)*sin(2*polyval([1/p(2),p(1)],t)) ... # 2mp
  + p(14)*sin(polyval([1/p(4),p(3)],t)) ... # D
  + p(15)*sin(2*(polyval([1/p(6),p(5)],t)-polyval([1/p(4),p(3)],t))) ... # 2(F-D)
  + p(16)*sin(2*(polyval([1/p(4),p(3)],t)-polyval([1/p(2),p(1)],t))) ... # 2(D-mp)
  + p(17)*sin(2*polyval([1/p(6),p(5)],t)) ... # 2F
  + p(18)*sin(2*polyval([1/p(4),p(3)],t)-polyval([1/p(8),p(7)],t)-polyval([1/p(2),p(1)],t)) ... # 2D - M - mp
  + p(19)*sin(2*polyval([1/p(4),p(3)],t)+polyval([1/p(2),p(1)],t)) ... # 2D + mp
  + p(20)*sin(2*polyval([1/p(4),p(3)],t)-polyval([1/p(8),p(7)],t)) ... # 2D - M
  + p(21)*sin(polyval([1/p(8),p(7)],t)-polyval([1/p(2),p(1)],t)) ... # M - mp
  + p(22)*sin(polyval([1/p(8),p(7)],t)+polyval([1/p(2),p(1)],t)) ... # M + mp
  + p(23)*sin(polyval([1/p(2),p(1)],t)+2*polyval([1/p(6),p(5)],t)) ... # mp + 2F
  + p(24)*sin(polyval([1/p(2),p(1)],t)-2*polyval([1/p(6),p(5)],t)) ... # mp - 2F
  + p(25)*sin(4*polyval([1/p(4),p(3)],t)-polyval([1/p(2),p(1)],t)) ... # 4D - mp
  + p(26)*sin(3*polyval([1/p(2),p(1)],t)) ... # 3mp
  + p(27)*sin(4*polyval([1/p(4),p(3)],t)-2*polyval([1/p(2),p(1)],t)) ... # 4D-2mp
  + p(28)*sin(polyval([1/p(6),p(5)],t)) ... # F
  + p(29)*sin(polyval([1/p(2),p(1)],t)+polyval([1/p(6),p(5)],t)) ... # mp+F
  + p(30)*sin(polyval([1/p(2),p(1)],t)-polyval([1/p(6),p(5)],t)) ... # mp-F
  + p(31)*sin(2*polyval([1/p(4),p(3)],t)-polyval([1/p(8),p(7)],t)+polyval([1/p(2),p(1)],t)) ... # 2D - M + mp
  + p(32)*sin(3*polyval([1/p(4),p(3)],t)-2*polyval([1/p(6),p(5)],t)) ... # 3D - 2F
  + p(33)*sin(polyval([1/p(4),p(3)],t)+2*polyval([1/p(6),p(5)],t)) ... # D + 2F
  + p(34)*sin(polyval([1/p(4),p(3)],t)-2*polyval([1/p(6),p(5)],t)) ... # D - 2F
  + p(35)*sin(2*polyval([1/p(4),p(3)],t)+polyval([1/p(8),p(7)],t)) ... # 2D + M
  + p(36)*sin(polyval([1/p(4),p(3)],t)-polyval([1/p(2),p(1)],t)) ... # D - mp
  + p(37)*sin(polyval([1/p(4),p(3)],t)+polyval([1/p(8),p(7)],t)) ... # D + M
  + p(38)*sin(2*polyval([1/p(4),p(3)],t)+polyval([1/p(8),p(7)],t)-polyval([1/p(2),p(1)],t)) ... # 2D + M - mp
  + p(39)*sin(2*(polyval([1/p(4),p(3)],t)+polyval([1/p(2),p(1)],t))) ... # 2D + 2mp
  + p(40)*sin(4*polyval([1/p(4),p(3)],t)) ... # 4D
  + p(41)*sin(2*polyval([1/p(4),p(3)],t)-3*polyval([1/p(2),p(1)],t)) ... # 2D - 3mp
  + p(42)*sin(polyval([1/p(8),p(7)],t)-2*polyval([1/p(2),p(1)],t)) ... # M - 2mp
  + p(43)*sin(6*polyval([1/p(4),p(3)],t)) ... # 6D
  ))/2 - phase));

fitfunc(args)

[result,fval]=fminsearch(fitfunc,args);
result
fval

