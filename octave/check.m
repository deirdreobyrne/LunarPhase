format long;

data = csvread("phase.csv");
n=524288;
phase = data([1:n],2);

data = csvread("dlambda.csv");
i = data([1:n],2);
i_phase=(1-cos(i))/2;

abs(max(i_phase-phase))

