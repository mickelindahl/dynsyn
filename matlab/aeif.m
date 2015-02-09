figure(1)
% SNR
g=5;
E_L=-53;
V_T=-53; % Spikeing above + a=3 gives -54
delta_T=2; % Gives critical voltage at -50 for which an action potantial is generatet
a=3;

f=@ (V) -g*(V-E_L)+g*delta_T.*exp((V-V_T)./delta_T)-a*(V-E_L);

v=-80:-45;
clf
plot(v,f(v))
hold all
plot(v,zeros(size(v)))
xlim([v(1),v(end)])


figure(2)
% GPE
g=3;
E_L=-53;
V_T=-53; % Spikeing above
delta_T=2; % Gives critical voltage at -43 for which an action potantial is generatet
a=1;

f=@ (V) -g*(V-E_L)+g*delta_T.*exp((V-V_T)./delta_T)-a*(V-E_L);

v=-80:-40;
clf
plot(v,f(v))
hold all
plot(v,zeros(size(v)))
xlim([v(1),v(end)])
