
Ntrials=100;

"call parpool"

tic
numCores = feature('numcores')
p = parpool(numCores);
toc

Msizes=[4 10 20 50];

for M_i = 1:length(Msizes)
    All_for_times(M_i)=skyrmions_fortest(Msizes(M_i));
    All_parfor_times(M_i)=skyrmions_parfortest(Msizes(M_i));
end

plot(Msizes,All_for_times,Msizes,All_parfor_times);
drawnow
saveas(gcf,"forvpar_times.fig")
loglog(Msizes,All_for_times,Msizes,All_parfor_times);
drawnow
saveas(gcf,"forvpar_log_times.fig")

delete(p);
