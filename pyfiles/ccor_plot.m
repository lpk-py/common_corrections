load lat.mat
load lon.mat
load tau.mat
scatter3(lon, lat, tau, 20, tau, 'filled')
view(0, 90)
colorbar

load ecorr
figure;
scatter3(lon, lat, ecorr, 20, ecorr, 'filled')
view(0, 90)
colorbar
