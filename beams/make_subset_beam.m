clear all;
close all;

load('beam_start_of_FFS_ey20nm_1e5.dat');
B = beam;

% select a random subset:
nSample = 10000;
nRows = max(size(B))
rndIDX = randperm(nRows); 
B10000 = B(rndIDX(1:nSample), :); 

% select a new random subset:
nSample = 20000;
nRows = max(size(B))
rndIDX = randperm(nRows); 
B20000 = B(rndIDX(1:nSample), :); 


[size(B); size(B20000); size(B10000)]
[mean(B); mean(B20000); mean(B10000)]
[std(B); std(B20000); std(B10000)]

figure
plot(B(:,4), B(:,1), 'b.', B20000(:,4), B20000(:,1), 'rx', B10000(:,4), B10000(:,1), 'gs') 

# save beams
beam = B20000;
save beam_start_of_FFS_ey20nm_2e4.dat beam;

beam = B10000;
save beam_start_of_FFS_ey20nm_1e4.dat beam;

