%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Check vertical grid size of ROMS-1D
%
% Patrick Marchesiello - June 2012
%
% Modified to output grid for GOTM
% Qing Li, 20180509
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

N=100;      % number of model levels
theta_s=2;  % surface stretching parameter
hmax=200;   % model depth

ods=N;ds=1/ods;
cff=hmax/sinh(theta_s);
for k=N:-1:1
  sc_w=ds*(k-N);
  z_w(k+1)=-cff*sinh(-theta_s*sc_w);
  sc_r=ds*(k-N-0.5);
  z_r(k)=-cff*sinh(-theta_s*sc_r);
end
z_w(1)=-hmax;
for k=1:N
  Hz(k)=z_w(k+1)-z_w(k);
end

figure
plot(Hz,z_r,'*');
title('Vertical grid size dz (m)')
xlabel('dz (m)')
ylabel('depth (m)')

% write to file
fout = ['dz_ROMS_D' sprintf('%d', hmax) 'N' sprintf('%d',N) '.dat'];
fileID = fopen(fout,'w');
fprintf(fileID,'%d\n', N);
for i=N:-1:1
    fprintf(fileID,'%12.8f\n', Hz(i));
end
fclose(fileID);

% test precision
fileID = fopen(fout, 'r');
Hz_r = textscan(fileID, '%12.8f\n', N, 'HeaderLines',1);
depth_r = sum(Hz_r{1});
fprintf('Total depth = %12.8e\n', depth_r);
fprintf('Diff. = %12.8e\n', depth_r-hmax);
fclose(fileID);