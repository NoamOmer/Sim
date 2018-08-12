% Loads an FID ASCII file and converts it into a complex vector  
%
% DO NOT ERASE LOG; ADD BELOW LOWEST LINE.
%
% Written By Yoav Shrot,     , July    2003
% Revised by Noam Ben-Eliezer, October 2006,
%  - Receive 'fullname' as an input parameter
%  - Commented function at the bottom is from BoazS and enables chooing between 2 different file formats.

function [DATA] = read_ASCII_fid(fullname)

fd = fopen(fullname,'r','native');
% if (fd <= 0) disp(sprintf('fd=%3.3f, fn=%s\n',fd,fullname)); end;
DATA_temp = transpose(fscanf(fd, '%f %f', [2,inf]));
fclose(fd);
DATA = DATA_temp(:,1) + i*DATA_temp(:,2);

return;

% function [DATA] = read_ASCII_fid(fullname)
% x=load(fullname);
% [a,b]=size(x);
% if  b==1
%     DATA = x(1:2:end)+i*x(2:2:end);
% elseif b==2
%     DATA = x(1:1:end,1)+i*x(1:1:end,2);
% end

