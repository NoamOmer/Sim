
function [pulse,N] = read_siemens_PTA_pulse(fn,plot_flag)

pulse = [];
if (exist(fn,'file') ~= 2)
	return;
end;

if (strcmp(fn(end-2:end),'pta') == 1)          % Siemens format file
	fid = fopen(fn);
	for idx = 1:10
		ln  = fgetl(fid);
	end;

	while (~isempty(ln) && isstr(ln))
		p = sscanf(ln,'%f %f ; %f');
		pulse(end+1) = p(1)*exp(1i*p(2));
		ln  = fgetl(fid);
	end;
	fclose(fid);

elseif (strcmp(fn(end-2:end),'txt') == 1)      % Siemens file, manually exported
	fid = fopen(fn);
	p = fscanf(fid,'%f %f',inf);
	fclose(fid);
	pulse = p(1:2:end-1) .* exp(1i*p(2:2:end));
else
	error('Unknown file format... Exiting');
end;

N = length(pulse);

if exist('plot_flag','var') && plot_flag
	figure;
	subplot(211); plot(  abs(pulse),'.-'); title('Pulse Amplitude');
	subplot(212); plot(angle(pulse),'.-'); title('Pulse phase [rad]');
end;

return;

