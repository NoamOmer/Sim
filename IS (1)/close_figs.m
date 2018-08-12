function close_figs(void)
for idx = 1:100
	if ishandle(idx)
		close(idx);
	end;
end;
return;

