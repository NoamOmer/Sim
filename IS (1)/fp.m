
function fp(vec)

vec = squeeze(vec);
if (ndims(vec) > 1)
	warning('Warning: fp does not support plotting nD vectors');
end;
	
figure; plot(vec,'.-');

return;
