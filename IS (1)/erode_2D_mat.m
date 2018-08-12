
function [out_mat] = erode_2D_mat(in_mat,nPx)

if (nPx <= 0)
	out_mat = in_mat;
	return;
end;

m = in_mat;
m(m~=0) = 1;

for idx = 1:nPx
	m1 = m;  m1(1:end-1,:) = m1(2:end  ,:);
	m2 = m;  m2(2:end  ,:) = m2(1:end-1,:);
	m3 = m;  m3(:,1:end-1) = m3(:,2:end  );
	m4 = m;  m4(:,2:end  ) = m4(:,1:end-1);

	m = m1 + m2 + m3 + m4;
	m(m  < 4) = 0;
	m(m ~= 0) = 1;
end;
out_mat = in_mat.*m;

return;

