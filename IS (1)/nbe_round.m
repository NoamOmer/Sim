function rounded_num = nbe_round(num,prec)
	if (prec == 0)
		rounded_num = num;
	else
		rounded_num = roundn(num,-prec);
	end;
return;

