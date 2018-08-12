% 2 functions file

% func_t_class = struct('one',@func_one,...
%                      'two',@func_two);

% fh_array = cellfun(@func_one, @func_two);
function mother_func(num)
if (num == 1)
eval('function func_one(void); a=5; return;');
else
eval('function func_one(void); a=5; return;');
end;
return;

%function func_one(void)
%disp('func one');
%return;

%function func_two(void)
%disp('func two');
%return;

