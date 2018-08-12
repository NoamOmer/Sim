function change_order_of_object(o,action)
p=get(o,'parent');
cp=get(p,'child');
k=find(cp==o);
if size(cp,2)<size(cp,1)
    cp=cp';
end

switch lower(action)
case 'bottom'
    cp(k)=[];
    cp=[cp o];
case 'front'
    cp(k)=[];
    cp=[o cp];
end
set(p,'child',cp);
    