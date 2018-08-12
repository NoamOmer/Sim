function set_pointer(pointer)
if nargin==0
    pointer='arrow';
else
    S=set(gcf,'pointer');
    K=strcmp(lower(S),lower(pointer));
    if sum(K)==0
        switch lower(pointer)
            case 'getpoints'
                [pointerShape, pointerHotSpot] = CreatePointer;
                set(gcf,'pointer','custom','PointerShapeCData', pointerShape, ...
                    'PointerShapeHotSpot', pointerHotSpot);    
        end
        return
    end
end
set(gcf,'pointer',pointer)

function [pointerShape, pointerHotSpot] = CreatePointer
pointerHotSpot = [8 8];
pointerShape = [ ...
        NaN NaN NaN NaN NaN   1   2 NaN   2   1 NaN NaN NaN NaN NaN NaN
    NaN NaN NaN NaN NaN   1   2 NaN   2   1 NaN NaN NaN NaN NaN NaN
    NaN NaN NaN NaN NaN   1   2 NaN   2   1 NaN NaN NaN NaN NaN NaN
    NaN NaN NaN NaN NaN   1   2 NaN   2   1 NaN NaN NaN NaN NaN NaN
    NaN NaN NaN NaN NaN   1   2 NaN   2   1 NaN NaN NaN NaN NaN NaN
    1   1   1   1   1   1   2 NaN   2   1   1   1   1   1   1   1
    2   2   2   2   2   2   2 NaN   2   2   2   2   2   2   2   2
    NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    2   2   2   2   2   2   2 NaN   2   2   2   2   2   2   2   2
    1   1   1   1   1   1   2 NaN   2   1   1   1   1   1   1   1
    NaN NaN NaN NaN NaN   1   2 NaN   2   1 NaN NaN NaN NaN NaN NaN
    NaN NaN NaN NaN NaN   1   2 NaN   2   1 NaN NaN NaN NaN NaN NaN
    NaN NaN NaN NaN NaN   1   2 NaN   2   1 NaN NaN NaN NaN NaN NaN
    NaN NaN NaN NaN NaN   1   2 NaN   2   1 NaN NaN NaN NaN NaN NaN
    NaN NaN NaN NaN NaN   1   2 NaN   2   1 NaN NaN NaN NaN NaN NaN
    NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];
