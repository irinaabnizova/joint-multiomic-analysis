function [second_max]=second_max(vec)
%vec=[1,2,3]
ma=max(vec);
vec1=vec;
for i=1:length(vec)
    if vec1(i)==ma,
        vec1(i)=min(vec);
    end
end

second_max=max(vec1);