function [m,n] = Normalization1s(m,mins,maxs)

    for i=1:size(m,2)

            m(:,i) = (m(:,i)-mins(i))/(maxs(i)-mins(i)+10e-8);

    end

