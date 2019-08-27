function [allfeatureauc cutoff flag2 posind negind] = AllAuc (allfeature,label,flag)


if nargin<3
    flag=0;
    if nargin<2
        label = allfeature(:,end);
        allfeature = allfeature(:,1:end-1);
    end

end
allfeatureauc = zeros(size(allfeature,2),1);
cutoff = zeros(size(allfeature,2),1);
posind = [];
negind = [];


if flag==0
    
    for fn = 1: size(allfeature,2)
%         fn
       feature = allfeature(:,fn);
       if sum(abs(feature))&length(unique(feature))>1
           try
            [allfeatureauc(fn) cutoff(fn) flag2 ] = SingleROC(feature,label);
           catch
             allfeatureauc(fn)=0;flag2=0;
           end

       else
           allfeatureauc(fn) = 0.5;
            cutoff(fn)= 0;
            flag2=3;
       end
       if flag==1
           posind = [posind;fn];
       else
           negind = [negind;fn];
       end

    end
else
    figure;
    hold on;
    colornum = size(allfeature,2);
    for fn = 1: size(allfeature,2)

       feature = allfeature(:,fn);
       if sum(abs(feature))&length(unique(feature))>1
            try
            [allfeatureauc(fn) cutoff(fn) flag2 ] = SingleROC(feature,label);
           catch
             allfeatureauc(fn)=0;
           end
       else
           allfeatureauc(fn) = 0.5;
            cutoff(fn)= 0;
       end
       if flag==1
           posind = [posind;fn];
       else
           negind = [negind;fn];
       end
       
       switch mod(fn,6)
           case 0
               color = [0 0 fn/colornum];
           case 1
               color = [0 fn/colornum 0];
           case 2
               color = [fn/colornum 0 0];
           case 3
               color = [0 fn/colornum fn/colornum];
           case 4
               color = [fn/colornum 0 fn/colornum];
           case 5
               color = [fn/colornum fn/colornum 0];
       end
       
       plot(xaxis, yaxis,'Color',color);
       

    end
    hold off
end