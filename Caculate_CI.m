function CI = Caculate_CI(Coef)


Coeft=[];
ind = [];

if size(Coef,1)>1
    for i=1:size(Coef,1)
        Coefs = Coef(i,:,:);
%         ind=find(abs(Coefs)>10*median(abs(Coefs)));
        Coefs(:,:,ind) = [];
        CI(i,:) = reshape(prctile(Coefs , [2.5, 97.5]),1,2);
    end

elseif size(Coef,2)>1
    
    for i=1:size(Coef,2)
        Coefs = Coef(:,i,:);
%         ind=find(abs(Coefs)>2*median(abs(Coefs)));
        Coefs(:,:,ind) = [];
        CI(i,:) = reshape(prctile(Coefs , [2.5, 97.5])',1,2);

    end


end