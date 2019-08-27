function eve = Estimation_boot(input)


radio = input(:,1);
labeltotal = input(:,3);
cutoff = unique(input(:,2));

pre = radio;
true = labeltotal;
% [~,cutoff,~,~] = AllAuc (pre,true);

ytest = true; ytest(ytest<=0)=-1;

rangepr = unique(pre);
if length(rangepr)==2
   yd = pre; yd(pre<=0)=-1;yd(pre>0)=1;
else
    yd = pre; yd(pre<=cutoff)=-1;yd(pre>cutoff)=1;
end




B= ytest-yd;
tp = length(find((ytest==1)&(yd==1)));
tptn=length(find(~B));
tn = tptn-tp;
fp = length(find(B==-2));
fn = length(find(B==2));

correctrate(1) = (tptn/(length(ytest)+ 10e-7))*100;
correctrate(2) = tp/(tp+fp+10e-7);

record = zeros(length(ytest),2);
record(:,1) = yd;
record(:,2) = ytest;
[~,~,~, auc]=perfcurve(record(:,2),pre, 1);

fpr = fp /(fp + tn + 10e-7);
fnr = fn /(tp + fn + 10e-7); 
tpr = tp /(tp + fn + 10e-7);
tnr = tn / (tn + fp + 10e-7);

accuracy = correctrate(:,1);
precision = correctrate(:,2);
Fscore = 2*tp/(2*tp+fp+fn+10e-7);

% eve = [accuracy auc fpr fnr tpr tnr];

eve = [accuracy auc tpr/fpr fnr/tnr tpr tnr];