clear 
clc
close all
load('allfeatures.mat');

trainfeature =[trainfeature trainclincal];
testfeature =[testfeature testclincal]; 
valfeature = [valfeature valclincal]; 
labeltrain = trainclincal(:,end-1);
labeltest = testclincal(:,end-1);
labelval = valclincal(:,end-1);

load( 'models.mat');
load('trainnorm.mat');
featurenum=790;
trainfeature (:,1:featurenum) = Normalization1s(trainfeature(:,1:featurenum),minms,maxms);
testfeature(:,1:featurenum)  = Normalization1s(testfeature(:,1:featurenum),minms,maxms);
valfeature(:,1:featurenum)  = Normalization1s(valfeature(:,1:featurenum),minms,maxms);


indtrain1 = find(trainclincal(:,3)==1);
indtrain2 = find(trainclincal(:,3)==2);
trainfeature1 = trainfeature(indtrain1,:);
trainfeature2 = trainfeature(indtrain2,:);

indtest1 = find(testclincal(:,3)==1);
indtest2 = find(testclincal(:,3)==2);
testfeature1 = testfeature(indtest1,:);
testfeature2 = testfeature(indtest2,:);



indval1 = find(valclincal(:,3)==1);
indval2 = find(valclincal(:,3)==2);
valfeature1 = valfeature(indval1,:);
valfeature2 = valfeature(indval2,:);

radiotrainpet1 =ones(size(trainfeature1,1),1)*modelcoefpet1(end,2);
radiotestpet1 = ones(size(testfeature1,1),1)*modelcoefpet1(end,2);
radiovalpet1 = ones(size(valfeature1,1),1)*modelcoefpet1(end,2);
for i = 1:size(modelcoefpet1,1)-1
    radiotrainpet1 = radiotrainpet1 + modelcoefpet1(i,2)*trainfeature1(:,modelcoefpet1(i,1));
    radiotestpet1 = radiotestpet1 + modelcoefpet1(i,2)*testfeature1(:,modelcoefpet1(i,1));
    radiovalpet1 = radiovalpet1 + modelcoefpet1(i,2)*valfeature1(:,modelcoefpet1(i,1));
end

radiotrainpet2 =ones(size(trainfeature2,1),1)*modelcoefpet2(end,2);
radiotestpet2 = ones(size(testfeature2,1),1)*modelcoefpet2(end,2);
radiovalpet2 = ones(size(valfeature2,1),1)*modelcoefpet2(end,2);
for i = 1:size(modelcoefpet2,1)-1
    radiotrainpet2 = radiotrainpet2 + modelcoefpet2(i,2)*trainfeature2(:,modelcoefpet2(i,1));
    radiovalpet2 = radiovalpet2 + modelcoefpet2(i,2)*valfeature2(:,modelcoefpet2(i,1));
    radiotestpet2 = radiotestpet2 + modelcoefpet2(i,2)*testfeature2(:,modelcoefpet2(i,1));
end

radiotrainct1 =ones(size(trainfeature1,1),1)*modelcoefct1(end,2);
radiotestct1 =ones(size(testfeature1,1),1)*modelcoefct1(end,2);
radiovalct1 =ones(size(valfeature1,1),1)*modelcoefct1(end,2);
for i = 1:size(modelcoefct1,1)-1
    radiotrainct1 = radiotrainct1 + modelcoefct1(i,2)*trainfeature1(:,modelcoefct1(i,1));
    radiotestct1 = radiotestct1 + modelcoefct1(i,2)*testfeature1(:,modelcoefct1(i,1));
    radiovalct1 = radiovalct1 + modelcoefct1(i,2)*valfeature1(:,modelcoefct1(i,1));
end

radiotrainct2 =ones(size(trainfeature2,1),1)*modelcoefct2(end,2);
radiotestct2 =ones(size(testfeature2,1),1)*modelcoefct2(end,2);
radiovalct2 =ones(size(valfeature2,1),1)*modelcoefct2(end,2);
for i = 1:size(modelcoefct2,1)-1
    radiotrainct2 = radiotrainct2 + modelcoefct2(i,2)*trainfeature2(:,modelcoefct2(i,1));
    radiotestct2 = radiotestct2 + modelcoefct2(i,2)*testfeature2(:,modelcoefct2(i,1));
    radiovalct2 = radiovalct2 + modelcoefct2(i,2)*valfeature2(:,modelcoefct2(i,1));
end

radiotrainpetct1 =ones(size(trainfeature1,1),1)*modelcoefpetct1(end,2);
radiotestpetct1 = ones(size(testfeature1,1),1)*modelcoefpetct1(end,2);
radiovalpetct1 = ones(size(valfeature1,1),1)*modelcoefpetct1(end,2);
for i = 1:size(modelcoefpetct1,1)-1
    radiotrainpetct1 = radiotrainpetct1 + modelcoefpetct1(i,2)*trainfeature1(:,modelcoefpetct1(i,1));
    radiotestpetct1 = radiotestpetct1 + modelcoefpetct1(i,2)*testfeature1(:,modelcoefpetct1(i,1));
     radiovalpetct1 = radiovalpetct1 + modelcoefpetct1(i,2)*valfeature1(:,modelcoefpetct1(i,1));
end

radiotrainpetct2 =ones(size(trainfeature2,1),1)*modelcoefpetct2(end,2);
radiotestpetct2 = ones(size(testfeature2,1),1)*modelcoefpetct2(end,2);
radiovalpetct2 = ones(size(valfeature2,1),1)*modelcoefpetct2(end,2);
for i = 1:size(modelcoefpetct2,1)-1
    radiotrainpetct2 = radiotrainpetct2 + modelcoefpetct2(i,2)*trainfeature2(:,modelcoefpetct2(i,1));
    radiotestpetct2 = radiotestpetct2 + modelcoefpetct2(i,2)*testfeature2(:,modelcoefpetct2(i,1));
    radiovalpetct2 = radiovalpetct2 + modelcoefpetct2(i,2)*valfeature2(:,modelcoefpetct2(i,1));
end


radiotrainfuse1 =ones(size(trainfeature1,1),1)*modelcoeffuse1(end,2);
radiotestfuse1 = ones(size(testfeature1,1),1)*modelcoeffuse1(end,2);
radiovalfuse1 = ones(size(valfeature1,1),1)*modelcoeffuse1(end,2);
for i = 1:size(modelcoeffuse1,1)-1
    radiotrainfuse1 = radiotrainfuse1 + modelcoeffuse1(i,2)*trainfeature1(:,modelcoeffuse1(i,1));
    radiotestfuse1 = radiotestfuse1 + modelcoeffuse1(i,2)*testfeature1(:,modelcoeffuse1(i,1));
    radiovalfuse1 = radiovalfuse1 + modelcoeffuse1(i,2)*valfeature1(:,modelcoeffuse1(i,1));
end


radiotrainfuse2 =ones(size(trainfeature2,1),1)*modelcoeffuse2(end,2);
radiotestfuse2 = ones(size(testfeature2,1),1)*modelcoeffuse2(end,2);
radiovalfuse2 = ones(size(valfeature2,1),1)*modelcoeffuse2(end,2);
for i = 1:size(modelcoeffuse2,1)-1
    radiotrainfuse2 = radiotrainfuse2 + modelcoeffuse2(i,2)*trainfeature2(:,modelcoeffuse2(i,1));
    radiotestfuse2 = radiotestfuse2 + modelcoeffuse2(i,2)*testfeature2(:,modelcoeffuse2(i,1));
    radiovalfuse2 = radiovalfuse2 + modelcoeffuse2(i,2)*valfeature2(:,modelcoeffuse2(i,1));
end



traindata1 = [ radiotrainpet1 radiotrainct1 radiotrainpetct1 radiotrainfuse1   trainfeature1(:,featurenum+1:end);];
traindata2 = [ radiotrainpet2 radiotrainct2 radiotrainpetct2  radiotrainfuse2 trainfeature2(:,featurenum+1:end);];
testdata1 =  [ radiotestpet1 radiotestct1 radiotestpetct1 radiotestfuse1  testfeature1(:,featurenum+1:end);];
testdata2 =  [ radiotestpet2 radiotestct2 radiotestpetct2  radiotestfuse2  testfeature2(:,featurenum+1:end);];
valdata1 =  [ radiovalpet1 radiovalct1 radiovalpetct1 radiovalfuse1  valfeature1(:,featurenum+1:end);];
valdata2 =  [ radiovalpet2 radiovalct2 radiovalpetct2  radiovalfuse2  valfeature2(:,featurenum+1:end);];

traindata = [traindata1;traindata2];
testdata = [testdata1;testdata2];
valdata = [valdata1;valdata2];
valdata(indval1,:)=valdata1;
valdata(indval2,:)=valdata2;



[X1,Y1,~,auc1]=perfcurve(traindata(:,end-1),traindata(:,1), 1);
[X2,Y2,~,auc2]=perfcurve(traindata(:,end-1),traindata(:,2), 1);
[X3,Y3,~,auc3]=perfcurve(traindata(:,end-1),traindata(:,3), 1);
[X4,Y4,~,auc4]=perfcurve(traindata(:,end-1),traindata(:,4), 1);

[X1t,Y1t,~,auc1t]=perfcurve(testdata(:,end-1),testdata(:,1), 1);
[X2t,Y2t,~,auc2t]=perfcurve(testdata(:,end-1),testdata(:,2), 1);
[X3t,Y3t,~,auc3t]=perfcurve(testdata(:,end-1),testdata(:,3), 1);
[X4t,Y4t,~,auc4t]=perfcurve(testdata(:,end-1),testdata(:,4), 1);

[X1v,Y1v,~,auc1v]=perfcurve(valdata(:,end-1),valdata(:,1), 1);
[X2v,Y2v,~,auc2v]=perfcurve(valdata(:,end-1),valdata(:,2), 1);
[X3v,Y3v,~,auc3v]=perfcurve(valdata(:,end-1),valdata(:,3), 1);
[X4v,Y4v,~,auc4v]=perfcurve(valdata(:,end-1),valdata(:,4), 1);

[~,cutoff]  = AllAuc(traindata(:,1:4),labeltrain);
[~,cutofft] = AllAuc(testdata(:,1:4),labeltest);
[~,cutoffv] = AllAuc(valdata(:,1:4),labelval);
cutofft=cutoff;
% Eve1 = EvaluationModel(traindata(:,1),traindata(:,end-1),1,cutoff(1));
% Eve1t = EvaluationModel(testdata(:,1),testdata(:,end-1),1,cutofft(1));
% Eve2= EvaluationModel(traindata(:,2),traindata(:,end-1),1,cutoff(2));
% Eve2t = EvaluationModel(testdata(:,2),testdata(:,end-1),1,cutofft(2));
% Eve3 = EvaluationModel(traindata(:,3),traindata(:,end-1),1,cutoff(3));
% Eve3t = EvaluationModel(testdata(:,3),testdata(:,end-1),1,cutofft(3));
% Eve4 = EvaluationModel(traindata(:,4),traindata(:,end-1),1,cutoff(4));
% Eve4t = EvaluationModel(testdata(:,4),testdata(:,end-1),1,cutofft(4));
% Eve1v = EvaluationModel(valdata(:,1),valdata(:,end-1),1,cutoffv(1));
% Eve2v= EvaluationModel(valdata(:,2),valdata(:,end-1),1,cutoffv(2));
% Eve3v = EvaluationModel(valdata(:,3),valdata(:,end-1),1,cutoffv(3));
% Eve4v = EvaluationModel(valdata(:,4),valdata(:,end-1),1,cutoffv(4));

figure,plot(X1,Y1,'k',X1t,Y1t,'r',X1v,Y1v,'b',X1,X1,'k-.','LineWidth',1.5);legend(['Training cohort',sprintf('\n') ,'AUC = ',num2str(roundn(auc1,-2))],['Test cohort',sprintf('\n') ,'AUC = ',num2str(roundn(auc1t,-2)) ],['Prospective Test cohort',sprintf('\n') ,'AUC = ',num2str(roundn(auc1v,-2)) ],'Reference line', 'Location','SouthEast');xlabel('1-Specificity');ylabel('Sensitivity');title('ROC curves of PETRS');set(gca,'FontSize',14)
figure,plot(X2,Y2,'k',X2t,Y2t,'r',X2v,Y2v,'b',X1,X1,'k-.','LineWidth',1.5);legend(['Training cohort',sprintf('\n') ,'AUC = ',num2str(roundn(auc2,-2))],['Test cohort',sprintf('\n') ,'AUC = ',num2str(roundn(auc2t,-2))],['Prospective Test cohort',sprintf('\n') ,'AUC = ',num2str(roundn(auc2v,-2)) ],'Reference line', 'Location','SouthEast');xlabel('1-Specificity');ylabel('Sensitivity');title('ROC curves of CTRS');set(gca,'FontSize',14)
figure,plot(X3,Y3,'k',X3t,Y3t,'r',X3v,Y3v,'b',X1,X1,'k-.','LineWidth',1.5);legend(['Training cohort',sprintf('\n') ,'AUC = ',num2str(roundn(auc3,-2))],['Test cohort',sprintf('\n') ,'AUC = ',num2str(roundn(auc3t,-2)) ],['Prospective Test cohort',sprintf('\n') ,'AUC = ',num2str(roundn(auc3v,-2)) ],'Reference line', 'Location','SouthEast');xlabel('1-Specificity');ylabel('Sensitivity');title('ROC curves of PETCTRS');set(gca,'FontSize',14)
figure,plot(X4,Y4,'k',X4t,Y4t,'r',X4v,Y4v,'b',X1,X1,'k-.','LineWidth',1.5);legend(['Training cohort',sprintf('\n') ,'AUC = ',num2str(roundn(auc4,-2))],['Test cohort',sprintf('\n') ,'AUC = ',num2str(roundn(auc4t,-2))],['Prospective Test cohort',sprintf('\n') ,'AUC = ',num2str(roundn(auc4v,-2)) ],'Reference line', 'Location','SouthEast');xlabel('1-Specificity');ylabel('Sensitivity');title('ROC curves of mpRS on all patients');set(gca,'FontSize',14)


echonum=1000;
everadio1 = bstrap(echonum,1,'Estimation_boot',[traindata(:,1) repmat(cutoff(1),size(traindata,1),1) traindata(:,end-1) ]);
everadio2 = bstrap(echonum,1,'Estimation_boot',[traindata(:,2) repmat(cutoff(2),size(traindata,1),1) traindata(:,end-1) ]);
everadio3 = bstrap(echonum,1,'Estimation_boot',[traindata(:,3) repmat(cutoff(3),size(traindata,1),1) traindata(:,end-1) ]);
everadio4 = bstrap(echonum,1,'Estimation_boot',[traindata(:,4) repmat(cutoff(4),size(traindata,1),1) traindata(:,end-1) ]);
everadiom1 =  [Eve1' Caculate_CI(everadio1)];
everadiom2 =  [Eve2' Caculate_CI(everadio2)];
everadiom3 =  [Eve3' Caculate_CI(everadio3)];
everadiom4 =  [Eve4' Caculate_CI(everadio4)];
everadiot1 = bstrap(echonum,1,'Estimation_boot',[testdata(:,1) repmat(cutoff(1),size(testdata,1),1) testdata(:,end-1) ]);
everadiot2 = bstrap(echonum,1,'Estimation_boot',[testdata(:,2) repmat(cutoff(2),size(testdata,1),1) testdata(:,end-1) ]);
everadiot3 = bstrap(echonum,1,'Estimation_boot',[testdata(:,3) repmat(cutoff(3),size(testdata,1),1) testdata(:,end-1) ]);
everadiot4 = bstrap(echonum,1,'Estimation_boot',[testdata(:,4) repmat(cutoff(4),size(testdata,1),1) testdata(:,end-1) ]);
everadiomt1 =  [Eve1t' Caculate_CI(everadiot1)];
everadiomt2 =  [Eve2t' Caculate_CI(everadiot2)];
everadiomt3 =  [Eve3t' Caculate_CI(everadiot3)];
everadiomt4 =  [Eve4t' Caculate_CI(everadiot4)];
everadiov1 = bstrap(echonum,1,'Estimation_boot',[valdata(:,1) repmat(cutoffv(1),size(valdata,1),1) valdata(:,end-1) ]);
everadiov2 = bstrap(echonum,1,'Estimation_boot',[valdata(:,2) repmat(cutoffv(2),size(valdata,1),1) valdata(:,end-1) ]);
everadiov3 = bstrap(echonum,1,'Estimation_boot',[valdata(:,3) repmat(cutoffv(3),size(valdata,1),1) valdata(:,end-1) ]);
everadiov4 = bstrap(echonum,1,'Estimation_boot',[valdata(:,4) repmat(cutoffv(4),size(valdata,1),1) valdata(:,end-1) ]);
everadiomv1 =  [Eve1v' Caculate_CI(everadiov1)];
everadiomv2 =  [Eve2v' Caculate_CI(everadiov2)];
everadiomv3 =  [Eve3v' Caculate_CI(everadiov3)];
everadiomv4 =  [Eve4v' Caculate_CI(everadiov4)];


numtrain = length(labeltrain);
numtest = length(labeltest);
numval = length(labelval);

grouptrainpet = zeros(numtrain,1);grouptrainct = zeros(numtrain,1);grouptrainpetct = zeros(numtrain,1);grouptrainfuse = zeros(numtrain,1);
grouptestpet = zeros(numtest,1);grouptestct = zeros(numtest,1);grouptestpetct = zeros(numtest,1);grouptestfuse = zeros(numtest,1);
groupvalpet = zeros(numval,1);groupvalct = zeros(numval,1);groupvalpetct = zeros(numval,1);groupvalfuse = zeros(numval,1);

n=2;

grouptrainpet(traindata(:,1)>=cutoff(1))=1;
grouptrainct(traindata(:,2)>=cutoff(2))=1;
grouptrainpetct(traindata(:,3)>=cutoff(3))=1;
grouptrainfuse(traindata(:,4)>=cutoff(4))=1;

grouptestpet(testdata(:,1)>=cutoff(1))=1;
grouptestct(testdata(:,2)>=cutoff(2))=1;
grouptestpetct(testdata(:,3)>=cutoff(3))=1;
grouptestfuse(testdata(:,4)>=cutoff(4))=1;

groupvalpet(valdata(:,1)>=cutoff(1))=1;
groupvalct(valdata(:,2)>=cutoff(2))=1;
groupvalpetct(valdata(:,3)>=cutoff(3))=1;
groupvalfuse(valdata(:,4)>=cutoff(4))=1;


traingroup = [ grouptrainpet grouptrainct grouptrainpetct grouptrainfuse  trainfeature(:,featurenum+1:end);];
testgroup = [ grouptestpet grouptestct grouptestpetct grouptestfuse  testfeature(:,featurenum+1:end);];
valgroup = [ groupvalpet groupvalct groupvalpetct groupvalfuse  valfeature(:,featurenum+1:end);];

indh1 = traingroup(:,4)==1;
indh2 = traingroup(:,4)==0;

X1 = [traingroup(indh1,21)/30 1-traingroup(indh1,22)];
X2 = [traingroup(indh2,21)/30 1-traingroup(indh2,22)];
figure, logrank(X1,X2)

X1 = [traingroup(indh1,24)/30 1-traingroup(indh1,25)];
X2 = [traingroup(indh2,24)/30 1-traingroup(indh2,25)];
figure, logrank(X1,X2)


indh1 = testgroup(:,4)==1;
indh2 = testgroup(:,4)==0;

X1 = [testgroup(indh1,21)/30 1-testgroup(indh1,22)];
X2 = [testgroup(indh2,21)/30 1-testgroup(indh2,22)];
figure, logrank(X1,X2)

X1 = [testgroup(indh1,24)/30 1-testgroup(indh1,25)];
X2 = [testgroup(indh2,24)/30 1-testgroup(indh2,25)];
figure, logrank(X1,X2)


indh1 = valgroup(:,4)==1;
indh2 = valgroup(:,4)==0;

X1 = [valgroup(indh1,21)/30 1-valgroup(indh1,22)];
X2 = [valgroup(indh2,21)/30 1-valgroup(indh2,22)];
figure, logrank(X1,X2)

X1 = [valgroup(indh1,24)/30 1-valgroup(indh1,25)];
X2 = [valgroup(indh2,24)/30 1-valgroup(indh2,25)];
figure, logrank(X1,X2)





