function [auc cutoff flag xaxis yaxis] = SingleROC (feature,label)
totalN = 1000;
opmin = 10e-9;

ll = sort(unique(label));
label(label==ll(1))=-1;
label(label==ll(2))=1;

[xaxis,yaxis,Th,auc] = perfcurve(label,feature,1) ;
cuts = yaxis-xaxis+1;
[maxtemp,ind] = max(cuts);
cutoff = Th(ind);

posf = feature(label==1);
negf = feature(label==-1);
if mean(posf)>=mean(negf)
    flag = 1;
else 
    flag=0;
end

% if flag==0
%     auc=1-auc;
% end

% if max(feature)>0
%     if max(feature)>=1
%         Digitnump = floor(log10(max(feature)))+1;
%     else
%         Digitnump = 0;
%     end
% else
%     Digitnump = 0;
% end
% if min(feature)<0
%     if min(feature)<=-1
%         Digitnumn = floor(log10(-min(feature)))+1;
%     else
%         Digitnumn = 0;
%     end
% else
%     Digitnumn = 0;
% end
% 
% if (Digitnump+Digitnumn)==0
%     number = totalN;
% else
%     number = totalN/(Digitnump+Digitnumn);
% end
% 
% 
% featurep = feature(feature>=0);
% featuren = feature(feature<0);
% 
% thresholdp = [];
% thresholdn = [];
% if (~isempty(featurep))
% 
%     if Digitnump==0
%         inc = (max(featurep)-min(featurep))/number;
%         thresholdp = [thresholdp min(featurep):inc:max(featurep)];
%     else  
%         if min(featurep)<1
%           for num = 0:Digitnump
%               if num==0
%                   minf = min(featurep);      
%               else
%                   minf = 10^(num-1);
%               end
%               
%               if num==Digitnump
%                   maxf = max(featurep);
%               else
%                   maxf = 10^num;
%               end
%                            
%               inc = (maxf - minf)/number;
%               thresholdp = [thresholdp minf:inc:maxf]; 
%               
%           end
%         else
%             for num = 1:Digitnump
%               if num==1
%                   minf = min(featurep);
%               else
%                   minf = 10^(num-1);
%               end
%               
%               if num==Digitnump
%                   maxf = max(featurep);
%               else
%                   maxf = 10^num;
%               end
%                            
%               inc = (maxf - minf)/number;
%               thresholdp = [thresholdp minf:inc:maxf]; 
%               
%             end
%         end
%     end
% end
% 
% 
% if (~isempty(featuren))
% 
%     if Digitnumn==0
%         inc = (max(featuren)-min(featuren))/number;
%         thresholdn = [thresholdn min(featuren):inc:max(featuren)];
%     else  
%         if max(featuren)>-1
%           for num = 0:Digitnumn
%               if num==0
%                   maxf = max(featuren);
%               else
%                   maxf = -10^(num-1);
%               end
%               
%               if num==Digitnumn
%                   minf = min(featuren);
%               else
%                   minf = -10^num;
%               end
%                            
%               inc = (maxf - minf)/number;
%               thresholdn = [thresholdn minf:inc:maxf]; 
%               
%           end
%         else
%             for num = 1:Digitnumn
%               if num==1
%                   maxf = max(featuren);
%               else
%                   maxf = -10^(num-1);
%               end
%               
%               if num==Digitnumn
%                   minf = min(featuren);
%               else
%                   minf = -10^num;
%               end
%                            
%               inc = (maxf - minf)/number;
%               thresholdn = [thresholdn minf:inc:maxf]; 
%               
%             end
%         end
%     end
% end
% 
% threshold =[thresholdn thresholdp];
% threshold = sort(threshold);
% 
% posf = feature(label==1);
% negf = feature(label==-1);
% 
% xaxis = zeros(length(threshold),1);
% yaxis = zeros(length(threshold),1);
% sensitivitys = zeros(length(threshold),1);
% specificitys = zeros(length(threshold),1);
% 
% 
% 
% 
% if mean(posf)>=mean(negf)
%     flag = 1;
%     y = zeros(length(feature),length(threshold));
%     thrs = repmat(threshold,length(feature),1);
%     features = repmat(feature,1,length(threshold));
%     if size(features)~=size(thrs)
%         disp('anyting');
%     end
%     try    
%     temp = features - thrs;
%     temp(abs(temp)<opmin)=0;
%     catch
%         disp('anyting');
%     end
%     y(temp>=0)=1;
%     y(temp<0)=-1;
% 
%     B = y - repmat(reshape(label,length(label),1),1,length(threshold));
%     labels =  repmat(reshape(label,length(label),1),1,length(threshold));
%     Btemp = B;
%     Btemp(Btemp~=0)=10;
%     Btemp(Btemp==0)=1;
%     Btemp(Btemp==10)=0;
%     Btemp2 = Btemp;
%     Btemp(labels~=1)=0;
%     Btemp2(labels~=-1)=0;
%     tp = sum(Btemp);
%     tn = sum(Btemp2);
% 
%     Btemp = B;
%     Btemp(Btemp~=2)=0;
%     Btemp(Btemp==2)=1;
%     fp = sum(Btemp);
%     Btemp = B;
%     Btemp(Btemp~=-2)=0;
%     Btemp(Btemp==-2)=1;
%     fn = sum(Btemp);
% 
%     sensitivitys = tp./(tp+fn+10e-7);
%     specificitys = tn./(tn+fp+10e-7);
%     xaxis = 1-specificitys;
%     yaxis = sensitivitys;
%     xaxis = [xaxis 0];
%     yaxis = [yaxis 0];
% 
% else
%     flag=0;
%     y = zeros(length(feature),length(threshold));
%     thrs = repmat(threshold,length(feature),1);
%     features = repmat(feature,1,length(threshold));
%     temp = features - thrs;
%     temp(abs(temp)<opmin)=0;
%     y(temp>=0)=-1;
%     y(temp<0)=1;
% 
%     B = y - repmat(reshape(label,length(label),1),1,length(threshold));
%     labels =  repmat(reshape(label,length(label),1),1,length(threshold));
%     Btemp = B;
%     Btemp(Btemp~=0)=10;
%     Btemp(Btemp==0)=1;
%     Btemp(Btemp==10)=0;
%     Btemp2 = Btemp;
%     Btemp(labels~=1)=0;
%     Btemp2(labels~=-1)=0;
%     tp = sum(Btemp);
%     tn = sum(Btemp2);
% 
%     Btemp = B;
%     Btemp(Btemp~=2)=0;
%     Btemp(Btemp==2)=1;
%     fp = sum(Btemp);
%     Btemp = B;
%     Btemp(Btemp~=-2)=0;
%     Btemp(Btemp==-2)=1;
%     fn = sum(Btemp);
% 
%     sensitivitys = tp./(tp+fn+10e-7);
%     specificitys = tn./(tn+fp+10e-7);
%     xaxis = 1-specificitys;
%     yaxis = sensitivitys;
% 
%     xaxis = [xaxis 1];
%     yaxis = [yaxis 1];
% 
%     
% end
% 
% 
% % xaxis = fliplr(xaxis);
% % yaxis = fliplr(yaxis);
% auctemp = abs(trapz(xaxis,yaxis));
% % [~,~,~, auctemp]=perfcurve(label,feature, 1);
% 
% % [X,Y,T,auc] = perfcurve(label,feature,1) ;
% 
% cuts = sensitivitys+specificitys;
% [maxtemp,ind] = max(cuts);
% 
% if ~isempty(cuts)
%    cutofftemp= threshold(ind); 
% else
%     cutofftemp = min(feature);
% end
% 
% auc = auctemp;
% cutoff = cutofftemp;


% 
% %%%%%%%%%%%%%%%%%%%%%%
% newmaxf = cutofftemp + inc/2;
% newminf = cutofftemp-inc/2;
% ranges = newmaxf-newminf;
% incs = ranges/100;
% thresholds = newminf:incs:newmaxf;
% xaxiss = zeros(length(thresholds),1);
% yaxiss = zeros(length(thresholds),1);
% sensitivityss = zeros(length(thresholds),1);
% specificityss = zeros(length(thresholds),1);
% 
% if mean(posf)>=mean(negf)
% 
%     for num = 1:length(thresholds)
% 
%         thr = thresholds(num);
%         y = zeros(length(feature),1);
%         y(feature>=thr) = 1;
%         y(feature<thr) = -1;
%         B = y - label;
%         tp = length(find(B==0&y==1));
%         tn = length(find(B==0&y==-1));
%         fp = length(find(B==2));
%         fn = length(find(B==-2));
%         sensitivityss(num) = tp/(tp+fn+10e-200);
%         specificityss(num) = tn/(tn+fp+10e-200);
%         xaxiss(num) = 1-specificityss(num);
%         yaxiss(num) = sensitivityss(num);
% 
%     end
% 
% else
%     for num = 1:length(thresholds)
% 
%         thr = thresholds(num);
%         y = zeros(length(feature),1);
%         y(feature<thr) = 1;
%         y(feature>=thr) = -1;
%         B = y - label;
%         tp = length(find(B==0&y==1));
%         tn = length(find(B==0&y==-1));
%         fp = length(find(B==2));
%         fn = length(find(B==-2));
%         sensitivityss(num) = tp/(tp+fn+10e-200);
%         specificityss(num) = tn/(tn+fp+10e-200);
%         xaxiss(num) = 1-specificityss(num);
%         yaxiss(num) = sensitivityss(num);
% 
%     end
% 
% end
% 
% cutss = sensitivityss+specificityss;
% [maxtemp2,ind] = max(cutss);
% 
% if ~isempty(cutss)
%    cutofftemp2= thresholds(ind); 
% else
%     cutofftemp2 = min(feature);
% end
% 
% 
% 
% cutoff = cutofftemp2;
% xaxisf = [xaxis;xaxiss];
% yaxisf = [yaxis;yaxiss];
% 
% [xaxisf,inds] =sort(xaxisf,'descend');
% yaxisf = yaxisf(inds);
% auc= abs(trapz(xaxisf,yaxisf));