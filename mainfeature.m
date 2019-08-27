clc;
clear 
str1 = 'D:/Code/LungImmu/tumor_ctmask/';
filestruct=dir(fullfile(str1,'*.mha*'));

str2 ='D:/Data/LungImmunewCase_fold/Mat/';
str ='D:/Data/LungImmunewCase_fold/';

% normPETspace = [2.734375 2.734375 2];
% normPETspace = [3.645833254 3.645833254  3.27];
normPETspace = [1 1 1];
VoxelVolumePET = normPETspace(1)*normPETspace(2)*normPETspace(3)*0.001;

normCTspace = [1 1 1];
VoxelVolumeCT = normCTspace(1)*normCTspace(2)*normCTspace(3)*0.001;


D = 7;
numlevel = 2^D;
warning off
for num=1:length(filestruct)
     close all
     filename = filestruct(num).name;
     filename2 = filename(1:end-17)
     dataset = filename(1:6);
     
     pathname_pet = sprintf(strcat(str2,dataset,'/pet.mat'));
     pathname_ct = sprintf(strcat(str2,dataset,'/ct.mat'));
     load(pathname_pet);
     load(pathname_ct);
          
   str1s = 'D:/Code/LungImmu/tumor_ctmask/';
   str3 = '_lungnodulect.mat';
   pathname = sprintf(strcat(str1s,filename2,str3));
   load(pathname);
   nRange = range;
   str1s = 'D:/Code/LungImmu/tumor_mask/';
   str3 = '_lungnodule.mat';
   pathname = sprintf(strcat(str1s,filename2,str3));
   load(pathname);
   sRange = srange;
     
    sourcedir0 = sprintf(strcat(str,dataset,'/ST0/SE0/'));
    filestruct0=dir(fullfile(sourcedir0,'*.*'));
    filenum0 = size(filestruct0,1);

    sourcedir1 = sprintf(strcat(str,dataset,'/ST0/SE1/'));
    filestruct1=dir(fullfile(sourcedir1,'*.*'));
    filenum1 = size(filestruct1,1);

    filename0=[sourcedir0,filestruct0(3).name];
    dataheadct = dicominfo(filename0);

    filename2=[sourcedir0,filestruct0(4).name];
    dataheadct2 = dicominfo(filename0);

    filename1=[sourcedir1,filestruct1(3).name];
    dataheadpet = dicominfo(filename1);    
        
   try
         middleoff = (dataheadpet.PixelSpacing./dataheadct.PixelSpacing).*double([dataheadpet.Width;dataheadpet.Height]);
         middleoff = floor((middleoff - double([dataheadct.Width;dataheadct.Height]))/2);
         offct = [(dataheadct.ImagePositionPatient(2)-dataheadpet.ImagePositionPatient(2))/dataheadct.PixelSpacing(2); (dataheadct.ImagePositionPatient(1)-dataheadpet.ImagePositionPatient(1))/dataheadct.PixelSpacing(1)]-middleoff;
         poffct = round(abs(offct))

    catch
         middleoff = (dataheadpet.PixelSpacing./dataheadct2.PixelSpacing).*double([dataheadpet.Width;dataheadpet.Height]);
         middleoff = floor((middleoff - double([dataheadct2.Width;dataheadct2.Height]))/2);
         offct = [(dataheadct2.ImagePositionPatient(2)-dataheadpet.ImagePositionPatient(2))/dataheadct2.PixelSpacing(2); (dataheadct2.ImagePositionPatient(1)-dataheadpet.ImagePositionPatient(1))/dataheadct2.PixelSpacing(1)]-middleoff;
         poffct = round(abs(offct))
    end
        
        

    dataheadpet = PETInfo;
    pet = suv;
    
    clear suv;
    clear rescalepet;
    clear PETInfo;
    range.srangex = sRange.rangex;
    range.srangey = sRange.rangey;
    range.srangez = sRange.rangez;
    range.lrangex = nRange.rangex;
    range.lrangey = nRange.rangey;
    range.lrangez = nRange.rangez;

    
    
    dataheadct = CTInfo;
    ct = rescalect;
    clear rescalect
   petmaskful = zeros(size(pet,1),size(pet,2),size(pet,3));
   petmaskful(range.srangex,range.srangey,range.srangez)=lungnodulemask;
   ctmaskful = zeros(size(ct,1),size(ct,2),size(ct,3));
   ctmaskful(range.lrangex,range.lrangey,range.lrangez)=lungnodulectmask;
       
   if poffct(1)|poffct(2)
%%
       ctnew = zeros(size(ct));
       ctmasknew = zeros(size(ctmaskful));
       if offct(1)<=0
             ctnew(1:size(ctnew,1)-poffct(1),:,:) = ct(poffct(1)+1:end,:,:);
             ctmasknew (1:size(ctnew,1)-poffct(1),:,:) = ctmaskful(poffct(1)+1:end,:,:);
           else
             ctnew(poffct(1)+1:end,:,:) = ct(1:size(ctnew,1)-poffct(1),:,:);  
             ctmasknew(poffct(1)+1:end,:,:) = ctmaskful(1:size(ctnew,1)-poffct(1),:,:);  
       end

       ct = ctnew;
       ctnew = zeros(size(ct));

       ctmaskful = ctmasknew;
       ctmasknew = zeros(size(ctmaskful));


       if offct(2)<=0
         ctnew(:,1:size(ctnew,2)-poffct(2),:) = ct(:,poffct(2)+1:end,:);
         ctmasknew(:,1:size(ctnew,2)-poffct(2),:) = ctmaskful(:,poffct(2)+1:end,:);
       else
         ctnew(:,poffct(2)+1:end,:) = ct(:,1:size(ctnew,2)-poffct(2),:); 
         ctmasknew(:,1:size(ctnew,2)-poffct(2),:) = ctmaskful(:,poffct(2)+1:end,:);
       end

       ct = ctnew;
       ctmaskful = ctmasknew; 
   end
    


   range.srangez = min(range.srangez(1),range.lrangez(1)):max(range.srangez(end),range.lrangez(end));
   range.lrangez = range.srangez;
       
   pet = pet(:,:,range.srangez);
   ct = ct(:,:,range.lrangez);  
   petmaskful = petmaskful(:,:,range.srangez);
   ctmaskful = ctmaskful(:,:,range.lrangez);  

       
   pet = resize_m(pet,[dataheadpet.spacex dataheadpet.spacey dataheadpet.spacez],normPETspace);
   ct = resize_m(ct,[dataheadct.spacex dataheadct.spacey dataheadct.spacez],normCTspace);       
   petmaskful = resize_m(petmaskful,[dataheadpet.spacex dataheadpet.spacey dataheadpet.spacez],normPETspace);
   temppetnodule = pet.*petmaskful;
       
   petmaskful(pet<min(0.25*max(temppetnodule(:)),2.5))=0;

%        petmaskful(petmaskful>0.8)=1;     
%        petmaskful(petmaskful<=0.8)=0;   
   ctmaskful = resize_m(ctmaskful,[dataheadct.spacex dataheadct.spacey dataheadct.spacez],normCTspace);
   ctmaskful(ctmaskful>0.6)=1;
   ctmaskful(ctmaskful<=0.6)=0;
%        petmaskful_ctr = resize_m(petmaskful,normPETspace, normCTspace);
%        pet_ctr = resize_m(pet,normPETspace, normCTspace);
   petmaskful_ctr = magnify3(petmaskful, normPETspace,normCTspace,size(petmaskful),size(ctmaskful));
   pet_ctr = magnify3(pet, normPETspace,normCTspace,size(pet),size(ct));
       
   [x,y,z] = ind2sub(size(petmaskful),find(petmaskful));
   srangex = max((min(x)-1),1):min((max(x)+1),size(petmaskful,1));
   srangey = max((min(y)-1),1):min((max(y)+1),size(petmaskful,2));
   srangez = max((min(z)-1),1):min((max(z)+1),size(petmaskful,3));
%        srangex = max((min(x)-(max(x)-min(x))-1),1):min((max(x)+(max(x)-min(x))+1),size(petmaskful,1));
%        srangey = max((min(y)-(max(y)-min(y))-1),1):min((max(y)+(max(y)-min(y))+1),size(petmaskful,2));
%        srangez = max((min(z)-(max(z)-min(z))-1),1):min((max(z)+(max(z)-min(z))+1),size(petmaskful,3));

   range2.srangex = srangex;
   range2.srangey = srangey;
   range2.srangez = srangez;

   [x,y,z] = ind2sub(size(ctmaskful),find(ctmaskful));
   lrangex = max((min(x)-1),1):min((max(x)+1),size(ctmaskful,1));
   lrangey = max((min(y)-1),1):min((max(y)+1),size(ctmaskful,2));
   lrangez = max((min(z)-1),1):min((max(z)+1),size(ctmaskful,3));
%        
%        lrangex = max((min(x)-(max(x)-min(x))-1),1):min((max(x)+(max(x)-min(x))+1),size(ctmaskful,1));
%        lrangey = max((min(y)-(max(y)-min(y))-1),1):min((max(y)+(max(y)-min(y))+1),size(ctmaskful,2));
%        lrangez = max((min(z)-(max(z)-min(z))-1),1):min((max(z)+(max(z)-min(z))+1),size(ctmaskful,3));

   range2.lrangex = lrangex;
   range2.lrangey = lrangey;
   range2.lrangez = lrangez;


   partpet = pet(range2.srangex,range2.srangey,range2.srangez);
   partct = ct(range2.lrangex,range2.lrangey,range2.lrangez);
   petmask = petmaskful(range2.srangex,range2.srangey,range2.srangez);
   ctmask = ctmaskful(range2.lrangex,range2.lrangey,range2.lrangez);
   partpet_ctr = pet_ctr(range2.lrangex,range2.lrangey,range2.lrangez);
   petmask_ctr = petmaskful_ctr(range2.lrangex,range2.lrangey,range2.lrangez);

           
       sigma = 0.5;%10
       filsize = 6*sigma;
       partpet = Gaussian3D(filsize,sigma,partpet);
       partct = Gaussian3D(filsize,sigma,partct);
       partfuse = normalization(partpet_ctr) +0.6*normalization(partct);
       petnodule = partpet.*petmask;
       petvolume = length(find(petmask==1))*VoxelVolumePET;
       ctnodule = partct.*ctmask;  
       ctvolume = length(find(ctmask==1))*VoxelVolumeCT;
       petnodule_ctr = partpet_ctr.*petmask_ctr;
       fusemask_ctr = petmask_ctr+ctmask;
       fusemask_ctr(fusemask_ctr~=0)=1;
       fusenodule = partfuse.*fusemask_ctr;
       fusevolume = length(find(fusemask_ctr==1))*VoxelVolumeCT;
       
       petctmask = petmask_ctr;%.*ctmask;
       petctmask = Refinesegment3(petctmask);
       petctnodule = petctmask.*partpet_ctr;       
       petctnodule(petctmask==0) = nan;
       nSep = 2; tolX = 0.001;
       [IDX, threh, sep] = otsu(petctnodule,nSep);
       petmask_ctr_part1 = zeros(size(IDX));
       petmask_ctr_part1(IDX==1)=1;
       petmask_ctr_part2 = zeros(size(IDX));
       petmask_ctr_part2(IDX==2)=1;
       petnodule_ctr_part1 = petmask_ctr_part1.*petctnodule;
       petnodule_ctr_part2 = petmask_ctr_part2.*petctnodule;
       meanvpart1 = mean(petnodule_ctr_part1(petmask_ctr_part1~=0));
       meanvpart2 = mean(petnodule_ctr_part2(petmask_ctr_part2~=0));
       if meanvpart1>meanvpart2
           petmask_ctr_highpart = petmask_ctr_part1;
       else
           petmask_ctr_highpart = petmask_ctr_part2;
       end
       for f = 1:nSep
            Img = petmask_ctr_highpart;
            [L,numsep] = bwlabeln(Img,6);
            count = length(L(L==1));
            index = 1;
            for i = 2:numsep
                if length(L(L==i)) > count
                    count = length(L(L==i));
                    index = i;
                end
            end
            Img(find(L~=index)) = 0;
            petmask_ctr_highpart = Img;
       end
       
       petmask_ctr_highparttemp = Refinesegment_hl(petmask_ctr_highpart);
       if sum(petmask_ctr_highparttemp(:))
           petmask_ctr_highpart = petmask_ctr_highparttemp;
       end
       petmask_ctr_lowpart = petctmask;
       petmask_ctr_lowpart(petmask_ctr_highpart~=0)=0;
       petmask_ctr_lowparttemp = Refinesegment_hl(petmask_ctr_lowpart);
       if sum(petmask_ctr_lowparttemp(:))
           petmask_ctr_lowpart = petmask_ctr_lowparttemp;
       end
       petmask_ctr_highpart = petmask_ctr_highpart.*ctmask;
       petmask_ctr_lowpart = petmask_ctr_lowpart.*ctmask;
       
       ctnodule_highpart = ctnodule.*petmask_ctr_highpart;
       ctnodule_lowpart = ctnodule.*petmask_ctr_lowpart;
       cthvolume = length(find(petmask_ctr_highpart==1))*VoxelVolumeCT;
       ctlvolume = length(find(petmask_ctr_lowpart==1))*VoxelVolumeCT;
       
       petnoduletemp = double(petnodule);
       petnoduletemp(petmask==0) = nan;
       nSep = 2; tolX = 0.001;
       [IDX, threh, sep] = otsu(petnoduletemp,nSep);
       petmask_part1 = zeros(size(IDX));
       petmask_part1(IDX==1)=1;
       petmask_part2 = zeros(size(IDX));
       petmask_part2(IDX==2)=1;
       petnodule_part1 = petmask_part1.*petnodule;
       petnodule_part2 = petmask_part2.*petnodule;
       meanvpart1 = mean(petnodule_part1(petmask_part1~=0));
       meanvpart2 = mean(petnodule_part2(petmask_part2~=0));
       if meanvpart1>meanvpart2
           petmask_highpart = petmask_part1;
       else
           petmask_highpart = petmask_part2;
       end
       for f = 1:nSep
            Img = petmask_highpart;
            [L,numsep] = bwlabeln(Img,6);
            count = length(L(L==1));
            index = 1;
            for i = 2:numsep
                if length(L(L==i)) > count
                    count = length(L(L==i));
                    index = i;
                end
            end
            Img(find(L~=index)) = 0;
            petmask_highpart = Img;
       end
       
       petmask_lowpart = petmask;
       petmask_lowpart(petmask_highpart~=0)=0;
       petnodule_highpart = petnodule.*petmask_highpart;
       petnodule_lowpart = petnodule.*petmask_lowpart;
       pethvolume = length(find(petmask_highpart==1))*VoxelVolumePET;
       petlvolume = length(find(petmask_lowpart==1))*VoxelVolumePET;
       
%    

         [PETL,PETS ]=LongShortDiameter(petmask);
         [CTL,CTS ]=LongShortDiameter(ctmask);
         [PETHL,PETHS ]=LongShortDiameter(petmask_highpart);
         [PETLL,PETLS ]=LongShortDiameter(petmask_lowpart);
         [CTHL,CTHS ]=LongShortDiameter(petmask_ctr_highpart);
         [CTLL,CTLS ]=LongShortDiameter(petmask_ctr_lowpart);
         [FuseL,FuseS ]=LongShortDiameter(fusemask_ctr);
          Diameter(num,:)  = [PETL PETS CTL CTS PETHL PETHS PETLL PETLS CTHL CTHS CTLL CTLS FuseL FuseS];
       
       disp('PET features...');
       try
            PET_Feature_Values = [CaculateFeature(petnodule,numlevel) petvolume];
       catch
            PET_Feature_Values = zeros(1,65);
                
       end
       disp('CT features...');
       try
            [Status1, Feature_names, CT_Feature_Values2] = CalculateFeatures(partct,ctmask,normCTspace);
       catch
           CT_Feature_Values2 = zeros(1,364);
       end
       if ~Status1
           CT_Feature_Values2 = zeros(1,364);
       end
       
       try
            CT_Feature_Values = [CaculateFeature(ctnodule,numlevel) ctvolume];
       catch
           
           CT_Feature_Values = zeros(1,65);
           
       end
       disp('PETH features...');
       try
            PETH_Feature_Values = [CaculateFeature(petnodule_highpart,numlevel) pethvolume];     
       catch
            PETH_Feature_Values = zeros(1,65);    
           
       end
            
       disp('PETL features...');
       
       try
            PETL_Feature_Values = [CaculateFeature(petnodule_lowpart,numlevel) petlvolume];
       catch
            PETL_Feature_Values = zeros(1,65); 
           
       end
       disp('CTH features...');
       try
          CTH_Feature_Values = [CaculateFeature(ctnodule_highpart,numlevel) cthvolume];
       catch
          CTH_Feature_Values = zeros(1,65);
       end
       disp('CTL features...');
       try
          CTL_Feature_Values = [CaculateFeature(ctnodule_lowpart,numlevel) ctlvolume];
       catch
          CTL_Feature_Values = zeros(1,65);
       end
       disp('FUSE features...');
       
       try
           FUSE_Feature_Values = [CaculateFeature(fusenodule,numlevel) fusevolume];
       catch
           FUSE_Feature_Values = zeros(1,65);
           
       end
       Status=1;
  
        %%Mahmoud features
%        
%        [Status1, Feature_names, PET_Feature_Values] = CalculateFeatures(partpet,petmask,normPETspace);
%        [Status2, Feature_names, CT_Feature_Values] = CalculateFeatures(partct,ctmask,normCTspace);
%        [Status3, Feature_names, PETH_Feature_Values] = CalculateFeatures(petnodule,petmask_highpart,normPETspace);
%        [Status4, Feature_names, PETL_Feature_Values] = CalculateFeatures(petnodule,petmask_lowpart,normPETspace);
%        if sum(petmask_ctr_highpart(:))
%           [Status5, Feature_names, CTH_Feature_Values] = CalculateFeatures(ctnodule,petmask_ctr_highpart,normCTspace);
%        else
%            CTH_Feature_Values = zeros(1,65);Status5=1;
%        end
%        if sum(petmask_ctr_lowpart(:))
%           [Status6, Feature_names, CTL_Feature_Values] = CalculateFeatures(ctnodule,petmask_ctr_lowpart,normCTspace);
%        else
%            CTH_Feature_Values = zeros(1,65);Status6=1;
%        end       
%        
%        [Status7, Feature_names, FUSE_Feature_Values] = CalculateFeatures(partfuse,fusemask_ctr,normCTspace);
%        Status =1;% Status1 && Status2 && Status3 && Status4 && Status5 && Status6 && Status7 ;
       if (Status==1)
           Feature_Values(num,:) = [CT_Feature_Values PET_Feature_Values  CT_Feature_Values2 PETH_Feature_Values PETL_Feature_Values CTH_Feature_Values CTL_Feature_Values FUSE_Feature_Values...
               PETL_Feature_Values.*CTH_Feature_Values PETH_Feature_Values.*CTH_Feature_Values PETH_Feature_Values.*CTL_Feature_Values PETL_Feature_Values.*CTL_Feature_Values...
               CT_Feature_Values.*CTH_Feature_Values CT_Feature_Values.*CTL_Feature_Values PET_Feature_Values.*PETH_Feature_Values PET_Feature_Values.*PETL_Feature_Values...
               FUSE_Feature_Values.*PETH_Feature_Values FUSE_Feature_Values.*PETL_Feature_Values FUSE_Feature_Values.*CTH_Feature_Values FUSE_Feature_Values.*CTL_Feature_Values];
       
       end

       
       if mod(num,20)==0
           str1t = 'D:/Code/LungImmu/Feature/';
           str2t = 'new2feature_128_1d_1re_diameter.mat';
           pathnamet = strcat(str1t,str2t);
           save(pathnamet,'Diameter');%Feature_Values
       end
%        if(Status==1)


end
% fclose(out);
       str1t = 'D:/Code/LungImmu/Feature/';
       str2t = 'new2feature_128_1d_1re_diameter.mat';
       pathnamet = strcat(str1t,str2t);
        save(pathnamet,'Diameter');%Feature_Values%,'PETfeature','CTfeature','PETHfeature','PETLfeature','CTHfeature','CTLfeature','FUSEfeature');

% warning on;
% clc;
clear;
load('Study_final.mat');
% % load('clinicalfeature2.mat');
% % % Study = [651189 651459 661609];
% % load('Noct.mat');
% % count=1;
% % clear Study2
% % for i = 1:length(Study)
% %     if ~ismember(i,Out)
% %         Study2(count,1) = Study(i);
% %         Effind(count,1) = i;
% %         count = count+1;
% %     end
% % end
% % Study = Study2;
% % clinicalfeature = clinicalfeature(Effind,:);
% %  ind2 = find(clinicalfeature(:,7)>0&clinicalfeature(:,7)<=4);
% %  
% % Study = Study(ind2);
% % 
% % ind3 = find(Study>661609);
% % Study = Study(ind3);
% % normPETspace = [2.734375 2.734375 2];
% 
% % normPETspace = [3.645833254 3.645833254  3.27];
% normPETspace = [1 1 1];
% VoxelVolumePET = normPETspace(1)*normPETspace(2)*normPETspace(3)*0.001;
% 
% % % % % % normPETspace = [5.46875 5.46875 3.27];
% % normCTspace = [0.976562 0.976562 3.27];
% % normCTspace = [0.976562 0.976562 2];
% normCTspace = [1 1 1];
% VoxelVolumeCT = normCTspace(1)*normCTspace(2)*normCTspace(3)*0.001;
% % % normCTspace = [1.3672 1.3672 3.27];
% % % normCTspace = [1.3672 1.3672 3.27];
% 
% 
D = [3 4 5 6 7 8];
D = 7;
% % fname = 'Features_Wei.csv';
% % out = fopen('Features_Wei.csv', 'w+');
% % flag=0;
% 
% 
% % 
warning off;
for dnum = 1:length(D)
    numlevel = 2.^D(dnum)

    for num =7:length(Study)
%         num = Studys(nums);
       dataset = num2str(Study(num))

        str1 = 'D:/Data/LungImTxn/';
        sourcedir0 = sprintf(strcat(str1,num2str(dataset),'/ST0/SE0/'));
        filestruct0=dir(fullfile(sourcedir0,'*.*'));
        filenum0 = size(filestruct0,1);

        sourcedir1 = sprintf(strcat(str1,num2str(dataset),'/ST0/SE1/'));
        filestruct1=dir(fullfile(sourcedir1,'*.*'));
        filenum1 = size(filestruct1,1);

        filename0=[sourcedir0,filestruct0(3).name];
        dataheadct = dicominfo(filename0);

        filename2=[sourcedir0,filestruct0(4).name];
        dataheadct2 = dicominfo(filename0);


        filename1=[sourcedir1,filestruct1(3).name];
        dataheadpet = dicominfo(filename1);

        try
             middleoff = (dataheadpet.PixelSpacing./dataheadct.PixelSpacing).*double([dataheadpet.Width;dataheadpet.Height]);
             middleoff = floor((middleoff - double([dataheadct.Width;dataheadct.Height]))/2);
             offct = [(dataheadct.ImagePositionPatient(2)-dataheadpet.ImagePositionPatient(2))/dataheadct.PixelSpacing(2); (dataheadct.ImagePositionPatient(1)-dataheadpet.ImagePositionPatient(1))/dataheadct.PixelSpacing(1)]-middleoff;
             poffct = round(abs(offct))

        catch
             middleoff = (dataheadpet.PixelSpacing./dataheadct2.PixelSpacing).*double([dataheadpet.Width;dataheadpet.Height]);
             middleoff = floor((middleoff - double([dataheadct2.Width;dataheadct2.Height]))/2);
             offct = [(dataheadct2.ImagePositionPatient(2)-dataheadpet.ImagePositionPatient(2))/dataheadct2.PixelSpacing(2); (dataheadct2.ImagePositionPatient(1)-dataheadpet.ImagePositionPatient(1))/dataheadct2.PixelSpacing(1)]-middleoff;
             poffct = round(abs(offct))
        end
    
       str1 = 'D:/Data/LungImTxn/Mat/';
       str2 = '/pet.mat';
       str3 = '/ct.mat';
       pathname_pet = sprintf(strcat(str1,dataset,str2));
       pathname_ct = sprintf(strcat(str1,dataset,str3));
       load(pathname_pet);
       load(pathname_ct);
       str1s = 'D:/Code/LungImmu/tumor_mask/';
       str2s = '_lungnodulepet_1.mat';
       pathname = sprintf(strcat(str1s,dataset,str2s));
       load(pathname);%lungnodulemask
       sRange = range;
       str1s = 'D:/Code/LungImmu/tumor_ctmask/';
       str2s = '_lungnodulect.mat';
       pathname = sprintf(strcat(str1s,dataset,str2s));
       load(pathname); %lungnodulectmask range.srange&lrange
       nRange = range;
       dataheadpet = PETInfo;
       pet = suv;
    
       clear suv;
       clear rescalepet;
       clear PETInfo;
       range.srangex = sRange.srangex;
       range.srangey = sRange.srangey;
       range.srangez = sRange.srangez;
       range.lrangex = nRange.lrangex;
       range.lrangey = nRange.lrangey;
       range.lrangez = nRange.lrangez;
       
       dataheadct = CTInfo;
       ct = rescalect;
       clear rescalect
       petmaskful = zeros(size(pet,1),size(pet,2),size(pet,3));
       petmaskful(range.srangex,range.srangey,range.srangez)=lungnodulemask;
       ctmaskful = zeros(size(ct,1),size(ct,2),size(ct,3));
       ctmaskful(range.lrangex,range.lrangey,range.lrangez)=lungnodulectmask;
       
       if poffct(1)|poffct(2)
   %%
           ctnew = zeros(size(ct));
           ctmasknew = zeros(size(ctmaskful));
           if offct(1)<=0
                 ctnew(1:size(ctnew,1)-poffct(1),:,:) = ct(poffct(1)+1:end,:,:);
                 ctmasknew (1:size(ctnew,1)-poffct(1),:,:) = ctmaskful(poffct(1)+1:end,:,:);
               else
                 ctnew(poffct(1)+1:end,:,:) = ct(1:size(ctnew,1)-poffct(1),:,:);  
                 ctmasknew(poffct(1)+1:end,:,:) = ctmaskful(1:size(ctnew,1)-poffct(1),:,:);  
           end

           ct = ctnew;
           ctnew = zeros(size(ct));

           ctmaskful = ctmasknew;
           ctmasknew = zeros(size(ctmaskful));


           if offct(2)<=0
             ctnew(:,1:size(ctnew,2)-poffct(2),:) = ct(:,poffct(2)+1:end,:);
             ctmasknew(:,1:size(ctnew,2)-poffct(2),:) = ctmaskful(:,poffct(2)+1:end,:);
           else
             ctnew(:,poffct(2)+1:end,:) = ct(:,1:size(ctnew,2)-poffct(2),:); 
             ctmasknew(:,1:size(ctnew,2)-poffct(2),:) = ctmaskful(:,poffct(2)+1:end,:);
           end

           ct = ctnew;
           ctmaskful = ctmasknew; 
       end
    
    
    
%        normPETspace = [dataheadpet.spacex  dataheadpet.spacey  dataheadpet.spacez ];
%        VoxelVolumePET = normPETspace(1)*normPETspace(2)*normPETspace(3)*0.001;
%        normCTspace = [dataheadct.spacex dataheadct.spacey dataheadct.spacez];
%        VoxelVolumeCT = normCTspace(1)*normCTspace(2)*normCTspace(3)*0.001;
       
       pet = pet(:,:,range.srangez);
       ct = ct(:,:,range.lrangez);  
       petmaskful = petmaskful(:,:,range.srangez);
       ctmaskful = ctmaskful(:,:,range.lrangez);  

       
       pet = resize_m(pet,[dataheadpet.spacex dataheadpet.spacey dataheadpet.spacez],normPETspace);
       ct = resize_m(ct,[dataheadct.spacex dataheadct.spacey dataheadct.spacez],normCTspace);       
       petmaskful = resize_m(petmaskful,[dataheadpet.spacex dataheadpet.spacey dataheadpet.spacez],normPETspace);
       temppetnodule = pet.*petmaskful;
       
        petmaskful(pet<min(0.25*max(temppetnodule(:)),2.5))=0;
%        petmaskful(petmaskful>0.8)=1;     
%        petmaskful(petmaskful<=0.8)=0;   
       ctmaskful = resize_m(ctmaskful,[dataheadct.spacex dataheadct.spacey dataheadct.spacez],normCTspace);
       ctmaskful(ctmaskful>0.6)=1;
       ctmaskful(ctmaskful<=0.6)=0;
%        petmaskful_ctr = resize_m(petmaskful,normPETspace, normCTspace);
%        pet_ctr = resize_m(pet,normPETspace, normCTspace);
       petmaskful_ctr = magnify3(petmaskful, normPETspace,normCTspace,size(petmaskful),size(ctmaskful));
       pet_ctr = magnify3(pet, normPETspace,normCTspace,size(pet),size(ct));
       
       [x,y,z] = ind2sub(size(petmaskful),find(petmaskful));
       srangex = max((min(x)-2),1):min((max(x)+2),size(petmaskful,1));
       srangey = max((min(y)-2),1):min((max(y)+2),size(petmaskful,2));
       srangez = max((min(z)-2),1):min((max(z)+2),size(petmaskful,3));
%        srangex = max((min(x)-(max(x)-min(x))-1),1):min((max(x)+(max(x)-min(x))+1),size(petmaskful,1));
%        srangey = max((min(y)-(max(y)-min(y))-1),1):min((max(y)+(max(y)-min(y))+1),size(petmaskful,2));
%        srangez = max((min(z)-(max(z)-min(z))-1),1):min((max(z)+(max(z)-min(z))+1),size(petmaskful,3));
       
       range2.srangex = srangex;
       range2.srangey = srangey;
       range2.srangez = srangez;
       
       [x,y,z] = ind2sub(size(ctmaskful),find(ctmaskful));
       lrangex = max((min(x)-2),1):min((max(x)+2),size(ctmaskful,1));
       lrangey = max((min(y)-2),1):min((max(y)+2),size(ctmaskful,2));
       lrangez = max((min(z)-2),1):min((max(z)+2),size(ctmaskful,3));
%        
%        lrangex = max((min(x)-(max(x)-min(x))-1),1):min((max(x)+(max(x)-min(x))+1),size(ctmaskful,1));
%        lrangey = max((min(y)-(max(y)-min(y))-1),1):min((max(y)+(max(y)-min(y))+1),size(ctmaskful,2));
%        lrangez = max((min(z)-(max(z)-min(z))-1),1):min((max(z)+(max(z)-min(z))+1),size(ctmaskful,3));
       
       range2.lrangex = lrangex;
       range2.lrangey = lrangey;
       range2.lrangez = lrangez;
       
                
       partpet = pet(range2.srangex,range2.srangey,range2.srangez);
       partct = ct(range2.lrangex,range2.lrangey,range2.lrangez);
       petmask = petmaskful(range2.srangex,range2.srangey,range2.srangez);
       ctmask = ctmaskful(range2.lrangex,range2.lrangey,range2.lrangez);
       partpet_ctr = pet_ctr(range2.lrangex,range2.lrangey,range2.lrangez);
       petmask_ctr = petmaskful_ctr(range2.lrangex,range2.lrangey,range2.lrangez);
       
%        slice_outline = bwperim(petmask,4);
%        maskctv = partpet;
%        maskctv(slice_outline) = 1.1 * max(maskctv(:));
%        nouse = SliceBrowser(maskctv,'maskpetv');
%        slice_outline = bwperim(ctmask,4);
%        maskctv = partct;
%        maskctv(slice_outline) = 1.1 * max(maskctv(:));
%        nouse = SliceBrowser(maskctv,'maskpetv');
% %      
%        if num==51||num==53||num==54||num==56
%            petmask = imdilate(petmask,strel('sphere',1));
%            petmask(partpet<2)=0;
%            petmask_ctr = imdilate(petmask_ctr,strel('sphere',2));
%            petmask_ctr(partpet_ctr<2)=0;
%        end
           
       sigma = 0.5;%10
       filsize = 6*sigma;
       partpet = Gaussian3D(filsize,sigma,partpet);
       partct = Gaussian3D(filsize,sigma,partct);
       partfuse = normalization(partpet_ctr) +0.6*normalization(partct);
       petnodule = partpet.*petmask;
       petvolume = length(find(petmask==1))*VoxelVolumePET;
       ctnodule = partct.*ctmask;  
       ctvolume = length(find(ctmask==1))*VoxelVolumeCT;
       petnodule_ctr = partpet_ctr.*petmask_ctr;
       fusemask_ctr = petmask_ctr+ctmask;
       fusemask_ctr(fusemask_ctr~=0)=1;
       fusenodule = partfuse.*fusemask_ctr;
       fusevolume = length(find(fusemask_ctr==1))*VoxelVolumeCT;
       
       petctmask = petmask_ctr;%.*ctmask;
       petctmask = Refinesegment3(petctmask);
       petctnodule = petctmask.*partpet_ctr;       
       petctnodule(petctmask==0) = nan;
       nSep = 2; tolX = 0.001;
       [IDX, threh, sep] = otsu(petctnodule,nSep);
       petmask_ctr_part1 = zeros(size(IDX));
       petmask_ctr_part1(IDX==1)=1;
       petmask_ctr_part2 = zeros(size(IDX));
       petmask_ctr_part2(IDX==2)=1;
       petnodule_ctr_part1 = petmask_ctr_part1.*petctnodule;
       petnodule_ctr_part2 = petmask_ctr_part2.*petctnodule;
       meanvpart1 = mean(petnodule_ctr_part1(petmask_ctr_part1~=0));
       meanvpart2 = mean(petnodule_ctr_part2(petmask_ctr_part2~=0));
       if meanvpart1>meanvpart2
           petmask_ctr_highpart = petmask_ctr_part1;
       else
           petmask_ctr_highpart = petmask_ctr_part2;
       end
       for f = 1:nSep
            Img = petmask_ctr_highpart;
            [L,numsep] = bwlabeln(Img,6);
            count = length(L(L==1));
            index = 1;
            for i = 2:numsep
                if length(L(L==i)) > count
                    count = length(L(L==i));
                    index = i;
                end
            end
            Img(find(L~=index)) = 0;
            petmask_ctr_highpart = Img;
       end
       
       petmask_ctr_highparttemp = Refinesegment_hl(petmask_ctr_highpart);
       if sum(petmask_ctr_highparttemp(:))
           petmask_ctr_highpart = petmask_ctr_highparttemp;
       end
       petmask_ctr_lowpart = petctmask;
       petmask_ctr_lowpart(petmask_ctr_highpart~=0)=0;
       petmask_ctr_lowparttemp = Refinesegment_hl(petmask_ctr_lowpart);
       if sum(petmask_ctr_lowparttemp(:))
           petmask_ctr_lowpart = petmask_ctr_lowparttemp;
       end
       petmask_ctr_highpart = petmask_ctr_highpart.*ctmask;
       petmask_ctr_lowpart = petmask_ctr_lowpart.*ctmask;
       
       ctnodule_highpart = ctnodule.*petmask_ctr_highpart;
       ctnodule_lowpart = ctnodule.*petmask_ctr_lowpart;
       cthvolume = length(find(petmask_ctr_highpart==1))*VoxelVolumeCT;
       ctlvolume = length(find(petmask_ctr_lowpart==1))*VoxelVolumeCT;
       
       petnoduletemp = petnodule;
       petnoduletemp(petmask==0) = nan;
       nSep = 2; tolX = 0.001;
       [IDX, threh, sep] = otsu(petnoduletemp,nSep);
       petmask_part1 = zeros(size(IDX));
       petmask_part1(IDX==1)=1;
       petmask_part2 = zeros(size(IDX));
       petmask_part2(IDX==2)=1;
       petnodule_part1 = petmask_part1.*petnodule;
       petnodule_part2 = petmask_part2.*petnodule;
       meanvpart1 = mean(petnodule_part1(petmask_part1~=0));
       meanvpart2 = mean(petnodule_part2(petmask_part2~=0));
       if meanvpart1>meanvpart2
           petmask_highpart = petmask_part1;
       else
           petmask_highpart = petmask_part2;
       end
       for f = 1:nSep
            Img = petmask_highpart;
            [L,numsep] = bwlabeln(Img,6);
            count = length(L(L==1));
            index = 1;
            for i = 2:numsep
                if length(L(L==i)) > count
                    count = length(L(L==i));
                    index = i;
                end
            end
            Img(find(L~=index)) = 0;
            petmask_highpart = Img;
       end
       
       petmask_lowpart = petmask;
       petmask_lowpart(petmask_highpart~=0)=0;
       petnodule_highpart = petnodule.*petmask_highpart;
       petnodule_lowpart = petnodule.*petmask_lowpart;
       pethvolume = length(find(petmask_highpart==1))*VoxelVolumePET;
       petlvolume = length(find(petmask_lowpart==1))*VoxelVolumePET;
       
%        slice_outline = bwperim(petmask,4);
%        maskctv = partpet;
%        maskctv(slice_outline) = 1.1 * max(maskctv(:));
%        nouse = SliceBrowser(maskctv,'maskpetv');
%        slice_outline = bwperim(ctmask,4);
%        maskctv = partct;
%        maskctv(slice_outline) = 1.1 * max(maskctv(:));
%        nouse = SliceBrowser(maskctv,'maskpetv');
%        slice_outline = bwperim(petmask_ctr_lowpart,4);
%        maskctv = partct;
%        maskctv(slice_outline) = 1.1 * max(maskctv(:));
%        nouse = SliceBrowser(maskctv,'maskpetv');

       
       
       
       %%Wei features
       [PETL,PETS ]=LongShortDiameter(petmask);
         [CTL,CTS ]=LongShortDiameter(ctmask);
         [PETHL,PETHS ]=LongShortDiameter(petmask_highpart);
         [PETLL,PETLS ]=LongShortDiameter(petmask_lowpart);
         [CTHL,CTHS ]=LongShortDiameter(petmask_ctr_highpart);
         [CTLL,CTLS ]=LongShortDiameter(petmask_ctr_lowpart);
         [FuseL,FuseS ]=LongShortDiameter(fusemask_ctr);
         
         Diameter(num,:)  = [PETL PETS CTL CTS PETHL PETHS PETLL PETLS CTHL CTHS CTLL CTLS FuseL FuseS];
         
%        
%        disp('PET features...');
%        try
%             PET_Feature_Values = [CaculateFeature(petnodule,numlevel) petvolume];
%        catch
%             PET_Feature_Values = zeros(1,65);
%                 
%        end
%        disp('CT features...');
%        try
%             [Status1, Feature_names, CT_Feature_Values2] = CalculateFeatures(partct,ctmask,normCTspace);
%        catch
%            CT_Feature_Values2 = zeros(1,364);
%        end
%        if ~Status1
%            CT_Feature_Values2 = zeros(1,364);
%        end
%        
%        try
%             CT_Feature_Values = [CaculateFeature(ctnodule,numlevel) ctvolume];
%        catch
%            
%            CT_Feature_Values = zeros(1,65);
%            
%        end
%        disp('PETH features...');
%        try
%             PETH_Feature_Values = [CaculateFeature(petnodule_highpart,numlevel) pethvolume];     
%        catch
%             PETH_Feature_Values = zeros(1,65);    
%            
%        end
%             
%        disp('PETL features...');
%        
%        try
%             PETL_Feature_Values = [CaculateFeature(petnodule_lowpart,numlevel) petlvolume];
%        catch
%             PETL_Feature_Values = zeros(1,65); 
%            
%        end
%        disp('CTH features...');
%        try
%           CTH_Feature_Values = [CaculateFeature(ctnodule_highpart,numlevel) cthvolume];
%        catch
%           CTH_Feature_Values = zeros(1,65);
%        end
%        disp('CTL features...');
%        try
%           CTL_Feature_Values = [CaculateFeature(ctnodule_lowpart,numlevel) ctlvolume];
%        catch
%           CTL_Feature_Values = zeros(1,65);
%        end
%        disp('FUSE features...');
%        
%        try
%            FUSE_Feature_Values = [CaculateFeature(fusenodule,numlevel) fusevolume];
%        catch
%            FUSE_Feature_Values = zeros(1,65);
%            
%        end
%        Status=1;
%   
        %%Mahmoud features
%        
%        [Status1, Feature_names, PET_Feature_Values] = CalculateFeatures(partpet,petmask,normPETspace);
%        [Status2, Feature_names, CT_Feature_Values] = CalculateFeatures(partct,ctmask,normCTspace);
%        [Status3, Feature_names, PETH_Feature_Values] = CalculateFeatures(petnodule,petmask_highpart,normPETspace);
%        [Status4, Feature_names, PETL_Feature_Values] = CalculateFeatures(petnodule,petmask_lowpart,normPETspace);
%        if sum(petmask_ctr_highpart(:))
%           [Status5, Feature_names, CTH_Feature_Values] = CalculateFeatures(ctnodule,petmask_ctr_highpart,normCTspace);
%        else
%            CTH_Feature_Values = zeros(1,65);Status5=1;
%        end
%        if sum(petmask_ctr_lowpart(:))
%           [Status6, Feature_names, CTL_Feature_Values] = CalculateFeatures(ctnodule,petmask_ctr_lowpart,normCTspace);
%        else
%            CTH_Feature_Values = zeros(1,65);Status6=1;
%        end       
%        
%        [Status7, Feature_names, FUSE_Feature_Values] = CalculateFeatures(partfuse,fusemask_ctr,normCTspace);
%        Status =1;% Status1 && Status2 && Status3 && Status4 && Status5 && Status6 && Status7 ;
%        if (Status==1)
%            Feature_Values(num,:) = [CT_Feature_Values PET_Feature_Values  CT_Feature_Values2 PETH_Feature_Values PETL_Feature_Values CTH_Feature_Values CTL_Feature_Values FUSE_Feature_Values...
%                PETL_Feature_Values.*CTH_Feature_Values PETH_Feature_Values.*CTH_Feature_Values PETH_Feature_Values.*CTL_Feature_Values PETL_Feature_Values.*CTL_Feature_Values...
%                CT_Feature_Values.*CTH_Feature_Values CT_Feature_Values.*CTL_Feature_Values PET_Feature_Values.*PETH_Feature_Values PET_Feature_Values.*PETL_Feature_Values...
%                FUSE_Feature_Values.*PETH_Feature_Values FUSE_Feature_Values.*PETL_Feature_Values FUSE_Feature_Values.*CTH_Feature_Values FUSE_Feature_Values.*CTL_Feature_Values];
%        
%        end

       
       if mod(num,20)==0
           str1 = 'D:/Code/LungImmu/Feature/';
           str2 = 'feature_128_1d_1re_diameter.mat';
           pathname = strcat(str1,str2);
%            save(pathname,'Feature_Values');Diameter
            save(pathname,'Diameter');
       end
%        if(Status==1)
%         %write the header if it is the first image/patient
% %         if(flag==0)
% %             fdispf(out,'Image Name,');
% %             fdispf(out,'Mask Name,');
% %             for count=1:19
% %                 for j=1:size(Feature_names,2)
% %                     fdispf(out,'%s,',char(Feature_names{1,j}));
% %                 end
% %             end
% %             fdispf(out,'\n');
% %             flag=1;
% %         end
% %         
%         % write the feature values
%         fprintf(out,strcat('Image_',num2str(Study(num)),','));
%         fprintf(out,strcat('Mask_',num2str(Study(num)),','));
%         
%         for j=1:size(Feature_Values,2)
%             fprintf(out,'%s,',num2str(Feature_Values(num,j)));
%         end
%         fprintf(out,'\n');
%         %------------------------------------------------------------------------------------
%         end
    
           
           
       
       
%        slice_outline = bwperim(fusemask_ctr,4);
%        maskctv = partpet_ctr+partct;
%        maskctv(slice_outline) = 1.1 * max(maskctv(:));
%        nouse = SliceBrowser(maskctv,'maskpetv');

    end
end
% fclose(out);
            str1 = 'D:/Code/LungImmu/Feature/';
           str2 = 'feature_128_1d_1re_diameter.mat';
           pathname = strcat(str1,str2);
%            save(pathname,'Feature_Values');Diameter
            save(pathname,'Diameter');
            %         str1 = 'J:/Wei/';
%          str2 = 'feature_128_1d_1re.mat';
%        pathname = strcat(str1,str2);
%        save(pathname,'Feature_Values');
% warning on;