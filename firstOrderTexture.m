%% Radiogenomics calculation engine -1 
% This radiogenomic engine calculates 1st order image texture
% Imageprocessing toolbox and Gray Level Run Length Matrix Toolbox are
% required. Gray Level Run Length Matrix Toolbox can be downloaded 
% from https://jp.mathworks.com/matlabcentral/fileexchange/17482-gray-level-run-length-matrix-toolbox


function finalOutPut = firstOrderTexture(path_image, path_voi, imageType, voiType)

cd(path_image)
imageDirListing = dir('*.nii');

cd(path_voi)
voiDirListing = dir('*.nii');

for dI = 1:length(imageDirListing)
    list1 = imageDirListing(dI).name;
    strlist1num = regexp(list1,'(\d{5})','tokens');
    list_one(dI) = str2double(strlist1num{1});
end


for dII = 1:length(voiDirListing)
    list2= voiDirListing(dII).name;
    strlist3num = regexp(list2,'(\d{5})','tokens');
    list_two(dII) = str2double(strlist3num{1});
end

listMatch = intersect(list_one,list_two);
listtemp1 = imageDirListing(1).name;
listtemp2 = voiDirListing(1).name;

clear imageDirListing;
clear voiDirListing;


h = waitbar(0,'Now calculating 1st order image texture');

for d = 1:length(listMatch)
    strList = sprintf('%05d',listMatch(d));
    imagelist (d,1) = listMatch(d);
    imageFileName = regexprep(listtemp1, '(\d{5})', strList);
    voiFileName = regexprep(listtemp2, '(\d{5})', strList);
    
    imageFileNameFull = fullfile(path_image,imageFileName);
    voiFileNameFull = fullfile(path_voi,voiFileName);
    
    imageNii = load_nii(imageFileNameFull);
    voiNii = load_nii(voiFileNameFull);
    
    [mx, my, mz]=size(imageNii.img);
    
     
    pcolor_max = max(max(max(imageNii.img)));
    pcolor_min = min(min(min(imageNii.img)));
    voiNii.img = double(voiNii.img);
    voiExt = (voiNii.img~=0);
    
    valueofvoi = double(double(imageNii.img).*double(voiExt));
    finalvalue = valueofvoi(valueofvoi~=0);
    
if sum(finalvalue)~=0
        %%%%%% 1st order texture values %%%%%%%%
        aveg(d,1) = mean(finalvalue);
        stdv(d,1) = std(finalvalue);
        maxx(d,1) = max(finalvalue);
        minn(d,1) = min(finalvalue);
        med(d,1)  = median(finalvalue);
        mod(d,1)  = mode(finalvalue);
        vari(d,1)  = var(finalvalue);
        rms(d,1)  = sqrt(1/length(finalvalue).*(sum(finalvalue.^2)));
        

        %%%%%% Entropy calculation %%%%%%%%
        % set range for histogram
        range = zeros(1,256);
        for n=1:256
            range(n)=pcolor_min+((pcolor_max-pcolor_min)/255*(n-1));
        end
        % calculate histogram
        [p, ~] = histcounts(finalvalue, range);
        % remove zero entries in p
        p(p==0) = [];
        % normalize p so that sum(p) is one.
        p = p ./ numel(finalvalue);
        % calculate entropy(entr)
        entr(d,1) = -sum(p.*log2(p));
        
        %%%%%% Kurtosis calculation %%%%%%%%
        krt(d,1) = (sum((finalvalue-mean(finalvalue)).^4)./length(finalvalue)) ./ (var(finalvalue,1).^2);
        
        %%%%%% Skewns calculation %%%%%%%%
        skw(d,1) = (sum((finalvalue-mean(finalvalue)).^3)./length(finalvalue)) ./ (var(finalvalue,1).^1.5);

 
        
        
        
      
        %%%%%% Axial GLCM calculation %%%%%%%%
        
        [r,c,v] = ind2sub(size(voiExt),find(voiExt == 1));
        vxmin = min(r);
        vxmax = max(r);
        vymin = min(c);
        vymax = max(c);
        vzmin = min(v);
        vzmax = max(v);
        
        cutout4GLCM = valueofvoi(vxmin:vxmax, vymin:vymax, vzmin:vzmax);
        mask4GLRLM = double(voiExt(vxmin:vxmax, vymin:vymax, vzmin:vzmax));
        
    
    if (my < mx) && (my < mz)
        cutout4GLCM = double(permute(cutout4GLCM,[3 1 2]));
        mask4GLRLM = double(permute(mask4GLRLM,[3 1 2]));
    end
    
    if (mx < my) && (mx < mz)
        cutout4GLCM = double(permute(cutout4GLCM,[2 3 1]));
        mask4GLRLM = double(permute(mask4GLRLM,[2 3 1]));
    end
        
   
    
        [~,~,z_dim] = size(cutout4GLCM);
        
%         clear arraysum
%         
%         for z=1:z_dim
%             arraysum(z) = sum(sum(valueofvoi(:,:,z)));
%         end
%         newz = find(arraysum);
%         valueofvoi = valueofvoi(:,:,newz);
%         valueofvoi(valueofvoi == 0) = NaN;
        
%         z_dim = size(valueofvoi,3);
        for z=1:z_dim
            offsets = [ 0 1; 0 2; 0 3;...
                -1 1; -2 2; -3 3;...
                -1 0; -2 0; -3 0;...
                -1 -1; -2 -2; -3 -3];
            glcm_mat = graycomatrix(cutout4GLCM(:,:,z),'GrayLimits',[minn(d,1) maxx(d,1)], 'NumLevels',16,'Offset',offsets);
            a = graycoprops(glcm_mat,'Contrast');
            contrast(z,:) = getfield(a, 'Contrast');
            b = graycoprops(glcm_mat,'Correlation');
            correlation(z,:) = getfield(b, 'Correlation');
            c = graycoprops(glcm_mat,'Energy');
            energy(z,:) = getfield(c, 'Energy');
            dd = graycoprops(glcm_mat,'Homogeneity');
            homogeniety(z,:) = getfield(dd, 'Homogeneity');
            
            Ig = cutout4GLCM(:,:,z);
%             mask = ones(size(Ig(:,:,1)));
            mask = mask4GLRLM(:,:,z);
            quantize = 16;
            
            if sum(sum(mask4GLRLM(:,:,z)))==0
                
            SRE(z) = NaN;
            LRE(z) = NaN;
            GLN(z) = NaN;
            RP(z) = NaN;
            RLN(z) = NaN;
            LGRE(z) = NaN;
            HGRE(z) = NaN;
            
            else
            
            [a1,a2,a3,a4,a5,a6,a7]  = glrlm(Ig,quantize,mask);
            SRE(z) = a1;
            LRE(z) = a2;
            GLN(z) = a3;
            RP(z) = a4;
            RLN(z) = a5;
            LGRE(z) = a6;
            HGRE(z) = a7;
            
            end
            
        end
        

        
%         meancontrast = NaN([length(listMatch), 12])
        
        for i=1:12
            meancontrast(d,i) = mean(contrast(:,i));
            meancorrelation(d,i) = mean(correlation(:,i));
            meanenergy(d,i) = mean(energy(:,i));
            meanhomogeniety(d,i) = mean(homogeniety(:,i));
        end
        
        GLCMcontrast_1(d,1) = mean(meancontrast(d,1:3:10),'omitnan');
        GLCMcontrast_2(d,1) = mean(meancontrast(d,2:3:11),'omitnan');
        GLCMcontrast_3(d,1) = mean(meancontrast(d,3:3:12),'omitnan');

        GLCMcorrelation_1(d,1) = mean(meancorrelation(d,1:3:10),'omitnan');
        GLCMcorrelation_2(d,1) = mean(meancorrelation(d,2:3:11),'omitnan');
        GLCMcorrelation_3(d,1) = mean(meancorrelation(d,3:3:12),'omitnan');
        
        GLCMenergy_1(d,1) = mean(meanenergy(d,1:3:10),'omitnan');
        GLCMenergy_2(d,1) = mean(meanenergy(d,2:3:11),'omitnan');
        GLCMenergy_3(d,1) = mean(meanenergy(d,3:3:12),'omitnan');
        
        GLCMhomogeniety_1(d,1) = mean(meanhomogeniety(d,1:3:10),'omitnan');
        GLCMhomogeniety_2(d,1) = mean(meanhomogeniety(d,2:3:11),'omitnan');
        GLCMhomogeniety_3(d,1) = mean(meanhomogeniety(d,3:3:12),'omitnan');
        
        GLRLMSre(d,1) = mean(SRE(:),'omitnan');
        GLRLMLre(d,1) = mean(LRE(:),'omitnan');
        GLRLMGln(d,1) = mean(GLN(:),'omitnan');
        GLRLMRp(d,1) = mean(RP(:),'omitnan');
        GLRLMRln(d,1) = mean(RLN(:),'omitnan');
        GLRLMLrge(d,1) = mean(LGRE(:),'omitnan');
        GLRLMHrge(d,1) = mean(HGRE(:),'omitnan');
        
        GLCMcontrast_1_SD(d,1) = std(meancontrast(d,1:3:10),'omitnan');
        GLCMcontrast_2_SD(d,1) = std(meancontrast(d,2:3:11),'omitnan');
        GLCMcontrast_3_SD(d,1) = std(meancontrast(d,3:3:12),'omitnan');

        GLCMcorrelation_1_SD(d,1) = std(meancorrelation(d,1:3:10),'omitnan');
        GLCMcorrelation_2_SD(d,1) = std(meancorrelation(d,2:3:11),'omitnan');
        GLCMcorrelation_3_SD(d,1) = std(meancorrelation(d,3:3:12),'omitnan');
        
        GLCMenergy_1_SD(d,1) = std(meanenergy(d,1:3:10),'omitnan');
        GLCMenergy_2_SD(d,1) = std(meanenergy(d,2:3:11),'omitnan');
        GLCMenergy_3_SD(d,1) = std(meanenergy(d,3:3:12),'omitnan');
        
        GLCMhomogeniety_1_SD(d,1) = std(meanhomogeniety(d,1:3:10),'omitnan');
        GLCMhomogeniety_2_SD(d,1) = std(meanhomogeniety(d,2:3:11),'omitnan');
        GLCMhomogeniety_3_SD(d,1) = std(meanhomogeniety(d,3:3:12),'omitnan');
        
        GLRLMSre_SD(d,1) = std(SRE(:),'omitnan');
        GLRLMLre_SD(d,1) = std(LRE(:),'omitnan');
        GLRLMGln_SD(d,1) = std(GLN(:),'omitnan');
        GLRLMRp_SD(d,1) = std(RP(:),'omitnan');
        GLRLMRln_SD(d,1) = std(RLN(:),'omitnan');
        GLRLMLrge_SD(d,1) = std(LGRE(:),'omitnan');
        GLRLMHrge_SD(d,1) = std(HGRE(:),'omitnan');

end
    waitbar(d/length(listMatch));
end
close(h);

indexName = {'Image_ID',[imageType,'_',voiType,'_Mean'],[imageType,'_',voiType,'_SD'],[imageType,'_',voiType,'_Var'],[imageType,'_',voiType,'_RMS'], ...
    [imageType,'_',voiType,'_Max'],[imageType,'_',voiType,'_Min'],[imageType,'_',voiType,'_Median'], ...
    [imageType,'_',voiType,'_Mode'],[imageType,'_',voiType,'_Entropy'],[imageType,'_',voiType,'_Kurtosis'],[imageType,'_',voiType,'_Skewness'], ...
    [imageType,'_',voiType,'_GLCMcontrast_1'], [imageType,'_',voiType,'_GLCMcontrast_2'], [imageType,'_',voiType,'_GLCMcontrast_3'], ...
    [imageType,'_',voiType,'_GLCMcorrelation_1'], [imageType,'_',voiType,'_GLCMcorrelation_2'], [imageType,'_',voiType,'_GLCMcorrelation_3'], ...
    [imageType,'_',voiType,'_GLCMenergy_1'], [imageType,'_',voiType,'_GLCMenergy_2'], [imageType,'_',voiType,'_GLCMenergy_3'], ...
    [imageType,'_',voiType,'_GLCMhomogeniety_1'], [imageType,'_',voiType,'_GLCMhomogeniety_2'], [imageType,'_',voiType,'_GLCMhomogeniety_3']...
    [imageType,'_',voiType,'_GLRLMSre'], [imageType,'_',voiType,'_GLRLMLre'], [imageType,'_',voiType,'_GLRLMGln'],  [imageType,'_',voiType,'_GLRLMRp'], ...
    [imageType,'_',voiType,'_GLRLMRln'], [imageType,'_',voiType,'_GLRLMLrge'], [imageType,'_',voiType,'_GLRLMHrge'],...
    [imageType,'_',voiType,'_GLCMcontrast_1_SD'], [imageType,'_',voiType,'_GLCMcontrast_2_SD'], [imageType,'_',voiType,'_GLCMcontrast_3_SD'], ...
    [imageType,'_',voiType,'_GLCMcorrelation_1_SD'], [imageType,'_',voiType,'_GLCMcorrelation_2_SD'], [imageType,'_',voiType,'_GLCMcorrelation_3_SD'], ...
    [imageType,'_',voiType,'_GLCMenergy_1_SD'], [imageType,'_',voiType,'_GLCMenergy_2_SD'], [imageType,'_',voiType,'_GLCMenergy_3_SD'], ...
    [imageType,'_',voiType,'_GLCMhomogeniety_1_SD'], [imageType,'_',voiType,'_GLCMhomogeniety_2_SD'], [imageType,'_',voiType,'_GLCMhomogeniety_3_SD']...
    [imageType,'_',voiType,'_GLRLMSre_SD'], [imageType,'_',voiType,'_GLRLMLre_SD'], [imageType,'_',voiType,'_GLRLMGln_SD'],  [imageType,'_',voiType,'_GLRLMRp_SD'], ...
    [imageType,'_',voiType,'_GLRLMRln_SD'], [imageType,'_',voiType,'_GLRLMLrge_SD'], [imageType,'_',voiType,'_GLRLMHrge_SD'],};


T = table(imagelist, aveg, stdv,vari, rms, maxx, minn, med, mod, entr, krt, skw,GLCMcontrast_1, GLCMcontrast_2, GLCMcontrast_3, ...
    GLCMcorrelation_1, GLCMcorrelation_2, GLCMcorrelation_3, GLCMenergy_1, GLCMenergy_2, GLCMenergy_3, ...
    GLCMhomogeniety_1, GLCMhomogeniety_2, GLCMhomogeniety_3, GLRLMSre, GLRLMLre, GLRLMGln, GLRLMRp, GLRLMRln, GLRLMLrge, GLRLMHrge,...
    GLCMcontrast_1_SD, GLCMcontrast_2_SD, GLCMcontrast_3_SD, ...
    GLCMcorrelation_1_SD, GLCMcorrelation_2_SD, GLCMcorrelation_3_SD, GLCMenergy_1_SD, GLCMenergy_2_SD, GLCMenergy_3_SD, ...
    GLCMhomogeniety_1_SD, GLCMhomogeniety_2_SD, GLCMhomogeniety_3_SD, GLRLMSre_SD, GLRLMLre_SD, GLRLMGln_SD, GLRLMRp_SD, GLRLMRln_SD, GLRLMLrge_SD, GLRLMHrge_SD);

T.Properties.VariableNames{'imagelist'} = char(indexName(1));
T.Properties.VariableNames{'aveg'} = char(indexName(2));
T.Properties.VariableNames{'stdv'} = char(indexName(3));
T.Properties.VariableNames{'vari'} = char(indexName(4));
T.Properties.VariableNames{'rms'} = char(indexName(5));
T.Properties.VariableNames{'maxx'} = char(indexName(6));
T.Properties.VariableNames{'minn'} = char(indexName(7));
T.Properties.VariableNames{'med'} = char(indexName(8));
T.Properties.VariableNames{'mod'} = char(indexName(9));
T.Properties.VariableNames{'entr'} = char(indexName(10));
T.Properties.VariableNames{'krt'} = char(indexName(11));
T.Properties.VariableNames{'skw'} = char(indexName(12));
T.Properties.VariableNames{'GLCMcontrast_1'} = char(indexName(13));
T.Properties.VariableNames{'GLCMcontrast_2'} = char(indexName(14));
T.Properties.VariableNames{'GLCMcontrast_3'} = char(indexName(15));
T.Properties.VariableNames{'GLCMcorrelation_1'} = char(indexName(16));
T.Properties.VariableNames{'GLCMcorrelation_2'} = char(indexName(17));
T.Properties.VariableNames{'GLCMcorrelation_3'} = char(indexName(18));
T.Properties.VariableNames{'GLCMenergy_1'} = char(indexName(19));
T.Properties.VariableNames{'GLCMenergy_2'} = char(indexName(20));
T.Properties.VariableNames{'GLCMenergy_3'} = char(indexName(21));
T.Properties.VariableNames{'GLCMhomogeniety_1'} = char(indexName(22));
T.Properties.VariableNames{'GLCMhomogeniety_2'} = char(indexName(23));
T.Properties.VariableNames{'GLCMhomogeniety_3'} = char(indexName(24));
T.Properties.VariableNames{'GLRLMSre'} = char(indexName(25));
T.Properties.VariableNames{'GLRLMLre'} = char(indexName(26));
T.Properties.VariableNames{'GLRLMGln'} = char(indexName(27));
T.Properties.VariableNames{'GLRLMRp'} = char(indexName(28));
T.Properties.VariableNames{'GLRLMRln'} = char(indexName(29));
T.Properties.VariableNames{'GLRLMLrge'} = char(indexName(30));
T.Properties.VariableNames{'GLRLMHrge'} = char(indexName(31));
T.Properties.VariableNames{'GLCMcontrast_1_SD'} = char(indexName(32));
T.Properties.VariableNames{'GLCMcontrast_2_SD'} = char(indexName(33));
T.Properties.VariableNames{'GLCMcontrast_3_SD'} = char(indexName(34));
T.Properties.VariableNames{'GLCMcorrelation_1_SD'} = char(indexName(35));
T.Properties.VariableNames{'GLCMcorrelation_2_SD'} = char(indexName(36));
T.Properties.VariableNames{'GLCMcorrelation_3_SD'} = char(indexName(37));
T.Properties.VariableNames{'GLCMenergy_1_SD'} = char(indexName(38));
T.Properties.VariableNames{'GLCMenergy_2_SD'} = char(indexName(39));
T.Properties.VariableNames{'GLCMenergy_3_SD'} = char(indexName(40));
T.Properties.VariableNames{'GLCMhomogeniety_1_SD'} = char(indexName(41));
T.Properties.VariableNames{'GLCMhomogeniety_2_SD'} = char(indexName(42));
T.Properties.VariableNames{'GLCMhomogeniety_3_SD'} = char(indexName(43));
T.Properties.VariableNames{'GLRLMSre_SD'} = char(indexName(44));
T.Properties.VariableNames{'GLRLMLre_SD'} = char(indexName(45));
T.Properties.VariableNames{'GLRLMGln_SD'} = char(indexName(46));
T.Properties.VariableNames{'GLRLMRp_SD'} = char(indexName(47));
T.Properties.VariableNames{'GLRLMRln_SD'} = char(indexName(48));
T.Properties.VariableNames{'GLRLMLrge_SD'} = char(indexName(49));
T.Properties.VariableNames{'GLRLMHrge_SD'} = char(indexName(50));

finalOutPut = T;





