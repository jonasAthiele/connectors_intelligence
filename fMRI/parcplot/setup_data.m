
%%

addpath(strcat(pwd,'/parcplot/src/'))
addpath(strcat(pwd,'/parcplot/data/'))
addpath(genpath(strcat(pwd,'/parcplot/src/external/')))

%% load surface info and save as mat

mkdir([pwd '/parcplot/data/fsaverage/mat/'])

%for fff = {'sphere','smoothwm','inflated_pre','inflated'} 
for fff = {'inflated'} 
   
    surfStruct = load_surfStruct([pwd '/parcplot/data/'],'fsaverage',fff{1}) ;
    fileName = [pwd '/parcplot/data/fsaverage/mat/fsaverage_',fff{1},'.mat' ] ;
    save(fileName,'surfStruct','-v7.3') ;
end

%% load annotations
    
% initialize a map
allAnnots = containers.Map ;


%% schaefer-yeo17

for iii = [ 200 ] 

    currName = ['schaefer' num2str(iii) '-yeo17']
    tmpAnnot = load_annotStruct([pwd '/parcplot/data/'],'fsaverage',currName) ;
     
    tmpAnnot.combo_table = [ tmpAnnot.LH.ct.table(2:end,:) ; tmpAnnot.RH.ct.table(2:end,:) ] ;
    tmpAnnot.roi_ids = [ tmpAnnot.LH.ct.table(2:end,5) ; tmpAnnot.RH.ct.table(2:end,5) ] ;
    tmpAnnot.combo_names = [ tmpAnnot.LH.ct.struct_names(2:end) ; tmpAnnot.RH.ct.struct_names(2:end) ] ;
    
    tmpAnnot.LH.border = get_parc_borders(...
        tmpAnnot.LH.labs,surfStruct.LH.nbrs,0) ;
    tmpAnnot.RH.border = get_parc_borders(...
        tmpAnnot.RH.labs,surfStruct.RH.nbrs,0) ;

    allAnnots(currName) = tmpAnnot ;

end



%% save it

fileName = [pwd '/parcplot/data/fsaverage/mat/fsaverage_annots.mat' ] ;
save(fileName,'allAnnots','-v7.3') ;
ls(fileName)
