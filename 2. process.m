%% Select .mat files to process.
% function process()
allpath=[];
rootpath='./';
allfolder=dir(fullfile(rootpath,'*clr4_MBC'));
for i=1:length(allfolder)
    folderfile=dir([rootpath,'/',allfolder(i).name,'/*.mat']);
%     folderfile=folderfile(3:end);
    folderpath=cellfun(@(x)fullfile(rootpath,allfolder(i).name,x),{folderfile.name},'UniformOutput',0);
    allpath=[allpath,folderpath];
end

allpath
%% Make 2D maps of pixel intensity and detect fluctuations.
% allfile=dir('*.mat');
% allpath={allfile.name};
parfor i=1:length(allpath)
    allpath{i}
    singleworker(allpath{i});
    
end
% end