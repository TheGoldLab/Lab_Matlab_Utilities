function pathname = goldLabDataSessionLocalPathnames(tag)
% function pathname = goldLabDataSessionLocalPathnames(tag)
%
% Resource file to define local pathnames used by goldLabDataSession
%   utilities

if strcmp(tag, 'baseDataDirectory')
    pathname = '/Users/jigold/Library/CloudStorage/Box-Box/GoldLab/Data/Physiology/';

elseif strcmp(tag, 'pythonSearchPath')
    pathname = '/Users/jigold/miniconda3/bin:/Users/jigold/miniconda3/condabin:';

elseif strcmp(tag, 'pyramidSearchPath')
    pathname = '/Users/jigold/GoldWorks/Mirror_jigold/Manuscripts/2022_dlPFC/mfiles/dataSession/pyramid/';
end

