function [rms_whole, rms_in, rms_out, diffmap] = calc_def_rms(defA, defB)
%return rms between two deformation field images
if nargin < 1
    [ui_na, ui_pa, ~] = uigetfile('*.nii', 'Choose deformation Image');
    defA = fullfile(ui_pa, ui_na);
    [ui_nb, ui_pb, ~] = uigetfile('*.nii', 'Choose deformation Image');
    defB = fullfile(ui_pb, ui_nb);
end
disp(defA);
disp(defB);
spmDir = spm('Dir');
brainmask = fullfile(spmDir,'tpm','mask_ICV.nii');
lesmask = '';
les_pattern_id = '_M';
if strfind(defA, les_pattern_id)
    pth = fileparts(defA);
    lesmaskfile = dir(fullfile(pth, 'cLesion*.nii'));
    if ~isempty(lesmaskfile)
        lesmask = fullfile(pth,lesmaskfile(1).name);
    end
elseif strfind(defB, les_pattern_id)
    pth = fileparts(defB);
    lesmaskfile = dir(fullfile(pth, 'cLesion*.nii'));
    if ~isempty(lesmaskfile)
        lesmask = fullfile(pth,lesmaskfile(1).name);
    end
end

%brainmask = fullfile(ui_pa,'BrainMask.nii'); % for ANTS
niiA = nii_tool('load', defA);
niiB = nii_tool('load', defB);
niiM = nii_tool('load', brainmask);
imgM = double(niiM.img);
if ~isempty(lesmask)
    niiLesM = nii_tool('load', lesmask);
    imgLesM = double(niiLesM.img);
else
    imgLesM = double(ones(size(niiM.img)));
end


sA = size(niiA.img);
sB = size(niiB.img);
pdA = round(niiA.hdr.pixdim(2:4),1);
pdB = round(niiB.hdr.pixdim(2:4),1);

if ~isequal(sA, sB)
    error('Images must be same dimensions');
end

if ~isequal(pdA, pdB)
    error('pixle dimensions do not match between images');
end

pixdim = round(niiA.hdr.pixdim(2:4),1);
%%%%%%%%%%%%%%%%%%% whole brain (lesion + outside lesion)
imgAx = niiA.img(:,:,:,:,1).*imgM;
imgAy = niiA.img(:,:,:,:,2).*imgM;
imgAz = niiA.img(:,:,:,:,3).*imgM;

imgBx = niiB.img(:,:,:,:,1).*imgM;
imgBy = niiB.img(:,:,:,:,2).*imgM;
imgBz = niiB.img(:,:,:,:,3).*imgM;

dab = sqrt((imgAx-imgBx).^2 + (imgAy-imgBy).^2 + (imgAz-imgBz).^2);
%rms = sqrt(sum((dab(:).^2)/length(dab(:) > 0)))
dab_grtz = dab(dab>0);
rms_whole = sqrt (mean (dab_grtz.^2));
%rmsmm = rmsvox * pixdim(1) * pixdim(2) * pixdim(3);
diffmap = dab;  %* pixdim(1) * pixdim(2) * pixdim(3); % in mm
niidiff = niiA;
niidiff.hdr.fname = fullfile(pth,['diffmap_' lesmaskfile(1).name]);
niidiff.img = diffmap;
nii_tool('save', niidiff, niidiff.hdr.fname);
%rmsmap = (dab.^2)/length(dab(:) > 0);

%%%%%%%%%%%%%%%%% only inside lesion
imgAx = niiA.img(:,:,:,:,1).*imgM.*imgLesM;
imgAy = niiA.img(:,:,:,:,2).*imgM.*imgLesM;
imgAz = niiA.img(:,:,:,:,3).*imgM.*imgLesM;

imgBx = niiB.img(:,:,:,:,1).*imgM.*imgLesM;
imgBy = niiB.img(:,:,:,:,2).*imgM.*imgLesM;
imgBz = niiB.img(:,:,:,:,3).*imgM.*imgLesM;

dab = sqrt((imgAx-imgBx).^2 + (imgAy-imgBy).^2 + (imgAz-imgBz).^2);
dab_grtz = dab(dab>0);
rms_in = sqrt (mean (dab_grtz.^2));
%rms = sqrt(sum((dab(:).^2)/length(dab(:))))
%rmsmm = rmsvox * pixdim(1) * pixdim(2) * pixdim(3);
diffmap = dab;  %* pixdim(1) * pixdim(2) * pixdim(3); % in mm
niidiff = niiA;
niidiff.hdr.fname = fullfile(pth,['diffmap_inside_lesion_' lesmaskfile(1).name]);
niidiff.img = diffmap;
nii_tool('save', niidiff, niidiff.hdr.fname);
%rmsmap = (dab.^2)/length(dab(:) > 0);

%%%%%%%%%%%%%%% only outside lesion
if isempty(lesmask)
    imgLesM = double(zeros(size(niiM.img)));
end
imgAx = niiA.img(:,:,:,:,1).*(imgM-imgLesM);
imgAy = niiA.img(:,:,:,:,2).*(imgM-imgLesM);
imgAz = niiA.img(:,:,:,:,3).*(imgM-imgLesM);

imgBx = niiB.img(:,:,:,:,1).*(imgM-imgLesM);
imgBy = niiB.img(:,:,:,:,2).*(imgM-imgLesM);
imgBz = niiB.img(:,:,:,:,3).*(imgM-imgLesM);

dab = sqrt((imgAx-imgBx).^2 + (imgAy-imgBy).^2 + (imgAz-imgBz).^2);
%rms = sqrt(sum((dab(:).^2)/length(dab(:) > 0)))
dab_grtz = dab(dab>0);
rms_out = sqrt (mean (dab_grtz.^2));
%rmsmm = rmsvox * pixdim(1) * pixdim(2) * pixdim(3);
diffmap = dab;  %* pixdim(1) * pixdim(2) * pixdim(3); % in mm
niidiff = niiA;
niidiff.hdr.fname = fullfile(pth,['diffmap_outside_lesion_' lesmaskfile(1).name]);
niidiff.img = diffmap;
nii_tool('save', niidiff, niidiff.hdr.fname);
%rmsmap = (dab.^2)/length(dab(:) > 0);




