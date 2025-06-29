clear all
close all

dataPartitionPath = '/data/'; %'D:/'
dataImagingPartitionPath = '/data_imaging/'; %'F:/'
% % Load data:
% loadData = 0;
% Overwrite Nifti:
overwriteNifti = 0; % 0: if nifti conversion exists, just reads the images and headers, 1: covnerts and overwrites again
%% ADD PATHS
addpath([dataPartitionPath 'UNSAM/Brain/dicm2nii/'])
unsamPetQuantificationToolsPath = '/data/UNSAM/PET/UNSAMPETQuantificationTools/';

%% PATHS AND FILENAMES
niftiExtension = '.nii.gz';
dcmHeadersFilename = 'dcmHeaders.mat';
dicomDataPath = [dataPartitionPath '/UNSAM/CovidProject/Estudio/PET/'];
%dicomDataPath = [dataPartitionPath '/UNSAM/CovidProject2/Imaging/PET/'];
preprocessedDataPath = [dataPartitionPath '/UNSAM/CovidProject/Estudio/PreprocessedPET/'];
preprocessedMriDataPath = [dataPartitionPath '/UNSAM/CovidProject/Estudio/PreprocessedMRI/'];
preprocessedDataPath = [dataImagingPartitionPath '/CovidProject/Estudio/PreprocessedPET/'];
preprocessedMriDataPath = [dataImagingPartitionPath '/CovidProject/Estudio/PreprocessedMRI/'];
%preprocessedDataPath = [dataImagingPartitionPath '/CovidProject/Estudio2/PreprocessedPET/'];
%preprocessedMriDataPath = [dataImagingPartitionPath '/CovidProject/Estudio2/PreprocessedMRI/'];
niftiDataPath = [preprocessedDataPath '/Nifti/'];
quantPath = [preprocessedDataPath '/Quantification/'];
filenameUPQTScript = [preprocessedDataPath 'runUnsamPetQuantTools.sh'];
if ~isdir(niftiDataPath)
    mkdir(niftiDataPath)
end
if ~isdir(quantPath)
    mkdir(quantPath)
end
%% NAME SEQUENCES
dfltNamePetCtAc = {'x_DetailBR_CTAC_Brain', 'x_BR_CTAC_Brain'};
dfltNamePetCtAc = {'x_DetailBR_CTAC_Brain', 'x_BR_CTAC_Brain'};
tagNamePetCtAc = 'CTAC'; % If the default name is not found, use this tag to identify the image.
dfltNameT1 = 't1_mprage_1x1x1_estric_';
nameFmri = 'funcional_3s_30te_2_4_iso';
nameFieldmappingMag1 = 'gre_field_mapping_2mm_e1';
nameFieldmappingMag2 = 'gre_field_mapping_2mm_e2';
nameFieldmappingPhase = 'gre_field_mapping_2mm_phase';
%% CASES TO PROCESS
casesToProcess = [];%{'CP0002', 'CP0006', 'CP0007', 'CP0008', 'CP0009', 'CP0010'};
indicesToRunUPQT = [];
%% PROCESS EACH CASE
if isempty(casesToProcess)
    dirPaths = dir(dicomDataPath);
    casesToProcess = {dirPaths(3:end).name};
end

for i = 1 : numel(casesToProcess)
    %% CONVERT DICOM TO NIFTI
    niftiPathThisSubject = [niftiDataPath casesToProcess{i} '/'];
    dicomPathThisSubject = [dicomDataPath casesToProcess{i}];
    if ~exist(niftiPathThisSubject) || overwriteNifti
        converterOut = dicm2nii(dicomPathThisSubject, niftiPathThisSubject, niftiExtension);
    end    
    dcmTags{i} = load([niftiPathThisSubject dcmHeadersFilename]);
%     %% CREATE FOLDERS FOR MELODIC
%     melodicPath = [pathThisSubject '/Nifti/'];
    % Images for this subject:
    imageNamesPerSubject{i} = fieldnames(dcmTags{i}.h);
    %% GET THE PET CT ATTENUATEC CORRECTED IMAGE
    for j = 1 : numel(dfltNamePetCtAc)
        indexPetCtAcImage = find(strcmp(imageNamesPerSubject{i}, dfltNamePetCtAc{j}));
        if ~isempty(indexPetCtAcImage)
            break;
        end
    end
    if ~isempty(indexPetCtAcImage)
        % If there are more than 1 reconstructed image, use the last one.
        if numel(indexPetCtAcImage) > 1 % There are two PET images, use the second one.
            indexPetCtAcImage = indexPetCtAcImage(1);
        end
        namePetImage = imageNamesPerSubject{i}{indexPetCtAcImage};
        dcmTagsPet{i} = getfield(dcmTags{i}.h, namePetImage);
        niftiPetCtAcFilenames{i} = [niftiPathThisSubject namePetImage '.nii.gz'];

        %fMRI_voxelSize_mm(i,:) = [dcmTagsRsFmri{i}.PixelSpacing' dcmTagsRsFmri{i}.SliceThickness];
        %fMRI_matrixSize_voxels(i,:) = [dcmTagsRsFmri{i}.AcquisitionMatrix(1) dcmTagsRsFmri{i}.AcquisitionMatrix(4) dcmTagsRsFmri{i}.LocationsInAcquisition]; 
        strTime = strsplit(dcmTagsPet{i}.AcquisitionDateTime, '.'); % In some cases a . is added at the end related to the series.
        acqDateTime(i) = datetime(strTime{1}, 'InputFormat', 'yyyyMMddHHmmss');
        corrections{i} = dcmTagsPet{i}.CorrectedImage;
        injectionDateTime(i) = datetime(dcmTagsPet{i}.RadiopharmaceuticalInformationSequence.Item_1.Unknown_0018_1078, 'InputFormat', 'yyyyMMddHHmmss');
        % Nor properly filled in in the dicom:
        if strcmp(casesToProcess{i}, 'CP0107')
            injectionDateTime(i) = datetime([dcmTagsPet{i}.RadiopharmaceuticalInformationSequence.Item_1.Unknown_0018_1078(1:8) '125500'], 'InputFormat', 'yyyyMMddHHmmss');
        end
        timeInjecttionToAcq(i) = acqDateTime(i)-injectionDateTime(i);
        % To check this field (activity or counts?):
        totalDose_Bq(i) = dcmTagsPet{i}.RadiopharmaceuticalInformationSequence.Item_1.Unknown_0018_1074;
        totalDose_mCi(i) = totalDose_Bq(i) ./37e6;
        frameDuration_msec(i) = dcmTagsPet{i}.ActualFrameDuration;
        frameDuration_sec(i) = frameDuration_msec(i)./1000;
        weight_kg(i) = dcmTagsPet{i}.PatientWeight;
        height_m(i) = dcmTagsPet{i}.PatientSize;
        promptsCounts(i) = dcmTagsPet{i}.Unknown_0054_1310;
        delayedCounts(i) = dcmTagsPet{i}.Unknown_0054_1311;
        scatterFraction(i) = dcmTagsPet{i}.Unknown_0054_1323;
        deadTimeFactor(i) = dcmTagsPet{i}.Unknown_0054_1324;
        info = niftiinfo([niftiPetCtAcFilenames{i}]);
        pet_imageSize_voxels(i,:) = info.ImageSize(1:3);
        pet_voxelSize_mm(i,:) = info.PixelDimensions(1:3);
        % Check if they have been already proccessed with the UNSAM PET
        % Quantification Tools:
        if ~exist([quantPath casesToProcess{i} '/'])
            indicesToRunUPQT = [indicesToRunUPQT i];
        end
        if strcmp(casesToProcess{i}, 'CP0107old')
            casesToProcess{i} = 'CP0107';
        end
    end
end
frameDuration_min = frameDuration_sec./60;
timeInjecttionToAcq_min = minutes(timeInjecttionToAcq);
f18HalfLife_min = 109.7;
activityAtAcquisitionTime_mCi = totalDose_mCi.*exp(-log(2)./f18HalfLife_min.*timeInjecttionToAcq_min);
activityAtAcquisitionTime_Bq = totalDose_Bq.*exp(-log(2)./f18HalfLife_min.*timeInjecttionToAcq_min);
%% SAVE DATA
save(strcat([preprocessedDataPath 'petInfo_' ], string(datetime('today','Format','y_MM_dd'))))

%% TABLE WITH ACQ PARAMETERS
tablePet = table(casesToProcess', pet_voxelSize_mm, pet_imageSize_voxels, ...
    totalDose_Bq', timeInjecttionToAcq_min', activityAtAcquisitionTime_Bq', frameDuration_sec', promptsCounts', delayedCounts',...
    scatterFraction', deadTimeFactor', weight_kg', height_m',...
    'VariableNames', {'ID', 'Voxel Size [mm]', 'Matrix Size [voxels]', 'Injected Activity [Bq]', 'Time Injection to Scan [min]',...
    'Activity at Scan Time [Bq]', 'Frame Duration [sec]', 'Prompts [counts]', 'Delayed [counts]', 'Scatter Fraction [%]',...
    'Dead Time Factor', 'Weight [kg]', 'Height [m]'});
writetable(tablePet, [preprocessedDataPath 'pet_parameters.csv' ]);

%% RUN UNSAM PET QUANTIFICATION TOOLS
% Call UNSAM PET quantification tools.
% python3 main.py path_PET path_MRI subject_name output_dir
if ~isempty(indicesToRunUPQT)
    % If in windows, save the script and run it later in Ubuntu for windows:
    fid = fopen(filenameUPQTScript,'w');
    fprintf(fid,'#!/bin/bash \n');
    % First lines to prepare the data:
    for i = indicesToRunUPQT
        %outputPathThisSubject = [quantPath  '/' casesToProcess{i} '/'];
        % The subfolder with the name is created by the tool
        % Check if the T1 weighted image is available:
        mriT1Filename = fullfile(preprocessedMriDataPath, 'Nifti', casesToProcess{i}, dfltNameT1);
        dirT1 = dir([mriT1Filename '*']);
        if ~isempty(dirT1)
            niftiMriT1Filename = fullfile(preprocessedMriDataPath, 'Nifti', ...
                casesToProcess{i},dirT1(end).name);
            % Check is freesurfer is available
            freesurferPath = fullfile(preprocessedMriDataPath, 'Freesurfer', casesToProcess{i});
            if exist(fullfile(freesurferPath, 'mri', 'aparc+aseg.mgz'))
                fprintf(fid, 'python3 %s/main.py -p %s -m %s -f %s -s %s -o %s\n', unsamPetQuantificationToolsPath, niftiPetCtAcFilenames{i}, ....
                    niftiMriT1Filename, freesurferPath, casesToProcess{i}, quantPath);
            else
                fprintf(fid, 'python3 %s/main.py -p %s -m %s -s %s -o %s\n', unsamPetQuantificationToolsPath, niftiPetCtAcFilenames{i}, ....
                    niftiMriT1Filename, casesToProcess{i}, quantPath);
            end
        else
            fprintf(fid, 'python3 %s/main.py -p %s -s %s -o %s\n', unsamPetQuantificationToolsPath, niftiPetCtAcFilenames{i}, ....
                casesToProcess{i}, quantPath);
        end
    end
    fclose(fid);
else
    disp('UNSAM Pet Quantification Tools analysis is already available for all the subjects.')
end
% %% COPY THE IMAGES AND THE MRI IMAGE TO A NEW SUBFOLDER
% mriSeqToCopy = {'t1_mprage_1x1x1_estric_', 't2_space_dark_fluid_sag_p3_iso_estric_', ...
%     'SWI_Images', 'Perfusion_Weighted','asl_3d_tra_iso_3_0_highres_estric_',...
%     'TENSOR_b1000_64_2x2x2_AP_estric_', 'TENSOR_b1000_64_2x2x2_AP_estric_ADC', 'TENSOR_b0_2x2x2_PA_estric_'};
% outputPetMriDataset = '/home/martin/data_imaging/CovidProject/Estudio/PET_MRI_Dataset/';
% for i = 1 : numel(casesToProcess)
%     outputPetMriDatasetThisSubject = [outputPetMriDataset casesToProcess{i}];
%     if ~exist(outputPetMriDatasetThisSubject)
%         mkdir(outputPetMriDatasetThisSubject)
%     end
%     % Copy mri
%     for j = 1 : numel(mriSeqToCopy)
%         mriFilename = fullfile(preprocessedMriDataPath, 'Nifti', casesToProcess{i}, [mriSeqToCopy{j} niftiExtension]);
%         if exist(mriFilename)
%             mriOutputFilename = fullfile(outputPetMriDatasetThisSubject, [casesToProcess{i}, '_', mriSeqToCopy{j}, niftiExtension]);
%             copyfile(mriFilename, mriOutputFilename);
%         end
%     end
%     % Copy t1:
%     petOutputFilename = fullfile(outputPetMriDatasetThisSubject, [casesToProcess{i}, '_fdg_pet', niftiExtension]);
%     copyfile(niftiPetCtAcFilenames{i}, petOutputFilename);
% end