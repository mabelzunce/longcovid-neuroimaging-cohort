% Reads dicom data exported on hard drive and with a number ID instad of
% name, and renames each folder based on name.
clear all
close all

dataPartitionPath = '/data/'; %'D:/'
imagesPartitionPath = '/data_imaging/'; %'F:/'
% % Load data:
% loadData = 0;
% Overwrite Nifti:
overwriteNifti = 1; % 0: if nifti conversion exists, just reads the images and headers, 1: covnerts and overwrites again
%% ADD PATHS
addpath([dataPartitionPath 'UNSAM/Tools/dicm2nii/'])
unsamPetQuantificationToolsPath = '/data/UNSAM/PET/UNSAMPETQuantificationTools/';

%% PATHS AND FILENAMES
niftiExtension = '.nii.gz';
dcmHeadersFilename = 'dcmHeaders.mat';
dicomDataPath = [imagesPartitionPath '/CEUNIM/PETdata/AllControlsAndCovid/DICOM/'];

%% NAME SEQUENCES
dfltNamePetCtAc = 'x_BR_CTAC_Brain';
tagNamePetCtAc = 'CTAC'; % If the default name is not found, use this tag to identify the image.
dfltNameT1 = 't1_mprage_1x1x1';
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
    dirPaths = dirPaths([dirPaths(:).isdir])
    casesToProcess = {dirPaths(3:end).name};
end

for i = 1 : numel(casesToProcess)
    % READ DICOM COLLECTION TO GET NAMES
    dicomPathThisSubject = [dicomDataPath casesToProcess{i}];
    tableDicomCollection = dicomCollection(dicomPathThisSubject);
    patientName = tableDicomCollection.PatientName{1};
    newDicomPathThisSubject = [dicomDataPath patientName];
    movefile(dicomPathThisSubject, newDicomPathThisSubject)
end

%% RENAME THE FOLDERS MANUALLY WITH THE IDS

%% ANONYMIZE DATA
dirPaths = dir(dicomDataPath);
dirPaths = dirPaths([dirPaths(:).isdir]);
casesToProcess = {dirPaths(3:end).name};
for i = 1 : 2%numel(casesToProcess)
    % READ DICOM COLLECTION TO GET NAMES
    values.PatientName = casesToProcess{i};
    dicomPathThisSubject = [dicomDataPath casesToProcess{i}];
    files = dir(dicomPathThisSubject);
    for p = 3:numel(files)
        if files(p).isdir
            files2 = dir(fullfile(dicomPathThisSubject, files(p).name));
            for q = 3:numel(files2)
                dicomanon(fullfile(dicomPathThisSubject, files(p).name, files2(q).name), ...
                    fullfile(dicomPathThisSubject, files(p).name, files2(q).name), ...
		            "update",values, "keep", ["PatientAge", "StudyID", "PatientSex","StudyDescription",...
                    "ImageIndex", "StudyInstanceUID", "SeriesInstanceUID"])
            end
        else
	        dicomanon(fullfile(dicomPathThisSubject, files(p).name), ...
                fullfile(dicomPathThisSubject, files(p).name), 'WritePrivate', true, ...
		        "update",values, "keep", ["PatientAge", "StudyID", "PatientSex","StudyDescription",...
                "ImageIndex", "StudyInstanceUID", "SeriesInstanceUID"])
        end
    end
end
