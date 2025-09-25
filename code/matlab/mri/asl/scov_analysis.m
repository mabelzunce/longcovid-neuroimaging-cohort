clear all
close all
dataPartitionPath = '/data/'; %'D:/'
imagingPartitionPath = '/data_imaging/'; %'F:/'
%% PATHS AND FILENAMES
studyDataPath = [dataPartitionPath '/UNSAM/CovidProject2/'];
resultsPath = [studyDataPath '/DataAnalysis/'];
responseExcelFilename = fullfile(studyDataPath, 'Respuestas.xlsx');
summaryExcelFilename = fullfile(studyDataPath, 'ResumenRespuestas.xlsx');
volunteersExcelFilename = fullfile(studyDataPath, 'VoluntariosProyectoCovidProlongado.xlsx');
%% LOAD ASL CSVs
aslCsvPaths = '/home/martin/data/UNSAM/CovidProject2/DataAnalysis/ASL/CleanData/';
% total scov:
scovGmFilename = fullfile(aslCsvPaths, 'scov_total_gm_clean.csv');

scovTable = readtable(scovGmFilename);

%% Match groups
summaryTable = readtable(summaryExcelFilename,'Sheet', 'resumenTotal','ReadVariableNames', true, 'ReadRowNames', false);
% Match the demographics with scov:
indicesScov = [];
indicesNonCov = [];
for i = 1 : numel(scovTable.participant_id)
    ind=find(strcmp(scovTable.participant_id{i}, summaryTable.ID)>0);
    if ~isempty(ind)
       indicesScov = [indicesScov ind];
    else
        warning(sprintf('Subjects %s not found', scovTable.ID{i}));
        indicesNonCov = [indicesNonCov i];
    end
end
summaryTable = summaryTable(indicesScov,:);
group = summaryTable.Grupo;
scovTable.group = group;
scovTable.sex = summaryTable.Genero;
scovTable.age = summaryTable.Edad;
%% SIMPLE OVERSAMPLING
indicesControl = strcmp(scovTable.group,'CONTROL');
matToConcate = repmat(scovTable(indicesControl,:),3,1);
scovTable = [scovTable;matToConcate];
indicesShuffle = randperm(size(scovTable,1));
% Sigo en python
% = scovTable(indicesShuffle,:);
%groupShuffled = scovTableShuffled.group;