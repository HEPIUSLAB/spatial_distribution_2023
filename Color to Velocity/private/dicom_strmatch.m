% Moved to this folder for function dicomread_debugged.m without edits
%
%
% Moved by: Denis Routkevitch (droutke1@jhmi.edu)
% Last accessed: 1/24/2021

function idx = dicom_strmatch(str, cellOfStrings)
%DICOM_STRMATCH   Find substrings at beginning of larger string.
%   IDX = DICOM_STRMATCH(STR, CELLOFSTRINGS) looks and acts like
%   STRMATCH but isn't STRMATCH.

% Copyright 2011 The MathWorks, Inc.

idx = find(strncmpi(str, cellOfStrings, numel(str)));
