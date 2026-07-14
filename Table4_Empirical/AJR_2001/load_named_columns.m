function D = load_named_columns(xlsxFile)
% LOAD_NAMED_COLUMNS  Read a data sheet and return a struct of named columns.
%
%   D = LOAD_NAMED_COLUMNS('Table4_Data.xlsx') reads the spreadsheet and
%   returns a struct D whose field names are the variable names in row 1 of
%   the sheet, e.g. D.logpgp95, D.avexpr, D.logem4, ...
%
%   This replaces the old "xlsread + assignin('base',...)" idiom, which
%   silently dumped variables into the base workspace. Collecting them in a
%   struct keeps the workspace clean and makes the specifications below
%   self-documenting (D.avexpr is clearly "the avexpr column").

    [Data, ~, rawData] = xlsread(xlsxFile);
    % Align the header names to the numeric data columns: xlsread drops the
    % leading text columns (a row index, and sometimes a country name) from
    % Data, so the variable names are the LAST size(Data,2) header cells.
    names = rawData(1, end - size(Data, 2) + 1 : end);
    D = struct();
    for i = 1:numel(names)
        v = strtrim(char(names{i}));
        if ~isempty(v)
            D.(v) = Data(:, i);
        end
    end
end
