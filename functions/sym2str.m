function T = sym2str(S);
% function T = sym2str(S);
%
% Converts a symbolic matrix (S) into a string expression (T),
% that can be evaluated using eval()
%
% Created by: Peter de Lange - pdldelange@gmail.com

T = char(S);                % convert symbolic expression to string
if  strmatch('matrix',T);   % check if string is matrix
    T = T(9:end-2);         % remove some begin/end characters
end
T = strrep(T,'], [',';');    % replace bracket/comma by semicolon