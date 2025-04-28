function annotations = convert_struct_to_annotations( s, pretag );
% annotations = convert_struct_to_annotations( s, pretag );
%
% Inputs
%  s =  struct with fields to convert to annotations
%  pretag = [optional] string with prefix for annotations
%
% (C) R. Das, 2023

if nargin == 0; help( mfilename ); return; end;
if ~exist('pretag','var') pretag = ''; end;
if length(pretag)>0 & ~strcmp(pretag(end),':'); pretag = [pretag,':']; end;

annotations = {};
if isempty(s); return; end;

assert(isstruct(s));
tags = fields( s );

for n = 1:length( tags )
    tag = tags{n};
    val = getfield(s,tag);
    if isnumeric(val) 
        if (length(val) == 1)
            if isa(val,'uint32')
                annotations = [annotations, {sprintf('%s%s:%d',pretag,tag,val)} ];
            else
                annotations = [annotations, {sprintf('%s%s:%f',pretag,tag,val)} ];
            end
        else
            assert( length(val) > 1);
            annotation = sprintf('%s%s:%f',pretag,tag,val(1));
            for n = 2:length(val)
                annotation = [annotation, sprintf(',%f',val(n)) ];
            end
            annotations = [annotations,annotation];
        end
    elseif ischar(val)
        annotations = [annotations, {sprintf('%s%s:%s',pretag,tag,val)} ];
    elseif iscell(val) 
        assert( ischar(val{1}) );
        annotations = [annotations, {sprintf('%s%s:%s',pretag,tag,strjoin(val,','))} ];
    elseif isstruct(val)
        annotations = [annotations, convert_struct_to_annotations(val,[pretag,tag,':'])];
    end
end