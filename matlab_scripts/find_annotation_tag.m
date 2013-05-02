function value = find_annotation_tag( annotation, tag )
% value = find_annotation_tag( tag )
%
% for use in parsing data_annotations or annotations cells.
%

value = {};
for m = 1:length( annotation )
  anot = annotation{m};
  idx = strfind( anot, tag );
  if ~isempty( idx ) & idx(1) == 1
    value = [value, anot( length(tag)+2: end )];
  end
end