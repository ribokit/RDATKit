function value = get_tag( tag, annotation )
% value = get_tag( tag, annotation )
%
% for use in parsing data_annotations or annotations cells.
%

vals = find_annotation_tag( annotation, tag );
if length( vals ) == 0; 
  value = ''; 
else
  value = vals{1};
end