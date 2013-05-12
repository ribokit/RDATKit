function values = get_tag( annotations, tag )
% value = get_tag( annotations, tag )
% value = get_tag( annotation, tag )
%
% for use in parsing data_annotations or annotations cells.
%

if ischar( annotations ) & iscell( tag ) % user has given in wrong order.
  tag_tmp = annotations;
  annotations = tag;
  tag = tag_tmp;
end
assert( ischar( tag ) );
assert( iscell( annotations ) );

if length( annotations ) > 0 & ischar( annotations{1} )
  values = get_tag_one_annotation( annotation, tag );
  return 
else
  for i = 1:length( annotations )
    values{i} = get_tag_one_annotation( annotations{i}, tag );
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%`
function value = get_tag_one_annotation( annotation, tag )

vals = find_annotation_tag( annotation, tag );
if length( vals ) == 0; 
  value = ''; 
else
  value = vals{1};
end