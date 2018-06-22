function rdat = fill_sequences_and_structures( rdat );

% let's try to fill the "sequences" field if it isn't there.
rdat = fill_sequences_if_empty( rdat );
rdat = fill_structures_if_empty( rdat );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rdat = fill_sequences_if_empty( rdat );

if isempty( rdat.data_annotations ); return; end;
for i = 1:size(  rdat.reactivity, 2 )
    if ( i > length( rdat.sequences )  |  isempty(rdat.sequences{i} ) ) & i <= length( rdat.data_annotations )
        rdat.sequences{i} = rdat.sequence;
        
        data_annotation = rdat.data_annotations{i};
        for m = 1:length( data_annotation )
            
            c = str2cell( data_annotation{m},':' );
            
            if ~isempty(c) & strcmp( c{1}, 'sequence' )
                rdat.sequences{i} = c{2};
                continue;
            end
            
            if ~isempty(c) & strcmp( c{1}, 'mutation' )
                [mutpos,mut_seq] = get_mutation_info_from_tag( data_annotation{m}, rdat );
                if isnan( mutpos ) continue; end;
                rdat.sequences{i} = [rdat.sequence(1: (min(mutpos)-1) ), mut_seq, rdat.sequence( max(mutpos)+1 : end )];
            end
            
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rdat = fill_structures_if_empty( rdat );

if size( rdat.reactivity, 2 ) == 0; return; end;

for i = 1:size(  rdat.reactivity, 2 )
    
    if i > length( rdat.structures )  |  isempty(rdat.structures{i} );
        rdat.structures{i} = rdat.structure;
        
        data_annotation = rdat.data_annotations{i};
        for m = 1:length( data_annotation )
            
            c = str2cell( data_annotation{m},':' );
            if ~isempty(c) & strcmp( c{1}, 'structure' );
                rdat.structures{i} = c{2};
                continue;
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = str2cell(s, delim)

if ~exist( 'delim', 'var' );
    tabchar = sprintf( '\t' );
    if ~isempty(strfind( s, tabchar ) )
        delim = tabchar;
    else
        delim = ' ';
    end
end;

rest = s;
i = 1;
c = {};
while ~isempty(rest);
    [t, rest] = strtok(rest, delim);
    c{i} = t;
    i = i + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = cell2str( c, delim )
if  ~exist('delim', 'var'); delim = ' '; end;
s = '';
if ~isempty(c)
  s = c{1};;
  for i = 2:length(c); s = [s,delim,c{i}]; end
end
