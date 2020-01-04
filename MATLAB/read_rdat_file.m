function rdat = read_rdat_file(filename)
%
%  rdat = parse_rdat(filename)
%
%   rdat is 'rdat' file format for high-throughput RNA capillary electrophoresis data
%
% Copyright P. Cordero, R. Das, Stanford University, 2010-2013.
%

if nargin==0; help( mfilename ); return; end;

rdat           = RDATFile;
rdat.version          = 0.0;
rdat.comments         = {};
rdat.name             = '';
rdat.sequence         = '';
rdat.structure        = '';
rdat.sequences        = {};
rdat.structures       = {};
rdat.offset           =  0;
rdat.seqpos           = [];
rdat.annotations      = {};
rdat.data_annotations = {};
rdat.reactivity        = [];
rdat.reactivity_error  = [];
rdat.trace            = [];
rdat.xsel             = [];
rdat.xsel_refine      = [];

legacy_data_anot_sequences  = {};
legacy_data_anot_structures = {};

sequence_seqpos = '';

fprintf( 1, 'Parsing file from rdat: %s\n', filename );
fid = fopen(filename);

while 1
    line = fgets(fid);
    if line == -1
        break;
    end
    line = strrep(line, sprintf('\n'), '');
    if strfind(line, 'VERSION') == 1 
      rdat.version = remove_tag( line, 'VERSION');
    elseif strfind(line, 'RDAT_VERSION') == 1 
      rdat.version = remove_tag( line, 'RDAT_VERSION');
    elseif strfind(line, 'COMMENT') == 1
      rdat.comments = [ rdat.comments, remove_tag(line, 'COMMENT') ];
    elseif ~isempty(strfind(line, 'ANNOTATION')) && isempty(strfind(line, 'ANNOTATION_DATA')) && isempty(strfind(line, 'DATA_ANNOTATION'))
      rdat.annotations = remove_empty_cells( str2cell( remove_tag(line,'ANNOTATION') ) );
    elseif strfind(line, 'NAME') == 1
      rdat.name = remove_tag(line, 'NAME');
    elseif strfind(line, 'SEQUENCE') == 1
      cols = str2cell( remove_tag( line, 'SEQUENCE' ) );
      if length( cols ) > 1          
          rdat.sequences{str2num(cols{1})} = strjoin(cols(2:end));
      else  
          rdat.sequence = strrep(cols{1}, ' ','');
      end
    elseif strfind(line, 'OFFSET') == 1
      rdat.offset = str2num( remove_tag(line,'OFFSET'));
    elseif strfind(line, 'SEQPOS') == 1
      %rdat.seqpos = strread( remove_tag(line, 'SEQPOS') );
      [ rdat.seqpos, sequence_seqpos ] = get_seqpos( remove_tag(line, 'SEQPOS') );
    elseif strfind(line, 'MUTPOS') == 1
      
      %rdat.mutpos = strread(remove_tag(strrep(line, 'WT', 'NaN'), 'MUTPOS'), '');
      fprintf( 'No longer reading in MUTPOS\n' );
    elseif strfind(line, 'STRUCTURE') == 1
      cols = str2cell( remove_tag( line, 'STRUCTURE' ) );
      if length( cols ) > 1
          rdat.structures{str2num(cols{1})} = strjoin(cols(2:end));
      else  
          rdat.structure = strrep(cols{1}, ' ','');
      end
    elseif strfind(line, 'ANNOTATION_DATA') == 1
      line = remove_tag( line, 'ANNOTATION_DATA' );
      rdat = get_data_annotation( rdat, line );
    elseif strfind(line, 'DATA_ANNOTATION') == 1
      line = remove_tag( line, 'DATA_ANNOTATION' );
      rdat = get_data_annotation( rdat, line );
    elseif strfind(line, 'REACTIVITY_ERROR') == 1
      line = remove_tag( line, 'REACTIVITY_ERROR' );
      line_read = strread( line );
      rdat.reactivity_error(:, line_read(1) ) = line_read(2:end);
    elseif strfind(line, 'REACTIVITY') == 1
      line = remove_tag( line, 'REACTIVITY' );
      line_read = strread( line );
      rdat.reactivity(:, line_read(1) ) = line_read(2:end);
    elseif strfind(line, 'DATA_ERROR') == 1  % backwards compatibility
      line = remove_tag( line, 'DATA_ERROR' );
      line_read = strread( line );
      rdat.reactivity_error(:, line_read(1) ) = line_read(2:end);
    elseif strfind(line, 'DATA') == 1  % backwards compatibility
      line = remove_tag( line, 'DATA' );
      line_read = strread( line );
      if length( line_read ) >  length( rdat.seqpos )
          rdat.reactivity(:, line_read(1) ) = line_read((end-length(rdat.seqpos)+1):end);
      else
          rdat.reactivity(1:length(line_read)-1, line_read(1) ) = line_read(2:end);
      end
      
    elseif strfind(line, 'AREA_PEAK_ERROR') == 1 % backwards compatibility
      line = remove_tag( line, 'AREA_PEAK_ERROR' );
      line_read = strread( line );
      rdat.reactivity_error(:, line_read(1) ) = line_read(2:end);
    elseif strfind(line, 'AREA_PEAK') == 1 % backwards compatibility
      line = remove_tag( line, 'AREA_PEAK' );
      line_read = strread( line );
      rdat.reactivity(:, line_read(1) ) = line_read(2:end);
    elseif strfind(line, 'TRACE') == 1
      line = remove_tag( line, 'TRACE' );
      line_read = strread( line );
      rdat.trace(:, line_read(1) ) = line_read(2:end);
    elseif strfind(line, 'XSEL_REFINE') == 1
      line = remove_tag( line, 'XSEL_REFINE' );
      line_read = strread( line );
      rdat.xsel_refine(:,line_read(1)) = line_read(2:end);
    elseif strfind(line, 'XSEL') == 1
      line = remove_tag( line, 'XSEL' );
      rdat.xsel = strread( line );
    else % might be a blank line
      [t,r] = strtok( line,' ' );
      if length(r) > 0
	fprintf('\nError parsing file %s\n', line);
      end
    end 
end
if ~isempty( legacy_data_anot_sequences ) rdat.sequences = legacy_data_anot_sequences; end
if ~isempty( legacy_data_anot_structures ) rdat.structures = legacy_data_anot_structures; end
rdat = fill_data_annotations_if_empty( rdat );
rdat = fill_sequences_and_structures( rdat );


if length(rdat.sequence) == 0 
 if isempty(rdat.sequences)
  fprintf( 'No sequences detected or sequence indices do not start at one' );
 else
  rdat.sequence = rdat.sequences{1};
 end
end
if length(rdat.structure) == 0
 if isempty(rdat.structures)
  fprintf( 'No structures detected or structure indices do not start at one' );
 else
  rdat.structure = rdat.structures{1};
 end
end
%if length(rdat.reactivity_error) == 0
%    rdat.reactivity_error = fill_reactivity_error_from_reads( rdat );
%end


% output a warning of the sequence characters in 'SEQPOS' don't match up with the given sequence...
check_sequence_seqpos( sequence_seqpos, rdat.seqpos, rdat.sequence, rdat.offset );

if size( rdat.trace, 2 ) > 0; fprintf( 'Number of traces         : %d \n', size( rdat.trace, 2 ) ); end;
fprintf( 'Number of reactivity lines: %d \n', size( rdat.reactivity, 2 ) );
fclose( fid );

check_rdat( rdat );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%show_rdat_data(rdat)
function c = str2cell(s, delim)

if ~exist( 'delim' ) 
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
while length(rest)
    [t, rest] = strtok(rest, delim);
    c{i} = t;
    i = i + 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = cell2str(c, delim)
if length(c) == 0; s = ''; return; end;
s = c{1};
for i = 2:length(c); s = [s, delim, c{i} ]; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  line = remove_tag( line, tag );

delim = ' ';

% does this line include tabs between fields or spaces?
tabchar = sprintf( '\t' );
if ~isempty(strfind( line, tabchar ) );  delim = tabchar; end;

if strfind( line, [tag,':'] ) % new format v0.23
  line = strrep(line, [tag,':'],'');
else
  line = strrep(line, [tag,delim],'');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rdat = fill_data_annotations_if_empty( rdat );

for i = 1:size(  rdat.reactivity, 2 )
  if i > length( rdat.data_annotations ) 
    rdat.data_annotations{i} = {};
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [seqpos, sequence_seqpos] = get_seqpos( seqpos_info );

seqpos_tags = str2cell( seqpos_info );
sequence_seqpos = '';

seqpos = [];

for i = 1:length( seqpos_tags )

  tag = seqpos_tags{i};

  if isempty( str2num( tag(1) ) ) & tag(1) ~= '-' % first letter is a character not a number of minus sign.
    sequence_seqpos = [sequence_seqpos, tag(1) ]; 
    tag = tag(2:end);
  end

  seqpos = [seqpos, str2num( tag )];

end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c_new = remove_empty_cells( c );
c_new = {};
for i = 1:length( c )
    if length( c{i} ) > 0 
        c_new = [ c_new, c{i} ];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c_new = remove_sequence_structure_fields_if_number( c );
c_new = {};
for i = 1:length( c )
    if length( c{i} ) > 0 
        cols = str2cell( c{i},':' );
        if length( cols ) > 1 & strcmp( cols{1}, 'sequence' ) & str2num( cols{2} ) ~= 0; continue; end;
        if length( cols ) > 1 & strcmp( cols{1}, 'structure' ) & str2num( cols{2} ) ~= 0; continue; end;
         c_new = [ c_new, c{i} ];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output a warning of the sequence characters in 'SEQPOS' don't match up with the given sequence...
function ok = check_sequence_seqpos( sequence_seqpos, seqpos, sequence, offset );

ok = 1;

if length( sequence_seqpos ) == 0; return; end;

if length( sequence_seqpos ) ~= length( seqpos ); 
  fprintf( 'Number of characters in sequence_seqpos %d, does not match number of seqpos %d\n', length( sequence_seqpos) , length( seqpos ) ); 
  ok = 0;
  return;
end

for i = 1:length( sequence_seqpos )
  c1 = lower( sequence_seqpos(i) );
  m = seqpos(i) - offset;
  if ( m < 1 | m > length( sequence ) )
    fprintf( 'Warning: seqpos %d is not inside sequence, given offset %d\n', seqpos(i), offset );
    ok = 0;
    continue;
  end
  c2 = lower( sequence( m ) );
  if ( c1 ~= 'X' & c2 ~= 'X' & c1 ~= c2 )
    fprintf( 'Warning: mismatch at seqpos %d, between SEQPOS nucleotide %s and SEQUENCE nucleotide %s\n', seqpos(i), sequence_seqpos(i),  sequence( seqpos(i) - offset ) );
    ok = 0;
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rdat = get_data_annotation( rdat, line )
cols = str2cell( line );
idx = str2num( cols{1} );
if isempty( idx ) % weird edge case where str2cell fails
    firstcols = str2cell(cols{1});
    cols = [firstcols,cols(2:end)];
    idx = str2num( cols{1} );
end
anot = cols(2:end);
rdat.data_annotations{idx} = remove_empty_cells( anot );
% look for a 'sequence' tag.
for j = 1:length( anot )
    if ~isempty( strfind( anot{j}, 'sequence:' ) )
        [dummy, r ] = strtok( anot{j}, 'sequence:' );
        if ~isnan( str2double( dummy ) )
            assert( length( rdat.sequences ) >= str2num(dummy) );
            legacy_data_anot_sequences{idx} = rdat.sequences{str2num(dummy)}; % old-style format.
        else
            rdat.sequences{idx} = dummy;
        end
    end
end
for j = 1:length( anot )
    if ~isempty( strfind( anot{j}, 'structure:' ) )
        [dummy, r ] = strtok( anot{j}, 'structure:' );
        if length( rdat.structures ) == 0 & length( rdat.structure ) > 0; rdat.structures{1} = rdat.structure; end;
        if ~isnan( str2double( dummy ) )
            assert( length( rdat.structures ) >= str2num(dummy) );
            legacy_data_anot_structures{idx} = rdat.structures{str2num(dummy)}; % old-style format.
        else
            rdat.structures{idx} = dummy;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function reactivity_error = fill_reactivity_error_from_reads( rdat )
% estimate errors based on Poisson [sqrt(reads)]
% this is not actually exact --
reactivity_error = [];
for i = 1:size( rdat.reactivity, 2)
    anot = rdat.data_annotations{i};
    reads = get_tag( anot, 'reads' );
    if isempty( reads ); 
        reactivity_error(1:size(rdat.reactivity,1),i) = NaN;
        continue;
    end
    reads = str2num(reads);
    if reads > 0
        rdat.reactivity(:,i) * sqrt(reads );
        reactivity_error = sqrt(rdat.reactivity(:,i))/ sqrt(reads);
        fprintf( 'WARNING! WARNING! Trying to fill in reactivity error from reads and reactivity -- this is not exact.' );
    else
        reactivity_error(1:size(rdat.reactivity,1),i) = 0.0;
    end
end






