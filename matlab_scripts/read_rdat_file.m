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
    elseif ~isempty(strfind(line, 'ANNOTATION')) && isempty(strfind(line, 'ANNOTATION_DATA'))
      rdat.annotations = str2cell( remove_tag(line,'ANNOTATION') );
    elseif strfind(line, 'NAME') == 1
      rdat.name = remove_tag(line, 'NAME');
    elseif strfind(line, 'SEQUENCE') == 1
      cols = str2cell( remove_tag( line, 'SEQUENCE' ) );
      if length( cols ) > 1
	rdat.sequences{str2num(cols{1})} = cols(2:end);
      else  
	rdat.sequence = strrep(cols{1}, ' ','');
      end
    elseif strfind(line, 'OFFSET') == 1
      rdat.offset = str2num( remove_tag(line,'OFFSET'));
    elseif strfind(line, 'SEQPOS') == 1
      %rdat.seqpos = strread( remove_tag(line, 'SEQPOS') );
      [ rdat.seqpos, sequence_seqpos ] = get_seqpos( remove_tag(line, 'SEQPOS') );
    %elseif strfind(line, 'MUTPOS') == 1
    %  rdat.mutpos = strread(remove_tag(strrep(line, 'WT', 'NaN'), 'MUTPOS'), '');
    elseif strfind(line, 'STRUCTURE') == 1
      cols = str2cell( remove_tag( line, 'STRUCTURE' ) );
      if length( cols ) > 1
	rdat.structures{str2num(cols{1})} = cols(2:end);
      else  
	rdat.structure = strrep(cols{1}, ' ','');
      end
    elseif strfind(line, 'ANNOTATION_DATA') == 1
      line = remove_tag( line, 'ANNOTATION_DATA' );
      cols = str2cell( line );
      idx = str2num( cols{1} );
      anot = cols(2:end);
      rdat.data_annotations{idx} = anot;
      % look for a 'sequence' tag.
      for j = 1:length( anot )
	if ~isempty( strfind( anot{j}, 'sequence:' ) )
	  [dummy, r ] = strtok( anot{j}, 'sequence:' );
	  rdat.sequences{idx} = dummy;
	end
      end
      for j = 1:length( anot )
	if ~isempty( strfind( anot{j}, 'structure:' ) )
	  [dummy, r ] = strtok( anot{j}, 'structure:' );
	  rdat.structures{idx} = dummy;
	end
      end
    elseif strfind(line, 'REACTIVITY_ERROR') == 1
      line = remove_tag( line, 'REACTIVITY_ERROR' );
      line_read = strread( line );
      rdat.reactivity_error(:, line_read(1) ) = line_read(2:end);
    elseif strfind(line, 'REACTIVITY') == 1
      line = remove_tag( line, 'REACTIVITY' );
      line_read = strread( line );
      rdat.reactivity(:, line_read(1) ) = line_read(2:end);
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

% let's try to fill the "sequences" field if it isn't there.
rdat = fill_sequences_if_empty( rdat );
rdat = fill_structures_if_empty( rdat );

if strcmp(rdat.sequence, '')
 if isempty(rdat.sequences)
  fprintf( 'No sequences detected or sequence indices do not start at one' );
 else
  rdat.sequence = rdat.sequences{1};
 end
end
if strcmp(rdat.structure, '')
 if isempty(rdat.structures)
  fprintf( 'No structures detected or structure indices do not start at one' );
 else
  rdat.structure = rdat.structures{1};
 end
end

% output a warning of the sequence characters in 'SEQPOS' don't match up with the given sequence...
check_sequence_seqpos( sequence_seqpos, rdat.seqpos, rdat.sequence, rdat.offset );

fprintf( 'Number of traces         : %d \n', size( rdat.trace, 2 ) );
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% no longer in use
function  d = fill_mutpos( d );

% need to figure out where mutations are based on tags like "mutation:G64C" in data_annotation
% Not actually in use...
d.mutpos = []
for k = 1:length( d.data_annotations )
  d.mutpos(k) = NaN;
  data_annotation = d.data_annotations{k};
  for m = 1:length( data_annotation )
    c = str2cell( data_annotation{m},':' );
    if length(c)> 0 & strcmp( c{1}, 'mutation' )
      num = str2num( remove_AGCTU( c{2} ) );
      if length(num)>0;  d.mutpos(k) =  num; end;
    end
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = remove_AGCTU( x )
x = strrep( x, 'A', '');
x = strrep( x, 'G', '');
x = strrep( x, 'C', '');
x = strrep( x, 'T', '');
x = strrep( x, 'U', '');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rdat = fill_sequences_if_empty( rdat );

if length( rdat.data_annotations ) == 0; return; end;

for i = 1:size(  rdat.reactivity, 2 )
  if ( i > length( rdat.sequences )  |  length(rdat.sequences{i} ) == 0 ) & i <= length( rdat.data_annotations )
    rdat.sequences{i} = rdat.sequence;
    
    data_annotation = rdat.data_annotations{i};
    for m = 1:length( data_annotation )
      c = str2cell( data_annotation{m},':' );
      if length(c) > 0  & strcmp( c{1}, 'mutation' )
      
	start_seq = '';
	mut_seq = '';
	mut_num = '';
	tag = cell2str( c(2:end), ':' );
      
	q = 1;	
	while ( q <= length(tag)  & isempty( str2num( tag(q) ) ) & tag(q)~='(' )
	  start_seq = [start_seq, tag(q) ];
	  q = q+1;
	end
	if q <= length(tag) & tag(q) == '('; q = q+1; end;
	while ( q <= length( tag ) &  (~isempty( str2num( tag(q) ) ) | tag(q)==':')  )
	  mut_num = [mut_num, tag(q) ];
	  q = q+1;
	end
	if q <= length(tag) & tag(q) == ')'; q = q+1; end;
	while q <= length( tag );
	  mut_seq = [mut_seq, tag(q) ];
	  q = q+1;
	end
	
	if  length( mut_num ) == 0 
	  if ~strcmp( start_seq, 'WT' )
	    fprintf( 'WARNING! Could not find mutation position in mutation annotation: %s\n', tag );
	  end
	  continue;
	end
	if  length( start_seq ) == 0 
	  fprintf( 'WARNING! Could not find starting nucleotide in mutation annotation: %s\n', tag );
	end
	if  length( mut_seq ) == 0 
	  fprintf( 'WARNING! Could not find mutation nucleotide in mutation annotation: %s\n', tag );
	end
	
	mutpos = str2num( mut_num ) - rdat.offset;
	if ( ~strcmp( rdat.sequence( mutpos ), start_seq ) )
	  mutpos = str2num( mut_num ); % perhaps specified without offset...
	    fprintf( 'WARNING! Mismatch between mutation nucleotides: %s (mutation annotation) vs. %s. (sequence)\n', rdat.sequence(mutpos), start_seq );
	  if ( strcmp( rdat.sequence( mutpos ), start_seq ) )
	    fprintf( 'OK, specified mutpos without taking into account offset...\n' );
	  end
	end
	rdat.sequences{i} = [rdat.sequence(1: (min(mutpos)-1) ), mut_seq, rdat.sequence( max(mutpos)+1 : end )];
      end
      
    end
  end
end
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rdat = fill_structures_if_empty( rdat );

if length( rdat.data_annotations ) == 0; return; end;

for i = 1:size(  rdat.reactivity, 2 )
  if i > length( rdat.structures )  |  length(rdat.structures{i} ) == 0
    rdat.structures{i} = rdat.structure;
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [seqpos, sequence_seqpos] = get_seqpos( seqpos_info );

seqpos_tags = str2cell( seqpos_info );
sequence_seqpos = '';

seqpos = [];

for i = 1:length( seqpos_tags )

  tag = seqpos_tags{i};

  if isempty( str2num( tag(1) ) ) % first letter is a character
    sequence_seqpos = [sequence_seqpos, tag(1) ]; 
    tag = tag(2:end);
  end

  seqpos = [seqpos, str2num( tag )];

end

return;


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
  c2 = lower( sequence( seqpos(i) - offset ) );
  if ( c1 ~= 'X' & c2 ~= 'X' & c1 ~= c2 )
    fprintf( 'Warning: mismatch at seqpos %d, between SEQPOS nucleotide %s and SEQUENCE nucleotide %s\n', seqpos(i), sequence_seqpos(i),  sequence( seqpos(i) - offset ) );
    ok = 0;
  end
end


