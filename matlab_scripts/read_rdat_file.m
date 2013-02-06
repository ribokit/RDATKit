function rdat = read_rdat_file(filename)
%
%  rdat = parse_rdat(filename)
%
%   rdat is 'rdat' file format for high-throughput RNA capillary electrophoresis data
%
% Copyright P. Cordero, R. Das, Stanford University, 2010-2013.
%

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
rdat.mutpos           = [];
rdat.annotations      = {};
rdat.data_annotations = {};
rdat.reactivity        = [];
rdat.reactivity_error  = [];
rdat.trace            = [];
rdat.xsel             = [];
rdat.xsel_refine      = [];

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
      rdat.seqpos = strread( remove_tag(line, 'SEQPOS') );
    elseif strfind(line, 'MUTPOS') == 1
      rdat.mutpos = strread(remove_tag(strrep(line, 'WT', 'NaN'), 'MUTPOS'), '');
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
      line = strrep(line, 'XSEL ', '');
      rdat.xsel = strread( line );
    else % might be a blank line
      [t,r] = strtok( line,' ' );
      if length(r) > 0
	fprintf('\nError parsing file %s\n', line);
      end
    end 

end

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
function  d = fill_mutpos( d );

% need to figure out where mutations are based on tags like "mutation:G64C" in data_annotation
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

