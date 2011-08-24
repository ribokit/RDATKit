function rdat = read_rdat_file(filename)
%
%  rdat = parse_rdat(filename)
%
%   rdat is 'rdat' file format for high-throughput RNA capillary electrophoresis data
%
% Copyright P. Cordero, R. Das, Stanford University, 2010.
%

rdat           = RDATFile;
rdat.version          = 0.0;
rdat.comments         = {};
rdat.name             = '';
rdat.sequence         = '';
rdat.structure        = '';
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
    line = strrep(line, '\n', '');
    if strfind(line, 'VERSION') > 0
      rdat.version = strrep( strrep(line(1:end-1), 'RDAT_VERSION ',''),  'VERSION ','');
     elseif strfind(line, 'COMMENT') > 0
      rdat.comments = [ rdat.comments, strrep(line(1:end-1), 'COMMENT ','') ];
    elseif ~isempty(strfind(line, 'ANNOTATION')) && isempty(strfind(line, 'ANNOTATION_DATA'))
      rdat.annotations = str2cell( strrep(line(1:end-1),'ANNOTATION','') );
    elseif strfind(line, 'NAME') > 0
      rdat.name = strrep(line, 'NAME ', '');
      rdat.name = rdat.name( 1:end-1 );% remove endline.
    elseif strfind(line, 'SEQUENCE') > 0
      rdat.sequence = strrep(strrep(line( 1:end-1 ), 'SEQUENCE ',''), ' ','');
    elseif strfind(line, 'OFFSET') > 0
      rdat.offset = str2num(strrep(line, 'OFFSET ',''));
    elseif strfind(line, 'SEQPOS') > 0
      rdat.seqpos = strread(strrep(line, 'SEQPOS ','') );
    elseif strfind(line, 'MUTPOS') > 0
      rdat.mutpos = strread(strrep(strrep(line, 'WT', 'NaN'), 'MUTPOS ',''), '');
    elseif strfind(line, 'STRUCTURE') > 0
      rdat.structure = strrep(strrep(line(1:end-1), 'STRUCTURE ',''), ' ', '');
    elseif strfind(line, 'ANNOTATION_DATA') > 0
      line = remove_tag( line, 'ANNOTATION_DATA' );
      [t,r] = strtok( line, ' ' );
      idx = str2num(t);
      rdat.data_annotations{idx} = str2cell( r );
    elseif strfind(line, 'REACTIVITY_ERROR') > 0
      line = remove_tag( line, 'REACTIVITY_ERROR' );
      line_read = strread( line );
      rdat.reactivity_error(:, line_read(1) ) = line_read(2:end);
    elseif strfind(line, 'REACTIVITY') > 0
      line = remove_tag( line, 'REACTIVITY' );
      line_read = strread( line );
      rdat.reactivity(:, line_read(1) ) = line_read(2:end);
    elseif strfind(line, 'AREA_PEAK_ERROR') > 0 % backwards compatibility
      line = remove_tag( line, 'AREA_PEAK_ERROR' );
      line_read = strread( line );
      rdat.reactivity_error(:, line_read(1) ) = line_read(2:end);
    elseif strfind(line, 'AREA_PEAK') > 0 % backwards compatibility
      line = remove_tag( line, 'AREA_PEAK' );
      line_read = strread( line );
      rdat.reactivity(:, line_read(1) ) = line_read(2:end);
    elseif strfind(line, 'TRACE') > 0
      line = remove_tag( line, 'TRACE' );
      line_read = strread( line );
      rdat.trace(:, line_read(1) ) = line_read(2:end);
    elseif strfind(line, 'XSEL_REFINE') > 0
      line = remove_tag( line, 'XSEL_REFINE' );
      line_read = strread( line );
      rdat.xsel_refine(:,line_read(1)) = line_read(2:end);
    elseif strfind(line, 'XSEL') > 0
      line = strrep(line, 'XSEL ', '');
      rdat.xsel = strread( line );
    else % might be a blank line
      [t,r] = strtok( line,' ' );
      if length(r) > 0
	fprintf('Error parsing file %s', line);
      end
    end 

end

fprintf( 'Number of traces         : %d \n', size( rdat.trace, 2 ) );
fprintf( 'Number of reactivity lines: %d \n', size( rdat.reactivity, 2 ) );
fclose( fid );

check_rdat( rdat );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%show_rdat_data(rdat)
function c = str2cell(s, delim)
if ~exist( 'delim' ) delim = ' '; end;
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
if strfind( line, [tag,':'] ) % new format v0.23
  line = strrep(line(1:end-1), [tag,':'],'');
else
  line = strrep(line(1:end-1), [tag,' '],'');
end
