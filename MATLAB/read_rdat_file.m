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
    elseif strfind(line, 'COMMENT') > 0
      rdat.comments = [ rdat.comments, remove_tag(line, 'COMMENT') ];
    elseif ~isempty(strfind(line, 'ANNOTATION')) & ...
            isempty(strfind(line, 'ANNOTATION_DATA')) & ....
            isempty(strfind(line, 'DATA_ANNOTATION'))
        if ~contains(line,sprintf('\t')) & contains(line,'pH '); line = strrep(line,'pH ','pH'); end; %edge case, RNASEP_SHP_0000.rdat
        rdat.annotations = strip(remove_empty_cells( str2cell( remove_tag(line,'ANNOTATION') )));
    elseif strfind(line, 'NAME') == 1
      rdat.name = remove_tag(line, 'NAME');
    elseif strfind(line, 'SEQUENCE') == 1
      cols = strip(str2cell( remove_tag( line, 'SEQUENCE' ) ));
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
      warning( 'No longer reading in MUTPOS' );
    elseif strfind(line, 'STRUCTURE') == 1
      cols = str2cell( remove_tag( line, 'STRUCTURE' ) );
      if length( cols ) > 1
          rdat.structures{str2num(cols{1})} = strjoin(cols(2:end));
      else  
          if length(cols) > 0
              rdat.structure = strrep(cols{1}, ' ','');
          else
              rdat.structure = repmat('.',1,length(rdat.sequence));
          end
      end
    elseif strfind(line, 'ANNOTATION_DATA') == 1
      line = remove_tag( line, 'ANNOTATION_DATA' );
      rdat = get_data_annotation( rdat, line );
    elseif strfind(line, 'DATA_ANNOTATION') == 1
      line = remove_tag( line, 'DATA_ANNOTATION' );
      rdat = get_data_annotation( rdat, line );
    elseif strfind(line, 'REACTIVITY_ERROR') == 1
      line = remove_tag( line, 'REACTIVITY_ERROR' );
      line_read = strread( line, '%f' );
      nval = size(rdat.reactivity_error,1);
      if nval & length(line_read)-1 ~= nval;  % edge case, e.g.,_H101_0001.rdat
          warning(sprintf('Number of values in reactivity_error line %d != reactivity_error length %d!',length(line_read)-1,nval));
          line_read = line_read(1:(nval+1)); 
      end 
      rdat.reactivity_error(:, line_read(1) ) = line_read(2:end);
    elseif strfind(line, 'REACTIVITY') == 1
      line = remove_tag( line, 'REACTIVITY' );
      line_read = strread( line,'%f');
      nval = size(rdat.reactivity,1);
      if nval & length(line_read)-1 ~= nval;  % edge case, e.g.,_H101_0001.rdat
          warning(sprintf('Number of values in reactivity line %d != reactivity length %d!',length(line_read)-1,nval));
          line_read = line_read(1:(nval+1)); 
      end  
      rdat.reactivity(:, line_read(1) ) = line_read(2:end);
    elseif strfind(line, 'DATA_ERROR') == 1  % backwards compatibility
      line = remove_tag( line, 'DATA_ERROR' );
      line_read = strread( line );
      rdat.reactivity_error(:, line_read(1) ) = line_read(2:end);
    elseif strfind(line, 'DATA') == 1  % backwards compatibility
      line = remove_tag( line, 'DATA' );
      line_read = strread( line );
      if length( line_read ) >  length( rdat.seqpos ) 
          n = line_read(1);
          rdat.reactivity(1:length(rdat.seqpos), n) = line_read((end-length(rdat.seqpos)+1):end);
      else
          rdat.reactivity(1:length(line_read)-1, line_read(1) ) = line_read(2:end);
          rdat.reactivity(length(line_read):length(rdat.seqpos), line_read(1) ) = 0;
      end      
      rdat = check_for_legacy_datatype_reads(rdat, line_read);
    elseif strfind(line, 'READS') == 1  % backwards compatibility
      line = remove_tag( line, 'READS' );
      line_read = strread( line );
      reads(:,line_read(1)) = line_read(2:end);
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

% backwards compatibility with old Lucks format. Could do a better job
%  by actually recovering reactivity from reads and tracking Poisson
%  errors.
if all(rdat.reactivity_error(:)==0) & exist('reads','var') & size(reads,1)==size(rdat.reactivity,1)+1 & size(reads,2)==size(rdat.reactivity,2)
    warning('Trying to figure out errors based on READS lines!');
    rdat.reactivity_error = 0*rdat.reactivity;
    for i = 1:2:size(reads,2) % assume lines alternate between signal and noise...
        rdat.reactivity_error(:,i) = estimate_error_from_reads(rdat.reactivity(:,i), reads(:,i),reads(:,i+1) );
    end
end

datatype = get_tag(rdat,'datatype');
if ~isempty(datatype) & ~isempty(datatype{1}) & any(contains(datatype,'REACTIVITY:rho')) & any(contains(datatype,'READS')) ;
    warning('Trying to figure out errors based on datatype:READS lines!');
    rdat.reactivity_error = 0*rdat.reactivity;
    reads = get_tag(rdat,'reads');
    for i = find(contains(datatype,'REACTIVITY:rho')) % assume rho lines are followed by READS lines.    
        reads_mod   = [str2num(reads{i+1})-sum(rdat.reactivity(2:end,i+1)); rdat.reactivity(:,i+1)];
        reads_nomod = [str2num(reads{i+2})-sum(rdat.reactivity(2:end,i+2)); rdat.reactivity(:,i+2)];
        rdat.reactivity_error(:,i) = estimate_error_from_reads(rdat.reactivity(:,i), reads_mod,reads_nomod );
    end
end

if length(rdat.seqpos)>1 & (rdat.seqpos(2) < rdat.seqpos(1))
    warning('SEQPOS has wrong ordering; will try to correct SEQPOS, REACTIVITY, and REACTIVITY_ERROR')
    rdat.seqpos = rdat.seqpos(end:-1:1);
    rdat.reactivity = rdat.reactivity(end:-1:1,:);
    rdat.reactivity_error = rdat.reactivity_error(end:-1:1,:);
    sequence_seqpos = sequence_seqpos(end:-1:1);
end

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
  if length(s)>0 & strcmp(s(1),tabchar); s = s(2:end); end; % edge case 23SRR_H101_0001
  if ~isempty(strfind( s, tabchar ) ) 
    delim = tabchar;
  else
    delim = ' '; 
  end
end;

c = remove_empty_cells(strsplit(s,delim));

% rest = s;
% i = 1;
% c = {};
% while length(rest)
%     [t, rest] = strtok(rest, delim);
%     if length(t)>0
%         c{i} = t;
%         i = i + 1;
%     end
% end

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
elseif strfind( line, [tag,delim] ) 
  line = strrep(line, [tag,delim],'');
else
  line = strrep(line,tag,''); % edge case
end
line = strip(line);

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
  if isempty(tag) continue; end; % edge case 5SRRNA_1M7_0007
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
    warning(sprintf( 'seqpos %d is not inside sequence, given offset %d\n', seqpos(i), offset ));
    ok = 0;
    continue;
  end
  c2 = lower( sequence( m ) );
  if ( c1 ~= 'X' & c2 ~= 'X' & c1 ~= c2 & ~(c1=='u' & c2=='t')  & ~(c1=='t' & c2=='u')  )
    warning(sprintf( 'Mismatch at seqpos %d, between SEQPOS nucleotide %s and SEQUENCE nucleotide %s\n', seqpos(i), sequence_seqpos(i),  sequence( seqpos(i) - offset ) ));
    ok = 0;
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rdat = get_data_annotation( rdat, line )
if contains(line,sprintf('\t')) & contains(line,' sequence:'); line = strrep(line,' ',sprintf('\t')); end; % edge case ETERNA_R44_0001
%if contains(line,sprintf('\t')) & contains(line,' modifier:'); line = strrep(line,' ',sprintf('\t')); end; % edge case ETERNA_R44_0001
cols = strip(str2cell( line ));
idx = str2num( cols{1} );
if isempty( idx ) % weird edge case where str2cell fails
    firstcols = str2cell(cols{1});
    cols = [firstcols,cols(2:end)];
    idx = str2num( cols{1} );
end
anot = cols(2:end);
anot = strrep(anot,'ligpos:','lig_pos:'); % old MOHCA like RNAPZ5_MCA_0002.rdat
rdat.data_annotations{idx} = remove_empty_cells( anot );
% look for a 'sequence' tag.
for j = 1:length( anot )
    if ~isempty( strfind( anot{j}, 'sequence:' ) )
        [dummy, r ] = strtok( anot{j}, 'sequence:' );
        if ~isnan( str2double( dummy ) )
            assert( length( rdat.sequences ) >= str2num(dummy) );
            legacy_data_anot_sequences{idx} = rdat.sequences{str2num(dummy)}; % old-style format.
        else
            if strcmp(dummy(end),'.'); dummy = dummy(1:end-1); end; %two sequences each in ETERNA_R83_0000.rdat and ETERNA_R83_0003.rdat?
            if contains(dummy,'('); % edge case NEIL1_DMS_0001.rdat
                cols = strsplit(dummy,' '); %AGCCUGCCCUCUGAUCUCUGCCUGUUCCUCUGUCCCACAGGGGGCAAUGGCUACGGGUGAGAGAGCGGGGAGGAGGAC (A48U, C59G)
                dummy = cols{1};
            end
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
function rdat = check_for_legacy_datatype_reads(rdat, line_read);
n = line_read(1);
if n <= length(rdat.data_annotations) & any(strcmp(rdat.data_annotations{n},'datatype:READS'))
    if length(line_read(2:end))==length(rdat.seqpos) % weird shift in early Lucks lab cases.
        rdat.reactivity(1:length(rdat.seqpos), n) = [line_read(3:end),0];
    end
    total_reads = sum(line_read(2:end));
    rdat.data_annotations{n} = [rdat.data_annotations{n},sprintf('reads:%d',total_reads)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function estimated_error = estimate_error_from_reads(reactivity,reads_mod,reads_nomod );
[profile_mod,profile_error_mod] = get_reactivity_from_reads(reads_mod,sqrt(reads_mod));
[profile_nomod,profile_error_nomod] = get_reactivity_from_reads(reads_nomod,sqrt(reads_nomod));
signal = profile_mod - profile_nomod;
err = sqrt(profile_error_mod.^2 + profile_error_nomod.^2);
positive_pos = find(reactivity>0);
negative_pos = find(reactivity<=0);
scalefactor = mean( reactivity(positive_pos))/mean(signal(positive_pos));
%clf; plot([reactivity,scalefactor*signal]); pause;
estimated_error = scalefactor * err;
%[reads_mod(1), reads_nomod(1), max(estimated_error)]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [profile,profile_error] = get_reactivity_from_reads(reads,reads_error);
x = reads./cumsum(reads);
x_err = reads_error./cumsum(reads);
profile = x(2:end);
profile_error = x_err(2:end);

