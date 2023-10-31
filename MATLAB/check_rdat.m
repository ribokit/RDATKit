function ok = check_rdat( rdat )
% ok = check_rdat( rdat )
%
% INPUT
%  rdat = RDAT structure with chemical mapping data
%
% OUTPUT
%   ok = passed all checks.
%
% (C) R. Das, 2011-2013, 2018.

ok = 0;
if nargin==0; help( mfilename ); return; end;

% A bunch of consistency checks...
if isempty( rdat.name ); 
    warning( 'WARNING! Must give a name!'); 
    return; 
end;

if isempty( rdat.sequence ); 
    warning( 'WARNING! Must supply sequence!'); 
    return;
end;

if strfind( rdat.sequence, 'T' ); warning( '\nWARNING! Warning: you have a T instead of a U in the sequence!!' ); end;

if ( (min(rdat.seqpos) - rdat.offset) < 1 ); 
    warning( 'WARNING! Offset/seqpos does not look right -- at least one index is too low for sequence' ); 
    return;
end;

if ( (max(rdat.seqpos) - rdat.offset > length(rdat.sequence)) ); 
    warning( 'WARNING! Offset/seqpos does not look right -- at least one index %d is too high for sequence length %d', max(rdat.seqpos) - rdat.offset, length(rdat.sequence) ); 
    return;
end;

%if ~exist( 'rdat.reactivity' )
%  warning(sprintf( '\nWARNING! No REACTIVITY data -- assuming these data are in AREA_PEAK'));
%  rdat.reactivity = rdat.area_peak;
%end

if ( size( rdat.reactivity, 1 ) ~= length( rdat.seqpos ) );
  warning(sprintf( 'WARNING! Number of bands in reactivity [%d] does not match length of seqpos [%d]', size( rdat.reactivity, 1), length( rdat.seqpos ) )); 
  return;
end;

if ( ~isempty( rdat.data_annotations ) && length( rdat.data_annotations ) ~= size( rdat.reactivity, 2 ) );
  warning(sprintf( 'WARNING! Number of bands in data_annotations [%d] does not match number of lanes in reactivity [%d]', length( rdat.data_annotations), size( rdat.reactivity, 2 ) ));
  return;
end;

if ( ~isempty( rdat.xsel ) ) ;
  if ( size( rdat.reactivity, 1) ~= length( rdat.xsel ) );
    warning(sprintf( 'WARNING! Number of bands in xsel  [%d] does not match number of bands in reactivity [%d]', length( rdat.xsel), size( rdat.reactivity,1) ));
  end;
end;
   
if ( ~isempty( rdat.xsel_refine ) );
  if ( size( rdat.reactivity, 2) ~= size( rdat.xsel_refine, 2 ) );
    warning(sprintf( 'WARNING! Number of lanes in xsel_refine  [%d] does not match number of lanes in reactivity [%d]', size( rdat.xsel_refine, 2), size( rdat.reactivity,2) ));
  end;
  if ( size( rdat.reactivity, 1) ~= size( rdat.xsel_refine, 1 ) );
    warning(sprintf( 'WARNING! Number of bands in xsel_refine  [%d] does not match number of bands in reactivity [%d]', size( rdat.xsel_refine, 1), size( rdat.reactivity,1) ));
  end;
end;

OUTPUT_WARNING_ON_BLANK_REACTIVITY_ERROR = 0;
if ( size( rdat.reactivity_error,1) == 0 & OUTPUT_WARNING_ON_BLANK_REACTIVITY_ERROR)
  warning(sprintf( 'WARNING! Reactivity_error is blank! ' )); 
  return;
end;

if ( all( rdat.reactivity_error(:) == 0.0 ) & OUTPUT_WARNING_ON_BLANK_REACTIVITY_ERROR)
  %warning(sprintf( 'WARNING! Reactivity_error is all 0 ' )); 
  %return;
end;

if ( size( rdat.reactivity_error, 1 ) > 0 & size( rdat.reactivity, 1 ) ~= size( rdat.reactivity_error, 1 ) );
  warning(sprintf( 'WARNING! Number of bands in reactivity [%d] does not match length of reactivity_error [%d]', size( rdat.reactivity, 1), size( rdat.reactivity_error, 1 ) )); 
  return;
end;


if ~isempty( rdat.annotations ); 
    if ~check_annotations( rdat.annotations ); 
        return;
    end
end;

if ~isempty( rdat.data_annotations ) 
  for i = 1:length( rdat.data_annotations )
    check_annotations( rdat.data_annotations{i} );
  end
end

ok = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ok = check_annotations( annotations )
ok = 0;
for j = 1:length( annotations ) ;
  if ~check_annotation( annotations{j} ) ;
    warning(sprintf( 'WARNING! Unrecognized annotation: %s', annotations{j} ) );
    return;
  end;
end;
ok = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ok = check_annotation( annotation )
ok_annotations = {'chemical','modifier','experimentType','temperature','chemical','mutation','processing','ERROR','warning',...
    'EteRNA','Eterna',...
    'sequence','structure',...
    'MAPseq','sequenceSource','signal_to_noise','feature','lig_pos','offset',...
    'scaling','normalization','reads','reverse_transcriptase','time','name',...
    'ID','replicate','experiment'};

ok = 0;

t = strtok( annotation, ':' );
for i = 1:length( ok_annotations );
  if ( strcmp( t, ok_annotations{i} ) );
    ok = 1; 
    break;
  end;
end;

