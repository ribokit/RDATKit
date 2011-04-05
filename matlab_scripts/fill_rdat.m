function rdat = fill_rdat( name, sequence, offset, seqpos, area_peak, mutpos, structure, ...
			   annotations, data_annotations, area_peak_error, trace_in, xsel, xsel_refine, comments );
% rdat = fill_rdat( name, sequence, offset, seqpos, area_peak, mutpos, structure, ...
%			   annotations, data_annotations, area_peak_error, trace, xsel, xsel_refine, comments );
%
% Copyright R. Das, P. Cordero, Stanford University, 2010,2011
%

rdat = RDATFile;
rdat.name = ''; % signal that we're not filled yet.

if ~exist( 'area_peak' );  fprintf( 'Must specify six variables: filename, name, sequence, offset, seqpos, area_peak'); end
if ~exist( 'mutpos' ); mutpos = []; end;
if ~exist( 'structure' ); structure=''; end;
if ~exist( 'annotations' ); annotations = {}; end;
if ~exist( 'data_annotations' ); data_annotations = {}; end;
if ~exist( 'area_peak_error' ); area_peak_error = {}; end;
if ~exist( 'trace_in' ); trace_in = []; end;
if ~exist( 'xsel' ); xsel = []; end;
if ~exist( 'xsel_refine' ); xsel_refine = []; end;
if ~exist( 'comments' ); comments = {}; end;

if length( xsel_refine ) > 0 & ( size( area_peak, 2) ~= size( xsel_refine, 2 ) )
  xsel_refine = xsel_refine';
end


% Current version of fill_rdat script... work in progress!
rdat.version = 0.22;
rdat.comments = comments;
rdat.name = name;
rdat.sequence = sequence;
rdat.structure = structure;
rdat.offset = offset;
rdat.seqpos = seqpos;
rdat.mutpos = mutpos;
rdat.annotations = annotations;
rdat.data_annotations = data_annotations;
rdat.area_peak = area_peak;
rdat.area_peak_error = area_peak_error;
rdat.xsel = xsel;
rdat.xsel_refine = xsel_refine;
rdat.trace = trace_in;

check_rdat( rdat );
