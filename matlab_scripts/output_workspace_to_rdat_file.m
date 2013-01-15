function rdat = output_workspace_to_rdat_file( filename, name, sequence, offset, seqpos, reactivity, ...
					       mutpos, structure, ...
					       annotations, data_annotations, ...
					       reactivity_error, ...
					       trace_in, xsel, xsel_refine, comments );
%
% output_workspace_to_rdat_file( filename, name, sequence, offset, seqpos, reactivity, ...
%					       mutpos, structure, ...
%					       annotations, data_annotations, ...
%					       reactivity_error, ...
%					       trace, xsel, xsel_refine, comments );
%
% Copyright R. Das, P. Cordero, Stanford University, 2010,2011
%

% set defaults.
if ~exist( 'reactivity' );  fprintf( 'Must specify six variables: filename, name, sequence, offset, seqpos, reactivity'); end
if ~exist( 'mutpos' ); mutpos = []; end;
if ~exist( 'annotations' ); annotations = {}; end;
if ~exist( 'data_annotations' ); data_annotations = {}; end;
if ~exist( 'reactivity_error' ); reactivity_error = {}; end;
if ~exist( 'xsel' ); xsel = []; end;
if ~exist( 'xsel_refine' ); xsel_refine = []; end;
if ~exist( 'comments' ); comments = {}; end;
if ~exist( 'trace_in' ); trace_in = []; end;
trace = trace_in; % this is necessary because matlab has a function called trace, of course

rdat = fill_rdat( name, sequence, offset, seqpos, reactivity, mutpos, structure, ...
		  annotations, data_annotations, reactivity_error, trace, xsel, xsel_refine, comments );

if length( rdat.name ) == 0; 
  fprintf( 'PROBLEM filling rdat!\n' );
  return;
end

output_rdat_to_file( filename, rdat );
