function [reactivities, seqpos ] = extract_reactivities( rdat,  outfile, idx, which_res );
% extract_reactivities( rdat,  outfile, idx, which_res );
%
%

if ischar( rdat ); rdat = read_rdat_file( rdat ); end;
if ~exist( 'idx','var') idx = 1; end;
if ~exist( 'which_res','var') which_res = []; end;

fid = fopen( outfile, 'w' );
for i = 1:length( rdat.seqpos );
  if length( which_res > 0 ) & isempty(find(rdat.seqpos(i)==which_res)); continue; end;
  fprintf( fid, '%d %9.5f\n', rdat.seqpos(i), rdat.reactivity(i,idx) );
end

fclose( fid );
fprintf( 'Created: %s\n', outfile );
