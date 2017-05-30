function [ mutpos, mut_seq ] = get_mutation_info_from_tag( annotation_tag, rdat );
% [mutpos,mut_seq] = get_mutation_info_from_tag( annotation_tag );
%
% annotation_tag = string like "mutation:A25C"
%
% (C) R. Das, Stanford University, 2017.

c = strsplit( annotation_tag,':' );

start_seq = '';
mut_seq = '';
mut_num = '';
tag = strjoin( c(2:end), ':' );

q = 1;
while ( q <= length(tag)  & isempty( str2num( tag(q) ) ) & tag(q)~='(' & tag(q)~='-' )
    start_seq = [start_seq, tag(q) ];
    q = q+1;
end
if q <= length(tag) & tag(q) == '('; q = q+1; end;
while ( q <= length( tag ) &  (~isempty( str2num( tag(q) ) ) | tag(q)==':' | tag(q)=='-' )  )
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
    mutpos = NaN;
    return;
end
if isempty( start_seq );
    fprintf( 'WARNING! Could not find starting nucleotide in mutation annotation: %s\n', tag );
end
if isempty( mut_seq );
    fprintf( 'WARNING! Could not find mutation nucleotide in mutation annotation: %s\n', tag );
end

mutpos = str2num( mut_num ) - rdat.offset;
if ( mutpos < 1 | ~strcmp( rdat.sequence( mutpos ), start_seq ) )
    mutpos = str2num( mut_num ); % perhaps specified without offset...
    fprintf( 'WARNING! Mismatch between mutation nucleotides: %s (mutation annotation) vs. %s. (sequence) in mutation %s.\n', rdat.sequence(mutpos), start_seq, tag );
    if ( strcmp( rdat.sequence( mutpos ), start_seq ) )
        fprintf( 'OK, specified mutpos without taking into account offset...\n' );
    else
        fprintf( 'WARNING! Mismatch between mutation nucleotides: %s (mutation annotation) vs. %s. (sequence)\n', rdat.sequence(mutpos), start_seq );
    end
end