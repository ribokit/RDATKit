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
if length( c ) == 0; return; end;
if strcmp(c{1},'mutation' )
    tag = strjoin( c(2:end), ':' );
else
    tag = annotation_tag;
end
if length(tag)>3 & strcmp(tag(1:3),'Lib') % edge case TRP4P6_SHP_0002
    assert(strcmp(tag(5),'-'));
    tag = tag(6:end);
end

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
    if ~contains( tag, 'WT' )
        warning(sprintf('Could not find mutation position in mutation annotation: %s\n', tag ));
    end
    mutpos = NaN;
    return;
end
if isempty( start_seq );
    warning(sprintf('Could not find starting nucleotide in mutation annotation: %s\n', tag ));
end
if isempty( mut_seq );
    warning(sprintf('Could not find mutation nucleotide in mutation annotation: %s\n', tag ));
end
mutpos = str2num( mut_num ) - rdat.offset;
if ( mutpos < 1 | ~strcmp( rdat.sequence( mutpos ), start_seq ) )
    mutpos_alt = str2num( mut_num ); % perhaps specified without offset...
    if ( strcmp( rdat.sequence( mutpos_alt ), start_seq ) )
        % Specified mutpos without taking into account offset...
        mutpos = mutpos_alt;
    else
        warning(sprintf('Mismatch between mutation nucleotides: %s (mutation annotation) vs. %s. (sequence) in mutation %s.\n', start_seq, rdat.sequence(mutpos), tag ));
    end
end