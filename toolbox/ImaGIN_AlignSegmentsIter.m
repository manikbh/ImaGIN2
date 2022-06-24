function [ alignments, alignedStimulations ] = ImaGIN_AlignSegmentsIter( d, stimulations, segmentHWidth, templateHWidth, corrHWidth, do_plot, iterMax )

if size(d,1) ~= 1
    error( 'Invalid d' )
end

fprintf( 1, [ 'alignSegments... ' ] ) ;

alignments          = nan( length(stimulations), iterMax ) ;
alignedStimulations = stimulations ;
segments            = nan( length(alignedStimulations), 2*segmentHWidth+1 ) ;
for i=1:iterMax

%     fprintf( 1, [ 'alignSegments iter ' num2str(i) '...' '\n' ] ) ;
    
    for s=1:length(alignedStimulations)
        segments(s,:)   = d( alignedStimulations(s) + (-segmentHWidth:segmentHWidth) ) ;
    end
    
    [ alignment, ~ ]    = alignSegments( segments, templateHWidth, corrHWidth, do_plot ) ;
    alignments(:,i)     = alignment ;
    
    if isempty(find(alignment,1))
%         fprintf( 1, [ 'alignSegments iter ' num2str(i) ' => BREAK' '\n' ] ) ;
        break
    else
        % it happens that aligned segment fall outside full range
        % in this case, the segment is not aligned
        for s=1:length(alignedStimulations)
            if alignedStimulations(s)+alignment(s)+segmentHWidth > length(d)
                alignment(s) = 0 ;
            end
        end
        alignedStimulations  = alignedStimulations + alignment' ;
    end
    
%     if i == iterMax
%         error( [ '>> WARNING: alignment not achieved after iterMax = ' num2str(iterMax) ] ) ;
% %         fprintf( 1, [ '>> WARNING: alignment not achieved after iterMax = ' num2str(iterMax) '\n' ] ) ;
%     end
    
end

fprintf( 1, [ '(' num2str(i) ' iterations)' '\n' ] ) ;

end



function [ alignment, alignedSegments ] = alignSegments( segments, templateHWidth, corrHWidth, do_plot )
  
if ~ ( corrHWidth > templateHWidth )
    error( 'Invalid corrHWidth' )
end

segmentCount    = size(segments,1) ;
segmentWidth    = size(segments,2) ;
segmentCenter   = round(segmentWidth/2);

templateWidth   = 2*templateHWidth+1 ;
corrWidth       = 2*corrHWidth+1 ;
corrCount       = corrWidth - templateWidth + 1 ;
corr            = nan( segmentCount, corrCount ) ;

template        = mean( segments(:,segmentCenter+(-templateHWidth:templateHWidth)), 1 ) ;

for s=1:segmentCount
    for c=1:corrCount
        segStart    = segmentCenter - corrHWidth -1 + c ;
        segEnd      = segStart + templateWidth - 1 ;
        r           = corrcoef( template, segments(s,segStart:segEnd) ) ;
        corr(s,c)  	= r(2,1) ;
    end
end

[ ~, i ]    = max( corr,[],2) ;
alignment  	= segmentCenter - corrHWidth -1 + i + templateHWidth - segmentCenter ;

alignedSegments = applyAlignment( segments, alignment ) ;

if do_plot
    plotSegmentsAlignment( segments, alignedSegments, segmentCenter ) ;
end

end



function alignedSegments = applyAlignment( segments, alignment )

alignedSegments = nan( size(segments) ) ;
for s=1:size(segments,1)
    if alignment(s) > 0
        alignedSegments(s,1:end-alignment(s)) = segments(s,1+alignment(s):end) ;
    elseif alignment(s) < 0
        alignedSegments(s,1+abs(alignment(s)):end) = segments(s,1:end-abs(alignment(s))) ;
    else
        alignedSegments(s,:) = segments(s,:);
    end
end

end



function plotSegmentsAlignment( segments, alignedSegments, segmentCenter )

figure ;
ax1 = subplot(2,1,1); hold on ;
plot( 1:size(segments,2), segments ) ;
plot( 1:size(segments,2), mean(segments,1), 'color', 'r', 'linewidth', 4 ) ;
line( [ segmentCenter; segmentCenter ], get( gca, 'YLim' )', 'color', 'm', 'linewidth', 2 ) ;
xlim( [ 1 size(segments,2) ] ) ;
title( 'templates' )
ax2 = subplot(2,1,2); hold on ;
plot( 1:size(segments,2), alignedSegments ) ;
plot( 1:size(segments,2), mean(alignedSegments,1), 'color', 'r', 'linewidth', 4 ) ;
line( [ segmentCenter; segmentCenter ], get( gca, 'YLim' )', 'color', 'm', 'linewidth', 2 ) ;
xlim( [ 1 size(segments,2) ] ) ;
title( 'templates aligned' )
linkaxes( [ ax1, ax2 ], 'xy' )

end


