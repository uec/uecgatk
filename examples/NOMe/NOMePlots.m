context='WCG';
filePresFull = { ...
    'ALIGN_2000_1000_mapq30_minct3_mincontext0.90_ctcf'
};
filePresStratified = { ...
    'ALIGN_2000_1000_mapq30_minct3_mincontext0.90_ctcfPairs-featAligner-WCG_GCH'
};
names = { ...
    'ctcf'
};
matches = labelNames(names,'ctcf'); % just to start list.
%matches = matches | labelNames(names,'promoter');

matches = 1:length(names);

if (1>0)
    n = length(filePres);
    for centralfrac = 1; %[0.5 1]
        for chrxonly = [0]; %0; %[0 1] % 0:1
            for i = [find(matches)'];
                DOFREQS = []; %1;
                SMOOTH=5;
                name = names{i};
                fprintf('Working on (centralfrac=%.2f, chrx=%d) %s (%s)\n',centralfrac, chrxonly, name,filePre);
                
                % First the complete WCG pattern
                filePre = filePresFull{i};
                %[a,b,c] = processGnomeAlignenmentsFull(filePre,sprintf('fullAlignments-%s',name),[],SMOOTH,0,'C');
                [hcg, gch, hcgMean, gchMean] = processGnomeAlignenments(filePre, sprintf('fullAlignments-%s',name), [], 1, [], [], 1, SMOOTH);

                
                % Do the stratified NOMe combos
                filePre = filePresStratified{i};
                [meanMat, totals, namesOut] = processGnomeAlignenmentsCombos(context, filePre, sprintf('comboAlignments-%s',name),SMOOTH,centralfrac,chrxonly);
            end
        end
    end
end
