clear; close all; clc;

%% Settings
figfile = 'fit_byom_metabolism_template_1.fig';
outfile = 'fit_byom_metabolism_template_1_modelFits_CIs.csv';

%% Open figure invisibly
fig = openfig(figfile,'invisible');
set(fig,'Visible','off');

%% --- Find the 3 real data axes (ignore legend/annotation axes) ---
axAll = findall(fig,'Type','axes');

% keep only axes that contain at least one line with numeric X/Y data
axKeep = [];
for i = 1:numel(axAll)
    ln = findall(axAll(i),'Type','line');
    ok = false;
    for k = 1:numel(ln)
        xd = get(ln(k),'XData'); yd = get(ln(k),'YData');
        if isnumeric(xd) && isnumeric(yd) && numel(xd) > 1 && numel(yd) > 1
            ok = true; break;
        end
    end
    if ok
        axKeep = [axKeep; axAll(i)]; %#ok<AGROW>
    end
end

% if more than 3 axes survive, take the 3 with most line objects (typical)
if numel(axKeep) > 3
    nLines = arrayfun(@(a) numel(findall(a,'Type','line')), axKeep);
    [~, idx] = sort(nLines,'descend');
    axKeep = axKeep(idx(1:3));
end

if numel(axKeep) ~= 3
    error('Expected 3 data axes, found %d. The figure structure is unusual.', numel(axKeep));
end

% stable ordering by axes position: top-to-bottom, left-to-right
pos = vertcat(axKeep.Position);
[~, ord] = sortrows([ -pos(:,2), pos(:,1) ]);
axKeep = axKeep(ord);

%% --- Extract fit (solid) + CI (dotted) per axes ---
fits = cell(3,1);
ciHi = cell(3,1);
ciLo = cell(3,1);

isNoMarker = @(h) isempty(get(h,'Marker')) || strcmpi(get(h,'Marker'),'none');

for a = 1:3
    ax = axKeep(a);
    lnAll = findall(ax,'Type','line');

    % keep only lines with usable numeric data
    ln = [];
    for k = 1:numel(lnAll)
        xd = get(lnAll(k),'XData'); yd = get(lnAll(k),'YData');
        if isnumeric(xd) && isnumeric(yd) && numel(xd) > 1 && numel(yd) > 1
            ln = [ln; lnAll(k)]; %#ok<AGROW>
        end
    end

    % candidates by style:
    %   Fit: solid line '-' and no marker
    %   CI : dotted line ':' and no marker
    fitCand = [];
    ciCand  = [];
    for k = 1:numel(ln)
        ls = get(ln(k),'LineStyle');
        if isNoMarker(ln(k))
            if strcmp(ls,'-')
                fitCand = [fitCand; ln(k)]; %#ok<AGROW>
            elseif strcmp(ls,':')
                ciCand  = [ciCand; ln(k)]; %#ok<AGROW>
            end
        end
    end

    if isempty(fitCand)
        error('No solid (fit) line detected in axes %d.', a);
    end
    if numel(ciCand) < 2
        error('Fewer than 2 dotted (CI) lines detected in axes %d.', a);
    end

    % choose fit line as the one with the longest YData (most common)
    fitLens = arrayfun(@(h) numel(get(h,'YData')), fitCand);
    [~, iFit] = max(fitLens);
    fitLine = fitCand(iFit);

    x = get(fitLine,'XData'); x = x(:);
    y = get(fitLine,'YData'); y = y(:);

    % choose two CI lines that match the fit length (robust against extra dotted lines)
    ciLens = arrayfun(@(h) numel(get(h,'YData')), ciCand);
    match = find(ciLens == numel(y));
    if numel(match) < 2
        error('Could not find two CI lines matching fit length in axes %d.', a);
    end
    ci1 = ciCand(match(1)); ci2 = ciCand(match(2));
    y1 = get(ci1,'YData'); y1 = y1(:);
    y2 = get(ci2,'YData'); y2 = y2(:);

    % assign hi/lo by mean value
    if mean(y1,'omitnan') >= mean(y2,'omitnan')
        yHi = y1; yLo = y2;
    else
        yHi = y2; yLo = y1;
    end

    fits{a} = [x, y];
    ciHi{a} = yHi;
    ciLo{a} = yLo;
end

%% --- Export ---
% if all x grids match exactly, export one time column; otherwise export time_1..time_3
sameX = true;
for a = 2:3
    sameX = sameX && numel(fits{a}(:,1))==numel(fits{1}(:,1)) && all(abs(fits{a}(:,1)-fits{1}(:,1)) < 1e-12);
end

if sameX
    M = [fits{1}(:,1), fits{1}(:,2), ciHi{1}, ciLo{1}, ...
                     fits{2}(:,2), ciHi{2}, ciLo{2}, ...
                     fits{3}(:,2), ciHi{3}, ciLo{3}];
    header = {'time','fit_1','ci_hi_1','ci_lo_1','fit_2','ci_hi_2','ci_lo_2','fit_3','ci_hi_3','ci_lo_3'};
else
    M = [fits{1}(:,1), fits{1}(:,2), ciHi{1}, ciLo{1}, ...
         fits{2}(:,1), fits{2}(:,2), ciHi{2}, ciLo{2}, ...
         fits{3}(:,1), fits{3}(:,2), ciHi{3}, ciLo{3}];
    header = {'time_1','fit_1','ci_hi_1','ci_lo_1','time_2','fit_2','ci_hi_2','ci_lo_2','time_3','fit_3','ci_hi_3','ci_lo_3'};
end

fid = fopen(outfile,'w');
fprintf(fid,'%s',header{1});
for i = 2:numel(header)
    fprintf(fid,',%s',header{i});
end
fprintf(fid,'\n');
fclose(fid);

writematrix(M, outfile);   % <-- FIXED

disp(['Wrote: ', outfile]);
