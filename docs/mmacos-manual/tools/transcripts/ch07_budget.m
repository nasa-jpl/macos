% ch07_budget.m -- Chapter 7 (Error budget): templates/50_sensitivities/e5hex1_budget, the library's
% worked budget (Redding 2014: a covariance budget, with a Monte-Carlo cross-check), on the public
% 7-segment hex fixture e5hex1.in, post-WFC (the 7 segments + the SM controlled).  NOT the thread's
% telescope: e2e has no budget stage, and a budget spec for it (its group sigmas) is a design choice
% for Dave (OPEN_QUESTIONS.md); the chapter says how run_budget takes the thread's s4_sens.mat.
%
% The runner writes its products BESIDE ITSELF (the committed .xlsx/.csv/.mat/.png/report) and finds
% its deck by a relative path, so the driver copies e5hex1_budget/ and e5hex1/ into the scratch
% directory under the same layout and runs the copy, unchanged.  The library is never written.
%
% Pinned as BOUNDS, one line each (the report's alpha and WFE columns, and its summary table, hold
% the residuals the controller projects out -- numerical zeros of order 1e-8 .. 1e-13):
%   ch07_floor   < 1e-6   every alpha or WFE entry below 1e-6 prints as '< 1e-06'; the pin is the
%                         largest of them (the sigma column is an input and is never touched)
% Signals (printed): the pre-WFC and post-WFC totals, the Monte-Carlo total, the three largest rows.

root = getenv('MMACOS_ROOT');
work = getenv('MMMAN_WORK');
tpl  = fullfile(root, 'templates', '50_sensitivities');
dst  = fullfile(work, 'templates', '50_sensitivities');
mkdir(dst);
copyfile(fullfile(tpl, 'e5hex1_budget'), fullfile(dst, 'e5hex1_budget'));
copyfile(fullfile(tpl, 'e5hex1'), fullfile(dst, 'e5hex1'));
addpath(fullfile(root, 'src'));  addpath(fullfile(root, 'design', 'runners'));
here = fullfile(dst, 'e5hex1_budget');

mm_listing('ch07_spec', fullfile(here, 'budget_spec_e5hex1.m'));

t0 = tic;
run(fullfile(here, 'run_budget_e5hex1.m'));
fprintf('budget wall time %.0f s\n', toc(t0));
work = getenv('MMMAN_WORK');  here = fullfile(work, 'templates', '50_sensitivities', 'e5hex1_budget');

rep = fileread(fullfile(here, 'e5hex1_budget_report.txt'));
% ---- the projected-out residuals as bounds: alpha + WFE columns of the group table, the summary ----
L = splitlines(string(rep));  big = 0;  NUM = '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?';
grp = ~cellfun(@isempty, regexp(L, '\s(angle|length_si|wfs|act)\s+\S+\s+\S+\s+\S+\s*$'));
inSum = false;
for i = 1:numel(L)
    if startsWith(strtrim(L(i)), "summary"), inSum = true; continue; end
    if startsWith(strtrim(L(i)), "method:"), inSum = false; end
    if grp(i)
        [tok, ext] = regexp(char(L(i)), ['(' NUM ')\s+(' NUM ')\s*$'], 'tokens', 'tokenExtents', 'once');
        if ~isempty(ext), ext = num2cell(ext, 2); end   % alpha, WFE (the sigma before them is an input)
        cols = [1 2];
    elseif inSum                                        % every number on a summary row
        [s0, e0, tok] = regexp(char(L(i)), NUM, 'start', 'end', 'match');
        ext = num2cell([s0(:) e0(:)], 2);
        cols = 1:numel(tok);
    else
        continue
    end
    if isempty(tok), continue; end
    s = char(L(i));
    for c = fliplr(cols)
        v = str2double(tok{c});
        if isfinite(v) && v ~= 0 && abs(v) < 1e-6
            big = max(big, abs(v));  e = ext{c};
            s = [s(1:e(1)-1) '< 1e-06' s(e(2)+1:end)];
        end
    end
    L(i) = string(s);
end
mm_pin('ch07_floor', big, 'abs:1e-06', '< 1e-06');
rep = strjoin(L, newline);
mm_transcript('ch07_report', strtrim(rep));

mm_figure('ch07_budget', fullfile(here, 'e5hex1_budget.png'));

t = regexp(rep, 'pre-WFC TOTAL\s+(\S+)', 'tokens', 'once');      mm_pin('ch07_pre', str2double(t{1}), 'printed', t{1});
t = regexp(rep, '\nTOTAL \(RSS\)\s+(\S+)', 'tokens', 'once');      mm_pin('ch07_post', str2double(t{1}), 'printed', t{1});
t = regexp(rep, 'Monte-Carlo \((\d+) draws\) total: (\S+) nm', 'tokens', 'once');
mm_pin('ch07_mc_draws', str2double(t{1}), 'sig5', t{1});
mm_pin('ch07_mc', str2double(t{2}), 'printed', t{2});
for r = ["PM Seg Piston Actuator" "PM Seg Tilt Actuator" "SM Tilt Actuator" "WFS Sensing Error"]
    t = regexp(rep, r + "\s+\S+\s+\S+\s+\S+\s+(\S+)", 'tokens', 'once');
    mm_pin("ch07_" + lower(regexprep(r, '\s+', '_')), str2double(t{1}), 'printed', t{1});
end
