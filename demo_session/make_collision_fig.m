% MAKE_COLLISION_FIG  Slide-13 figure: the collimator-in-feed-beam collision,
% drawn in the part's own plane -- the failure the layout ray traces cannot
% show (it lives in the UNION over the field box, not in any single field).
%
% Engine truth: afocal4_union supplies the body hulls and floors; the
% per-field beam patches are the same repointed traces the gate uses
% (field_union_ recipe), intersected with the worst pair's body plane.
% Left panel: committed 343 mm (body inside the beam union, -79.9 mm).
% Right panel: cleared -10 deg (union moved off the part, +37.8 mm).
%
% Run:  MACOS_HOME=~/dev/macos/macos_f90 matlab -batch "run('make_collision_fig.m')"
% Writes figs/fig_collision.png.  Reads the record; modifies nothing.

run('/home/dcr/dev/MACOS_res_dev/mmacos/mmacos_setup.m');
A4 = '/home/dcr/dev/MACOS_res_dev/mmacos/challenges/afocal4';
addpath(A4); addpath(fullfile(A4,'clearing'));

P = afocal4_params;   F = P.Fsolve;
macos.init(256);

decks = { fullfile(A4,'afocal4_b2long_343mm.in'), ...
          fullfile(A4,'clearing','afocal4_clear_343mm.in') };
names = { 'committed 343 mm', 'cleared (tilt -10\circ, re-solved)' };

ACC = [0.12 0.31 0.47];  LIGHT = [0.55 0.68 0.82];

fig = figure('Position',[80 80 1280 520], 'Color','w');
tl = tiledlayout(fig, 1, 2, 'Padding','compact', 'TileSpacing','compact');

legk = [];  obst = [];
for d = 1:2
    K  = afocal4_union(decks{d}, 'fields', F, 'quiet', true);
    Kb = afocal4_union(decks{d}, 'fields', F, 'body_k', 1.0, ...
                       'body_pad', 0.0, 'quiet', true);
    if d == 1   % the worst pair on the COMMITTED deck names the story
        w = K.pair(K.worst);  legk = w.leg;  obst = w.obst;
        fprintf('worst pair: %s\n', w.name);
    end
    ip = find([K.pair.leg] == legk & [K.pair.obst] == obst, 1);
    ib = find([Kb.pair.leg] == legk & [Kb.pair.obst] == obst, 1);
    d_body = K.pair(ip).d_m;   d_bare = Kb.pair(ib).d_m;
    fprintf('%s: declared %+.2f mm, bare %+.2f mm\n', names{d}, ...
            d_body*1e3, d_bare*1e3);

    B  = K.body(obst);   Bb = Kb.body(obst);

    % per-field crossings of the leg with the body plane (field_union_ recipe)
    txt = fileread(decks{d});
    cd0 = g3(txt,'ChfRayDir');  cp0 = g3(txt,'ChfRayPos');  ap = g3(txt,'ApStop');
    stand = dot(ap - cp0, cd0);
    bx0 = asin(cd0(1));  by0 = asin(cd0(2));
    tmp = [tempname '.in'];
    UV = cell(1, size(F,1));   allUV = [];
    ihl = ceil(size(F,1)/2);   uvo = [];   % highlight field: its OWN glass
    for i = 1:size(F,1)
        bx = bx0 + F(i,1);   by = by0 + F(i,2);
        cdir = [sin(bx); sin(by); sqrt(max(0,1-sin(bx)^2-sin(by)^2))];
        cpos = ap - stand*cdir;
        s = regexprep(txt, '(ChfRayDir=\s*)[^\n]*', ['$1' v3(cdir)]);
        s = regexprep(s,   '(ChfRayPos=\s*)[^\n]*', ['$1' v3(cpos)]);
        s = regexprep(s,   '(yGrid=\s*)[^\n]*',     ['$1' v3([0;cos(by);-sin(by)])]);
        fid = fopen(tmp,'w');  fprintf(fid,'%s',s);  fclose(fid);
        macos.load_rx(tmp);
        macos.ray_hist('on');  t = macos.trace();  hi = macos.ray_hist(t.nRays);
        macos.ray_hist('off');
        off = size(hi.P,3) - K.nElt;
        m  = logical(hi.ok(:,legk+off)) & logical(hi.ok(:,legk+1+off));
        A  = squeeze(hi.P(:,m,legk+off));   Bp = squeeze(hi.P(:,m,legk+1+off));
        den = B.n.'*(Bp - A);
        tt  = (B.n.'*(B.c - A))./den;
        k   = isfinite(tt) & tt >= 0 & tt <= 1;
        X   = A(:,k) + (Bp(:,k) - A(:,k)).*tt(k);
        uv  = [B.u.'; B.v.']*(X - B.c);
        UV{i} = uv;   allUV = [allUV, uv]; %#ok<AGROW>
        if i == ihl   % this field's own footprint ON the part
            mo  = logical(hi.ok(:,obst+off));
            Qo  = squeeze(hi.P(:,mo,obst+off));
            uvo = [B.u.'; B.v.']*(Qo - B.c);
        end
    end
    delete(tmp);
    macos.load_rx(decks{d});

    ax = nexttile(tl);  hold(ax,'on');  axis(ax,'equal');
    % body: declared allowance (light) over bare lit glass (darker)
    pd = B.poly*1e3;   pb = Bb.poly*1e3;
    fill(ax, pd(1,[1:end 1]), pd(2,[1:end 1]), [0.88 0.88 0.88], ...
         'EdgeColor',[0.45 0.45 0.45], 'LineWidth',1.0);
    fill(ax, pb(1,[1:end 1]), pb(2,[1:end 1]), [0.70 0.70 0.70], ...
         'EdgeColor',[0.30 0.30 0.30], 'LineWidth',1.0);
    % per-field beam patches
    for i = 1:numel(UV)
        q = UV{i}*1e3;
        if size(q,2) >= 3
            h = convhull(q(1,:).', q(2,:).');
            fill(ax, q(1,h), q(2,h), ACC, 'FaceAlpha',0.07, ...
                 'EdgeColor',LIGHT, 'LineWidth',0.7);
        end
    end
    % the union of the feed beam over the field box
    qa = allUV*1e3;
    hu = convhull(qa(1,:).', qa(2,:).');
    plot(ax, qa(1,hu), qa(2,hu), '--', 'Color',ACC, 'LineWidth',2.0);
    % one field, told whole: its own glass (green) vs its feed beam (bold)
    if d == 1 && ~isempty(uvo)
        qo = uvo*1e3;   ho = convhull(qo(1,:).', qo(2,:).');
        fill(ax, qo(1,ho), qo(2,ho), [0.25 0.58 0.32], 'FaceAlpha',0.18, ...
             'EdgeColor',[0.12 0.42 0.18], 'LineWidth',1.6);
        text(ax, mean(qo(1,:)), min(qo(2,:))-9, ...
             'the centre field''s own glass', 'FontSize',9, ...
             'Color',[0.10 0.35 0.15], 'HorizontalAlignment','center');
        qc = UV{ihl}*1e3;   hc = convhull(qc(1,:).', qc(2,:).');
        plot(ax, qc(1,hc), qc(2,hc), '-', 'Color',ACC, 'LineWidth',1.8);
        text(ax, mean(qc(1,:)), max(qc(2,:))+9, ...
             'the same field''s feed beam', 'FontSize',9, 'Color',ACC, ...
             'HorizontalAlignment','center');
    end
    title(ax, sprintf('%s:  floor %+.1f mm (declared body), %+.1f mm bare glass', ...
          names{d}, d_body*1e3, d_bare*1e3), 'FontSize',11);
    xlabel(ax, 'in the part''s plane  [mm]');
    if d == 1, ylabel(ax, '[mm]'); end
    grid(ax,'on');  set(ax,'GridAlpha',0.15);
end
title(tl, ['the M2 \rightarrow field-mirror feed beam at the collimator: ' ...
    'nine fields (light), their union (dashed), and the part ' ...
    '(bare glass + declared allowance)'], 'FontSize',12, 'FontWeight','bold');

exportgraphics(fig, 'figs/fig_collision.png', 'Resolution', 150);
fprintf('wrote figs/fig_collision.png\n');

function x = g3(txt, key)
    tok = regexp(txt, [key '=\s*([^\n]*)'], 'tokens', 'once');
    x = sscanf(tok{1}, '%f');   x = x(1:3);
end
function s = v3(x)
    s = sprintf('%+.15e %+.15e %+.15e', x(1), x(2), x(3));
end
