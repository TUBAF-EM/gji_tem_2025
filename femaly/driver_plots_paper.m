%% Modellantwort vs. empymod
% KI-generiert:

% 1. Dateien laden
fig1 = hgload('results/timings_dBdt.fig', struct('Visible','off'));
fig2 = hgload('results/timings_E.fig', struct('Visible','off'));

% 2. Nur echte Plot-Achsen finden (schließt Legenden-Handles aus)
ax1_all = flipud(findobj(fig1, 'Type', 'axes', '-not', 'Tag', 'legend', '-not', 'Tag', 'Colorbar'));
ax2_all = flipud(findobj(fig2, 'Type', 'axes', '-not', 'Tag', 'legend', '-not', 'Tag', 'Colorbar'));

% Hilfsfunktion: Kopiert beides, gibt aber nur das Achsen-Handle (1) an subplot zurück
copyAx = @(sourceAx, targetFig) copyobj([sourceAx, findobj(sourceAx.Parent, ...
           'Type', 'legend', 'Axes', sourceAx)], targetFig);

% 3. Neue Ziel-Abbildung
final_fig = figure(3);
clf;

% 4. Subplots kopieren
% WICHTIG: Das (1) am Ende extrahiert das Achsen-Handle aus dem kopierten [Achse, Legende] Array

% Fig 1
res1 = copyAx(ax1_all(1), final_fig);
set(res1(1).Children(1:2), 'Visible', 'off', 'HandleVisibility', 'off');
tmp1 = subplot(2, 3, 1, res1(1));
set(tmp1, 'XTick', [1e-7, 1e-5, 1e-3]);

res4 = copyAx(ax1_all(2), final_fig);
tmp4 = subplot(2, 3, 4, res4(1));
set(tmp4, 'YLim', [-10 10], 'XTick', [1e-7, 1e-5, 1e-3]);

res3 = copyAx(ax1_all(5), final_fig);
set(res3(1).Children(1:2), 'Visible', 'off', 'HandleVisibility', 'off');
tmp3 = subplot(2, 3, 3, res3(1));
set(tmp3, 'XTick', [1e-7, 1e-5, 1e-3]);

res6 = copyAx(ax1_all(6), final_fig);
tmp6 = subplot(2, 3, 6, res6(1));
set(tmp6, 'YLim', [-10 10], 'XTick', [1e-7, 1e-5, 1e-3]);

% Fig 2
res2 = copyAx(ax2_all(3), final_fig);
set(res2(1).Children(1), 'Visible', 'off', 'HandleVisibility', 'off');
tmp2 = subplot(2, 3, 2, res2(1));
set(tmp2, 'XTick', [1e-7, 1e-5, 1e-3]);
tmp2.Title.String{1}='x_s = [0, 0, 0], x_0 = [20, 0, 0]';

res5 = copyAx(ax2_all(4), final_fig);
tmp5 = subplot(2, 3, 5, res5(1));
set(tmp5, 'YLim', [-10 10], 'XTick', [1e-7, 1e-5, 1e-3]);
tmp5.Title.String = 'Relative error (%)';

% Legenden verschieben
set(findall(gcf, 'type', 'legend'), 'Location', 'northeast');

exportgraphics(gcf, 'timings_E-dBdt.pdf', 'ContentType', 'vector', 'BackgroundColor','none');

% 5. Aufräumen
close(fig1);
close(fig2);

%% Zeitdiskretisierung

add_t = load('results/timings_times.mat');
% Die Indizes bestimmen, an denen eine Änderung zum Vorgänger stattfindet
indices = [1, find(diff(add_t.dt) ~= 0) + 1];
num_t = diff([indices, length(add_t.t)]);
num_t = cell2mat(arrayfun(@(x) {1:x}, num_t));
num_t = [num_t, num_t(end)];

figure(1);
clf;
semilogx(add_t.t, num_t, '.');
xline(add_t.t_obs_);
legend('t^i', 't^g');
title('Time discretization');
xlabel('Time (s)');
ylabel('n_{t}');

exportgraphics(gcf, 'timings_time1D.pdf', 'ContentType', 'vector', 'BackgroundColor','none');

%% Timings / Speedup

cpu = [1,2,4,8,16,20];
cpu_times = [1850, 1310, 614, 383, 259, 215];
cpu_speedup = cpu_times(1) ./ cpu_times(2:end);

figure(1);
clf;
title('Calculation times S: n_{g} = 20, n_{d} = 6, n_{FBS} = 6050');
yyaxis left
plot(cpu, cpu_times, '.-');
ylabel('t_{solve}');
yyaxis right
plot(cpu(2:end), cpu_speedup, '.-');
ylim([0, 20]);
ylabel('speedup');
xlabel('num_{worker}');

exportgraphics(gcf, 'timings_S.pdf', 'ContentType', 'vector', 'BackgroundColor','none');

%%

figTt = hgload('results/timings_Taylor.fig');
figTt.Children(2).Title.String = 'Taylor test convergence';
figTt.Children(2).YLabel.String = 'Normalized error function';

exportgraphics(gcf, 'timings_Taylor.pdf', 'ContentType', 'vector', 'BackgroundColor','none');
