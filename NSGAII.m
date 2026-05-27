\begin{lstlisting}[language=Matlab,
    breaklines=true,
    columns=fullflexible,
    basicstyle=\ttfamily\small,
    commentstyle=\color{myblue}\ttfamily\small,
    frame=single,
    showspaces=false,
    showstringspaces=false
]
function perforatedElectrodeNSGA2()
% NSGA-II optimisation of perforated Ag/AgCl defibrillation electrode.
% Decision variables: hole radius [m], holes per side, r_perf [m]
% Objectives: minimise K (non-uniformity) and peak |J| in tissue
% Ref: Krasteva & Papazov (2002), BioMedical Engineering OnLine.

close all; clc;

fprintf('==========================================================\n');
fprintf(' NSGA-II OPTIMISATION: PERFORATED Ag/AgCl ELECTRODE\n');
fprintf('==========================================================\n\n');

%% 1. Parameters and geometry
p = defaultParameters();
[p.dl, p.externalEdges, p.elec1Edge, p.elec2Edge] = buildGeometry(p);

%% 2. Variable bounds
lb    = [0.0002,  2, 0.005];
ub    = [0.0025, 30, 0.045];
nVars = 3;

fprintf('Bounds:\n');
fprintf('  hole radius:    [%.1f, %.1f] mm\n', lb(1)*1e3, ub(1)*1e3);
fprintf('  holes per side: [%d, %d]\n',         lb(2),     ub(2));
fprintf('  r_perf:         [%.1f, %.1f] mm\n\n', lb(3)*1e3, ub(3)*1e3);

%% 3. NSGA-II setup
popSize = 60;
maxGen  = 40;
opts = optimoptions('gamultiobj', ...
    'PopulationSize', popSize, 'MaxGenerations', maxGen, ...
    'CrossoverFraction', 0.8, 'ParetoFraction', 0.35, ...
    'FunctionTolerance', 1e-4, 'Display', 'iter', 'UseParallel', false);
intcon = 2;

%% 4. Run optimisation
allX = [];  allF = [];

    function f = evaluateAndStore(x)
        f = evaluateDesign(x, p);
        allX = [allX; x];
        allF = [allF; f];
    end

fprintf('Starting NSGA-II (pop=%d, gen=%d) ...\n', popSize, maxGen);
fprintf('----------------------------------------------------------\n');
tStart = tic;
[xPareto, fPareto] = gamultiobj(@evaluateAndStore, nVars, ...
    [], [], [], [], lb, ub, [], intcon, opts);
tTotal = toc(tStart);
fprintf('----------------------------------------------------------\n');
fprintf('Done in %.1f min | %d Pareto, %d total evals\n\n', ...
    tTotal/60, size(xPareto,1), size(allX,1));

%% 5. Post-process
[~, si] = sort(fPareto(:,1));
xPareto = xPareto(si,:);
fPareto = fPareto(si,:);

validMask = allF(:,1) < 50 & allF(:,2) < 50;
allXv = allX(validMask,:);
allFv = allF(validMask,:);

fprintf('PARETO FRONT:\n');
fprintf('  %-4s  %-8s  %-6s  %-10s | %-8s  %-10s\n', ...
    '#','r(mm)','n/s','r_perf(mm)','K','peakJ');
for ii = 1:size(xPareto,1)
    fprintf('  %-4d  %-8.2f  %-6d  %-10.1f | %-8.2f  %-10.4f\n', ...
        ii, xPareto(ii,1)*1e3, round(xPareto(ii,2)), ...
        xPareto(ii,3)*1e3, fPareto(ii,1), fPareto(ii,2));
end

[bestKIdx, kneeIdx, bestJIdx] = identifyKeySolutions(fPareto);
labels = {'Best K','Knee','Best peakJ'};
for kk = 1:3
    ii = [bestKIdx, kneeIdx, bestJIdx]; ii = ii(kk);
    fprintf('  %-10s: K=%.2f peakJ=%.2f (r=%.2fmm n=%d rp=%.1fmm)\n', ...
        labels{kk}, fPareto(ii,1), fPareto(ii,2), ...
        xPareto(ii,1)*1e3, round(xPareto(ii,2)), xPareto(ii,3)*1e3);
end

%% 6. Solid baseline
fprintf('\nSolid baseline: ');
[~, ~, Ks, ~, pJs] = solveElectrode(0, 0, 0.025, p, false);
fprintf('K=%.2f, peakJ=%.2f\n\n', Ks, pJs);

%% 7. High-res re-solve
fprintf('Re-solving at fine mesh ...\n');
xBestJ = xPareto(bestJIdx,:);
[resSolid, rzSolid, KsHR,  ~, pJsHR]  = solveElectrode(0, 0, 0.025, p, true);
[resBestJ, rzBestJ, KbjHR, ~, pJbjHR]  = ...
    solveElectrode(xBestJ(1), round(xBestJ(2)), xBestJ(3), p, true);
fprintf('  Solid:  K=%.2f, peakJ=%.2f\n', KsHR, pJsHR);
fprintf('  BestJ:  K=%.2f, peakJ=%.2f\n\n', KbjHR, pJbjHR);

%% 8. Line profiles (0.5 mm below interface)
evalZ = p.d + 0.0005;
[distSolid, JlineSolid] = extractLineProfile(rzSolid, p, evalZ);
[distPerf,  JlinePerf]  = extractLineProfile(rzBestJ, p, evalZ);

%% 9. Figures
cfg = plotDefaults();
plotParetoFront(allXv, allFv, xPareto, fPareto, bestKIdx, kneeIdx, bestJIdx, Ks, pJs, cfg);
plotCurrentMaps(resSolid, resBestJ, p, KsHR, KbjHR, cfg);
plotZoomedMaps(rzSolid, rzBestJ, p, cfg);
plotLineProfiles(distSolid, JlineSolid, distPerf, JlinePerf, p, xBestJ, cfg);
plotSensitivity(allXv, allFv, xPareto, fPareto, Ks, pJs, cfg);

%% 10. Save
save('nsga2_perforated_results.mat', ...
    'xPareto','fPareto','allXv','allFv','p', ...
    'Ks','pJs','KsHR','pJsHR','KbjHR','pJbjHR', ...
    'xBestJ','bestKIdx','bestJIdx','kneeIdx');
fprintf('Results saved. All figures generated.\n');
end


%% ---- defaultParameters ----
function p = defaultParameters()
    p.domainWidth  = 0.50;       % m
    p.domainHeight = 0.25;       % m
    p.heartRadius  = 0.045;      % m
    p.heartX = 0.25;  p.heartZ = 0.095;
    p.probeX = 0.25;  p.probeZ = 0.05;
    p.sigma_tissue = 0.2;        % S/m
    p.sigma_heart  = 0.5;
    p.Vapp     = 10;             % V
    p.padWidth = 0.10;           % m
    p.pad1X = 0.25;  p.pad2X = 0.25;
    p.d = 0.001;                 % electrode thickness
    p.rho_metal   = 5;           % surface resistivity
    p.z_metal     = p.rho_metal * p.d;
    p.sigma_metal = p.d / p.z_metal;
    p.sigma_hole  = 1e-8;        % insulating
    p.pad1L = p.pad1X - p.padWidth/2;
    p.pad1R = p.pad1X + p.padWidth/2;
    p.pad2L = p.pad2X - p.padWidth/2;
    p.pad2R = p.pad2X + p.padWidth/2;
    p.R_pad = p.padWidth / 2;
end


%% ---- buildGeometry ----
function [dl, externalEdges, elec1Edge, elec2Edge] = buildGeometry(p)
% Rectangle (torso) + circle (heart). Bottom and top edges split at
% electrode pad boundaries so BCs can be applied per-edge.
    vertices = [0,0; p.pad1L,0; p.pad1R,0; p.domainWidth,0;
                p.domainWidth,p.domainHeight; p.pad2R,p.domainHeight;
                p.pad2L,p.domainHeight; 0,p.domainHeight];
    nv = size(vertices,1);
    pgon = [2; nv; vertices(:,1); vertices(:,2)];
    circ = [1; p.heartX; p.heartZ; p.heartRadius; zeros(6,1)];
    maxLen = max(length(pgon),length(circ));
    pgon(end+1:maxLen) = 0;  circ(end+1:maxLen) = 0;

    [dl,~] = decsg([pgon,circ], 'P1+C1', char('P1','C1')');

    numEdges = size(dl,2);
    externalEdges = [];
    for eID = 1:numEdges
        if dl(6,eID)==0 || dl(7,eID)==0
            externalEdges(end+1) = eID; %#ok<AGROW>
        end
    end

    elec1Edge = 2;  elec2Edge = 6;
    tol = 1e-6;
    for eID = externalEdges
        mx = (dl(2,eID)+dl(3,eID))/2;
        my = (dl(4,eID)+dl(5,eID))/2;
        if abs(my)<tol && mx>p.pad1L-tol && mx<p.pad1R+tol
            elec1Edge = eID;
        end
        if abs(my-p.domainHeight)<tol && mx>p.pad2L-tol && mx<p.pad2R+tol
            elec2Edge = eID;
        end
    end
    fprintf('Geometry: %d edges, elec1=%d, elec2=%d\n\n', numEdges, elec1Edge, elec2Edge);
end


%% ---- solveElectrode ----
function [res, resZoom, K, JAtProbe, peakJtissue] = ...
        solveElectrode(holeR, nH, rPerf, p, highRes)
% Solve -div(sigma*grad(V))=0 for a given perforation pattern.
% holeR=0 or nH=0 gives solid electrode.

    % Compute hole positions
    if holeR <= 0 || nH < 1
        useHoles = false;  holeLArr = [];  holeRArr = [];
        rPerfUsed = rPerf;
    else
        useHoles = true;
        rPerfUsed = max(0.002, min(rPerf, p.R_pad - 0.002));
        L = p.R_pad - rPerfUsed;
        if L < 0.002
            useHoles = false;  holeLArr = [];  holeRArr = [];
        else
            spacing = L / (nH + 1);
            actualR = min(holeR, spacing * 0.40);  % 40% cap prevents overlap
            centresR = p.pad1X + rPerfUsed + spacing*(1:nH)';
            centresL = p.pad1X - rPerfUsed - spacing*(1:nH)';
            allCentres = [centresL; centresR];
            holeLArr = allCentres - actualR;
            holeRArr = allCentres + actualR;
        end
    end

    function tf = inHole(xi)
        if ~useHoles || abs(xi - p.pad1X) < rPerfUsed
            tf = false;
        else
            tf = any((xi >= holeLArr) & (xi <= holeRArr));
        end
    end

    % PDE model
    mdl = createpde();
    geometryFromEdges(mdl, p.dl);
    if highRes
        generateMesh(mdl, 'Hmax',0.0012, 'Hmin',0.0003, 'GeometricOrder','linear');
    else
        generateMesh(mdl, 'Hmax',0.002, 'Hmin',0.0005, 'GeometricOrder','linear');
    end

    % Conductivity: tissue/heart/metal/hole depending on location
    function sig = condCoeff(region, ~)
        x = region.x;  z = region.y;
        sig = p.sigma_tissue * ones(size(x));
        for ii = 1:numel(x)
            xi = x(ii);  zi = z(ii);
            if (xi-p.heartX)^2 + (zi-p.heartZ)^2 <= p.heartRadius^2
                sig(ii) = p.sigma_heart;
            elseif zi <= p.d && xi >= p.pad1L && xi <= p.pad1R
                if inHole(xi), sig(ii) = p.sigma_hole;
                else,          sig(ii) = p.sigma_metal; end
            elseif zi >= p.domainHeight-p.d && xi >= p.pad2L && xi <= p.pad2R
                if inHole(xi), sig(ii) = p.sigma_hole;
                else,          sig(ii) = p.sigma_metal; end
            end
        end
    end

    % BCs: V=0 ground, V=Vapp active, dV/dn=0 elsewhere
    specifyCoefficients(mdl, 'm',0, 'd',0, 'c',@condCoeff, 'a',0, 'f',0);
    applyBoundaryCondition(mdl, 'dirichlet', 'Edge', p.elec1Edge, 'u', 0);
    applyBoundaryCondition(mdl, 'dirichlet', 'Edge', p.elec2Edge, 'u', p.Vapp);
    for eID = setdiff(p.externalEdges, [p.elec1Edge, p.elec2Edge])
        applyBoundaryCondition(mdl, 'neumann', 'Edge', eID, 'g', 0, 'q', 0);
    end

    % Solve
    result = solvepde(mdl);
    nodeX = mdl.Mesh.Nodes(1,:)';
    nodeZ = mdl.Mesh.Nodes(2,:)';
    [gradX, gradZ] = evaluateGradient(result, nodeX, nodeZ);

    % Nodal conductivities for J = -sigma*grad(V)
    sigN = p.sigma_tissue * ones(size(nodeX));
    for ii = 1:numel(nodeX)
        xi = nodeX(ii);  zi = nodeZ(ii);
        if (xi-p.heartX)^2+(zi-p.heartZ)^2 <= p.heartRadius^2
            sigN(ii) = p.sigma_heart;
        elseif zi <= p.d && xi >= p.pad1L && xi <= p.pad1R
            if inHole(xi), sigN(ii) = p.sigma_hole;
            else,          sigN(ii) = p.sigma_metal; end
        elseif zi >= p.domainHeight-p.d && xi >= p.pad2L && xi <= p.pad2R
            if inHole(xi), sigN(ii) = p.sigma_hole;
            else,          sigN(ii) = p.sigma_metal; end
        end
    end
    Jmag = sqrt((-sigN.*gradX).^2 + (-sigN.*gradZ).^2);

    % Interpolate onto regular grids
    nGx = 200; nGz = 100;
    if highRes, nGx = 250; nGz = 125; end
    [gX,gZ] = meshgrid(linspace(0,p.domainWidth,nGx), linspace(0,p.domainHeight,nGz));
    JG = griddata(nodeX,nodeZ,Jmag,gX,gZ);
    VG = griddata(nodeX,nodeZ,result.NodalSolution,gX,gZ);
    res = struct('J',JG, 'V',VG, 'X',gX, 'Z',gZ);

    if highRes
        [zX,zZ] = meshgrid(linspace(p.pad1L-0.005,p.pad1R+0.005,500), linspace(0,0.025,200));
        resZoom = struct('J',griddata(nodeX,nodeZ,Jmag,zX,zZ), ...
            'V',griddata(nodeX,nodeZ,result.NodalSolution,zX,zZ), 'X',zX, 'Z',zZ);
    else
        resZoom = struct('J',[], 'V',[], 'X',[], 'Z',[]);
    end

    % K = Jmax/Jmin along line 0.5mm below interface
    evalZ = p.d + 0.0005;
    [~,zi] = min(abs(gZ(:,1) - evalZ));
    xm = (gX(1,:) >= p.pad1L) & (gX(1,:) <= p.pad1R);
    Js = JG(zi,xm);  Js = Js(~isnan(Js) & Js>0);
    if ~isempty(Js), K = max(Js)/min(Js); else, K = NaN; end

    % Probe |J|
    JAtProbe = interp2(gX,gZ,JG,p.probeX,p.probeZ,'linear');
    if isnan(JAtProbe)
        dd = (gX-p.probeX).^2 + (gZ-p.probeZ).^2;
        dd(isnan(JG)) = inf;
        [~,idx] = min(dd(:));  JAtProbe = JG(idx);
    end

    % Peak |J| in tissue (excl. heart and electrode layers)
    heartM = ((gX-p.heartX).^2+(gZ-p.heartZ).^2) <= p.heartRadius^2;
    intfM  = (gZ <= p.d) | (gZ >= p.domainHeight-p.d);
    Jt = JG(~heartM & ~intfM & ~isnan(JG));
    if ~isempty(Jt), peakJtissue = max(Jt); else, peakJtissue = NaN; end
end


%% ---- evaluateDesign ----
function f = evaluateDesign(x, p)
% Objective wrapper for gamultiobj. Returns [100,100] on failure.
    try
        [~,~,K,~,peakJ] = solveElectrode(x(1), round(x(2)), x(3), p, false);
        if isnan(K) || isnan(peakJ), f = [100,100];
        else, f = [K, peakJ]; end
    catch
        f = [100, 100];
    end
end


%% ---- identifyKeySolutions ----
function [bestKIdx, kneeIdx, bestJIdx] = identifyKeySolutions(fPareto)
    [~, bestKIdx] = min(fPareto(:,1));
    [~, bestJIdx] = min(fPareto(:,2));
    % Knee: closest to origin in normalised space
    fNorm = (fPareto - min(fPareto)) ./ (max(fPareto) - min(fPareto) + eps);
    [~, kneeIdx] = min(sqrt(fNorm(:,1).^2 + fNorm(:,2).^2));
end


%% ---- extractLineProfile ----
function [dist, Jline] = extractLineProfile(resZoom, p, evalZ)
    [~,zIdx] = min(abs(resZoom.Z(:,1) - evalZ));
    xLine = resZoom.X(1,:);
    Jrow  = resZoom.J(zIdx,:);
    distFull = xLine - p.pad1X;
    mask = abs(distFull) <= p.R_pad + 0.005;
    dist = distFull(mask);  Jline = Jrow(mask);
end


%% ---- plotDefaults ----
function cfg = plotDefaults()
    cfg.fs = 10;  cfg.fn = 'Helvetica';
    cfg.colAll = [0.75 0.75 0.75];
    cfg.setAx = @(ax) set(ax, 'FontSize',10, 'FontName','Helvetica', ...
        'TickDir','out', 'LineWidth',0.5, 'Box','on');
end


%% ---- plotParetoFront ----
function plotParetoFront(allXv, allFv, xP, fP, bK, kn, bJ, Ks, pJs, cfg)
    figure('Units','cm','Position',[2 4 17.4 8],'Color','w','PaperPositionMode','auto');
    sgtitle('NSGA-II Results','FontSize',cfg.fs+1,'FontName',cfg.fn,'FontWeight','bold');

    subplot(1,2,1);
    plot(allFv(:,1), allFv(:,2), '.', 'Color',cfg.colAll, 'MarkerSize',4); hold on;
    plot(fP(:,1), fP(:,2), 'ko', 'MarkerSize',5, 'MarkerFaceColor','k');
    plot(fP(bK,1),fP(bK,2),'s','Color',[0 .5 0],'MarkerSize',10,'MarkerFaceColor',[0 .7 0],'LineWidth',1.2);
    plot(fP(bJ,1),fP(bJ,2),'d','Color',[.6 0 0],'MarkerSize',10,'MarkerFaceColor',[.9 .2 .2],'LineWidth',1.2);
    plot(fP(kn,1),fP(kn,2),'^','Color',[0 0 .6],'MarkerSize',10,'MarkerFaceColor',[.2 .2 .9],'LineWidth',1.2);
    plot(Ks, pJs, 'kx', 'MarkerSize',12, 'LineWidth',2.5); hold off;
    xlabel('K'); ylabel('Peak |J| (A/m^2)');
    title('(a) Objective space'); legend('All','Pareto','Best K','Best peakJ','Knee','Solid','Location','best','FontSize',7,'Box','off');
    cfg.setAx(gca);

    subplot(1,2,2);
    scatter3(allXv(:,1)*1e3, allXv(:,2), allXv(:,3)*1e3, 10, allFv(:,1), 'filled','MarkerFaceAlpha',0.3);
    colormap(gca,parula(64)); cb=colorbar; ylabel(cb,'K');
    hold on;
    scatter3(xP(:,1)*1e3, xP(:,2), xP(:,3)*1e3, 40, fP(:,1), 'filled','MarkerEdgeColor','k','LineWidth',0.8);
    hold off;
    xlabel('Hole radius (mm)'); ylabel('Holes/side'); zlabel('r_{perf} (mm)');
    title('(b) Decision space'); view(135,25); cfg.setAx(gca);
end


%% ---- plotSensitivity ----
function plotSensitivity(allXv, allFv, xP, fP, Ks, pJs, cfg)
    varNames  = {'Hole radius (mm)','Holes/side','r_{perf} (mm)'};
    varScales = [1000, 1, 1000];
    objNames  = {'K','Peak |J| (A/m^2)'};
    baselines = [Ks, pJs];

    figure('Units','cm','Position',[1 1 17.4 17],'Color','w','PaperPositionMode','auto');
    sgtitle('Parameter Sensitivity','FontSize',cfg.fs+1,'FontName',cfg.fn,'FontWeight','bold');

    for row = 1:2
        for col = 1:3
            subplot(2,3,(row-1)*3+col);
            scatter(allXv(:,col)*varScales(col), allFv(:,row), 8, cfg.colAll, 'filled'); hold on;
            scatter(xP(:,col)*varScales(col), fP(:,row), 20, 'k', 'filled');
            yline(baselines(row),'r--','LineWidth',0.8); hold off;
            xlabel(varNames{col}); ylabel(objNames{row});
            title(sprintf('(%s)', char('a'-1+(row-1)*3+col)));
            cfg.setAx(gca);
        end
    end
end


%% ---- plotCurrentMaps ----
function plotCurrentMaps(resSolid, resBestJ, p, KsHR, KbjHR, cfg)
    Jall = [resSolid.J(:); resBestJ.J(:)];
    clims = [0, prctile(Jall(~isnan(Jall)), 99)];
    th = linspace(0,2*pi,150);
    hx = p.heartX*100 + p.heartRadius*100*cos(th);
    hz = p.heartZ*100 + p.heartRadius*100*sin(th);

    figure('Units','cm','Position',[2 3 17.4 8],'Color','w','PaperPositionMode','auto');
    sgtitle('|J|: Solid vs Optimised','FontSize',cfg.fs+1,'FontName',cfg.fn,'FontWeight','bold');

    for kk = 1:2
        ax = subplot(1,2,kk);
        if kk==1, dat=resSolid; Kv=KsHR; lbl='(a) Solid';
        else,     dat=resBestJ; Kv=KbjHR; lbl='(b) Best peakJ'; end
        imagesc([0 50],[0 25],dat.J); colormap(ax,jet(256)); caxis(clims);
        hold on; plot(hx,hz,'w-','LineWidth',1);
        plot(p.probeX*100,p.probeZ*100,'w+','MarkerSize',6,'LineWidth',1.2); hold off;
        xlabel('x (cm)'); ylabel('z (cm)');
        title(sprintf('%s, K=%.2f',lbl,Kv));
        axis equal tight; xlim([0 50]); ylim([0 25]);
        set(ax,'YDir','normal'); cfg.setAx(ax);
    end
    cb = colorbar('Position',[0.925 0.17 0.018 0.70]);
    ylabel(cb,'|J| (A/m^2)');
end


%% ---- plotZoomedMaps ----
function plotZoomedMaps(rzSolid, rzBestJ, p, cfg)
    zAll = [rzSolid.J(:); rzBestJ.J(:)];
    zClims = [0, prctile(zAll(~isnan(zAll)), 98)];

    figure('Units','cm','Position',[2 2 17.4 8],'Color','w','PaperPositionMode','auto');
    sgtitle('Zoomed |J| Under Electrode','FontSize',cfg.fs+1,'FontName',cfg.fn,'FontWeight','bold');

    for kk = 1:2
        ax = subplot(1,2,kk);
        if kk==1, rz=rzSolid; lbl='(a) Solid';
        else,     rz=rzBestJ; lbl='(b) Best peakJ'; end
        imagesc(rz.X(1,:)*100, rz.Z(:,1)*100, rz.J);
        colormap(ax,jet(256)); caxis(zClims);
        hold on; plot([p.pad1L p.pad1R]*100,[0 0],'w-','LineWidth',2); hold off;
        xlabel('x (cm)'); ylabel('z (cm)'); title(lbl);
        xlim([p.pad1L-0.005,p.pad1R+0.005]*100); ylim([0 2.5]);
        set(ax,'YDir','normal'); cfg.setAx(ax);
    end
    cb = colorbar('Position',[0.925 0.17 0.018 0.70]);
    ylabel(cb,'|J| (A/m^2)');
end


%% ---- plotLineProfiles ----
function plotLineProfiles(distS, JlineS, distP, JlineP, p, xBJ, cfg)
    % Side-by-side
    figure('Units','cm','Position',[2 1 17.4 8],'Color','w','PaperPositionMode','auto');
    sgtitle('J Distribution Under Electrode','FontSize',cfg.fs+1,'FontName',cfg.fn,'FontWeight','bold');
    subplot(1,2,1);
    plot(distS, JlineS, '-', 'Color',[0 0 .7], 'LineWidth',1.2);
    xlabel('Distance from centre (m)'); ylabel('|J| (A/m^2)');
    title('(a) Solid'); xlim([-0.06 0.06]); cfg.setAx(gca); grid on;
    subplot(1,2,2);
    plot(distP, JlineP, '-', 'Color',[.7 0 0], 'LineWidth',1.2);
    xlabel('Distance from centre (m)'); ylabel('|J| (A/m^2)');
    title(sprintf('(b) r=%.1fmm n=%d rp=%.0fmm', xBJ(1)*1e3, round(xBJ(2)), xBJ(3)*1e3));
    xlim([-0.06 0.06]); cfg.setAx(gca); grid on;

    % Overlay
    figure('Units','cm','Position',[2 0 12 8],'Color','w','PaperPositionMode','auto');
    plot(distS*100, JlineS, '-','Color',[0 0 .7],'LineWidth',1.5); hold on;
    plot(distP*100, JlineP, '-','Color',[.7 0 0],'LineWidth',1.5);
    xline(-p.R_pad*100,':','Color',[.4 .4 .4],'LineWidth',0.7);
    xline( p.R_pad*100,':','Color',[.4 .4 .4],'LineWidth',0.7); hold off;
    xlabel('Distance from centre (cm)'); ylabel('J_S (A/m^2)');
    title('Solid vs Perforated');
    legend('Solid','Perforated','Location','south','FontSize',8,'Box','off');
    xlim([-6 6]); cfg.setAx(gca); grid on;
end

\end{lstlisting}