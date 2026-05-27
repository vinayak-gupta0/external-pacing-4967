\begin{lstlisting}[language=Matlab,
    breaklines=true,
    columns=fullflexible,
    basicstyle=\ttfamily\small,
    commentstyle=\color{myblue}\ttfamily\small,
    frame=single,
    showspaces=false,
    showstringspaces=false
]
function ringElectrodeSimulation()
% Simulate the concentric-ring perforated electrode and sweep hole diameters.
% 8 rings x 2 sides = 16 discrete holes in the cross-section.
% Ref: Krasteva & Papazov (2002), BioMedical Engineering OnLine.

close all; clc;

fprintf('==========================================================\n');
fprintf(' RING ELECTRODE SIMULATION: HOLE DIAMETER SWEEP\n');
fprintf('==========================================================\n\n');

%% 1. Parameters
p.domainWidth  = 0.50;
p.domainHeight = 0.25;
p.heartRadius  = 0.045;
p.heartX = 0.25;   p.heartZ = 0.095;
p.probeX = 0.25;   p.probeZ = 0.05;
p.sigma_tissue = 0.2;
p.sigma_heart  = 0.5;
p.Vapp     = 10;
p.padWidth = 0.10;
p.pad1X = 0.25;  p.pad2X = 0.25;
p.d        = 0.001;
p.rho_metal = 5;
p.z_metal     = p.rho_metal * p.d;
p.sigma_metal = p.d / p.z_metal;
p.sigma_hole  = 1e-8;
p.pad1L = p.pad1X - p.padWidth/2;
p.pad1R = p.pad1X + p.padWidth/2;
p.pad2L = p.pad2X - p.padWidth/2;
p.pad2R = p.pad2X + p.padWidth/2;
p.R_pad = p.padWidth / 2;

[p.dl, p.externalEdges, p.elec1Edge, p.elec2Edge] = buildRingGeometry(p);

%% 2. Ring layout
ringDiameters_cm = [6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5];
ringRadii_m      = (ringDiameters_cm / 2) * 1e-2;
nHolesPerRing    = [100, 125, 150, 175, 200, 225, 250, 275];
nRings           = length(ringRadii_m);

fprintf('Ring layout (8 rings, 5 mm spacing):\n');
for ii = 1:nRings
    fprintf('  Ring %d: d=%.1f cm, r=%.1f mm, %d holes\n', ...
        ii, ringDiameters_cm(ii), ringRadii_m(ii)*1e3, nHolesPerRing(ii));
end

% Hole centres: pad centre +/- each ring radius
holeCentres_m = zeros(2*nRings, 1);
for rr = 1:nRings
    holeCentres_m(2*rr-1) = p.pad1X - ringRadii_m(rr);
    holeCentres_m(2*rr)   = p.pad1X + ringRadii_m(rr);
end
fprintf('  16 holes in cross-section\n\n');

%% 3. Diameter sweep
holeDiameters_mm = [0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0];
nSweep = length(holeDiameters_mm);

fprintf('Sweep: [%s] mm\n\n', strjoin(arrayfun(@(x) sprintf('%.1f',x), ...
    holeDiameters_mm, 'Uni', false), ', '));

%% 4. Solid baseline
fprintf('Solid baseline: ');
[resSolid, rzSolid, Ks, ~, pJs] = solveWithHoles([], [], p, true);
fprintf('K=%.2f, peakJ=%.4f\n\n', Ks, pJs);

%% 5. Sweep
results = struct('diam_mm',{},'K',{},'peakJ',{},'res',{},'rz',{},'dist',{},'Jline',{});

for ss = 1:nSweep
    holeR_m = holeDiameters_mm(ss) * 1e-3 / 2;

    if holeR_m > 0.00125
        fprintf('  WARNING: d=%.1f mm may overlap (spacing=2.5 mm)\n', holeDiameters_mm(ss));
    end

    holeLArr = holeCentres_m - holeR_m;
    holeRArr = holeCentres_m + holeR_m;

    fprintf('  d=%.1f mm ... ', holeDiameters_mm(ss));
    [res, rz, K, ~, peakJ] = solveWithHoles(holeLArr, holeRArr, p, true);
    fprintf('K=%.2f, peakJ=%.4f\n', K, peakJ);

    evalZ = p.d + 0.0005;
    [dist, Jline] = extractProfile(rz, p, evalZ);

    results(ss) = struct('diam_mm',holeDiameters_mm(ss), 'K',K, 'peakJ',peakJ, ...
        'res',res, 'rz',rz, 'dist',dist, 'Jline',Jline);
end

%% 6. Baseline profile
[distSolid, JlineSolid] = extractProfile(rzSolid, p, p.d + 0.0005);

%% 7. Results table
fprintf('\n  %-10s  %-8s  %-12s  %-8s  %-8s\n', 'd(mm)','K','peakJ','dK(%)','dJ(%)');
fprintf('  %-10s  %-8.2f  %-12.4f\n', 'Solid', Ks, pJs);
for ss = 1:nSweep
    dK = (results(ss).K - Ks)/Ks*100;
    dJ = (results(ss).peakJ - pJs)/pJs*100;
    fprintf('  %-10.1f  %-8.2f  %-12.4f  %-+8.1f  %-+8.1f\n', ...
        results(ss).diam_mm, results(ss).K, results(ss).peakJ, dK, dJ);
end
fprintf('\n');

%% 8. Figures
fntSz = 10; fntName = 'Helvetica';
setAx = @(ax) set(ax, 'FontSize',fntSz, 'FontName',fntName, ...
    'TickDir','out', 'LineWidth',0.5, 'Box','on');

sweepCols = [linspace(0,0.85,nSweep); linspace(0.2,0.1,nSweep); linspace(0.7,0.1,nSweep)]';

% Fig 1: All line profiles
figure('Units','cm','Position',[2 4 17.4 10],'Color','w','PaperPositionMode','auto');
plot(distSolid*100, JlineSolid, 'k--', 'LineWidth',2.0); hold on;
leg = {'Solid'};
for ss = 1:nSweep
    lw = 1.0; if results(ss).diam_mm==0.4, lw=2.2; end
    plot(results(ss).dist*100, results(ss).Jline, '-', 'Color',sweepCols(ss,:), 'LineWidth',lw);
    leg{end+1} = sprintf('d=%.1fmm', results(ss).diam_mm); %#ok<AGROW>
end
xline(-p.R_pad*100,':','Color',[.4 .4 .4],'LineWidth',0.7);
xline( p.R_pad*100,':','Color',[.4 .4 .4],'LineWidth',0.7); hold off;
xlabel('Distance from pad centre (cm)'); ylabel('|J| (A/m^2)');
title('Hole Diameter Sweep: Line Profiles');
legend(leg,'Location','southoutside','NumColumns',3,'FontSize',7,'Box','off');
xlim([-6 6]); setAx(gca); grid on;

% Fig 2: K and peakJ vs diameter
diams = [results.diam_mm];  Kvec = [results.K];  pJvec = [results.peakJ];
bIdx = find(diams==0.4);

figure('Units','cm','Position',[2 3 17.4 8],'Color','w','PaperPositionMode','auto');
sgtitle('Effect of Hole Diameter','FontSize',fntSz+1,'FontName',fntName,'FontWeight','bold');

subplot(1,2,1);
plot(diams, Kvec, 'ko-','MarkerFaceColor',[.2 .4 .8],'MarkerSize',6,'LineWidth',1.2);
hold on; yline(Ks,'r--','LineWidth',1);
if ~isempty(bIdx), plot(0.4,Kvec(bIdx),'rs','MarkerSize',12,'MarkerFaceColor',[1 .6 .6],'LineWidth',1.5); end
hold off;
xlabel('Hole diameter (mm)'); ylabel('K');
title('(a) Non-uniformity'); legend('Ring','Solid','0.4mm','Location','best','FontSize',7,'Box','off');
setAx(gca); grid on;

subplot(1,2,2);
plot(diams, pJvec, 'ko-','MarkerFaceColor',[.8 .3 .2],'MarkerSize',6,'LineWidth',1.2);
hold on; yline(pJs,'r--','LineWidth',1);
if ~isempty(bIdx), plot(0.4,pJvec(bIdx),'rs','MarkerSize',12,'MarkerFaceColor',[1 .6 .6],'LineWidth',1.5); end
hold off;
xlabel('Hole diameter (mm)'); ylabel('Peak |J| (A/m^2)');
title('(b) Peak tissue |J|'); legend('Ring','Solid','0.4mm','Location','best','FontSize',7,'Box','off');
setAx(gca); grid on;

% Fig 3: Full-domain |J| maps (solid vs 0.4mm)
if ~isempty(bIdx)
    th = linspace(0,2*pi,150);
    hx = p.heartX*100 + p.heartRadius*100*cos(th);
    hz = p.heartZ*100 + p.heartRadius*100*sin(th);
    Jall = [resSolid.J(:); results(bIdx).res.J(:)];
    clims = [0, prctile(Jall(~isnan(Jall)),99)];

    figure('Units','cm','Position',[2 2 17.4 8],'Color','w','PaperPositionMode','auto');
    sgtitle('|J| Maps: Solid vs Ring (0.4mm)','FontSize',fntSz+1,'FontName',fntName,'FontWeight','bold');

    for kk = 1:2
        ax = subplot(1,2,kk);
        if kk==1, dat=resSolid; lbl=sprintf('(a) Solid, K=%.2f',Ks);
        else, dat=results(bIdx).res; lbl=sprintf('(b) Ring 0.4mm, K=%.2f',results(bIdx).K); end
        imagesc([0 50],[0 25],dat.J); colormap(ax,jet(256)); caxis(clims);
        hold on; plot(hx,hz,'w-','LineWidth',1);
        plot(p.probeX*100,p.probeZ*100,'w+','MarkerSize',6,'LineWidth',1.2); hold off;
        xlabel('x (cm)'); ylabel('z (cm)'); title(lbl);
        axis equal tight; xlim([0 50]); ylim([0 25]);
        set(ax,'YDir','normal'); setAx(ax);
    end
    cb=colorbar('Position',[0.925 0.17 0.018 0.70]); ylabel(cb,'|J| (A/m^2)');
end

% Fig 4: Zoomed maps
if ~isempty(bIdx)
    zAll = [rzSolid.J(:); results(bIdx).rz.J(:)];
    zClims = [0, prctile(zAll(~isnan(zAll)),98)];

    figure('Units','cm','Position',[2 1 17.4 8],'Color','w','PaperPositionMode','auto');
    sgtitle('Zoomed |J| Under Electrode','FontSize',fntSz+1,'FontName',fntName,'FontWeight','bold');

    for kk = 1:2
        ax = subplot(1,2,kk);
        if kk==1, rz=rzSolid; lbl='(a) Solid';
        else, rz=results(bIdx).rz; lbl='(b) Ring 0.4mm'; end
        imagesc(rz.X(1,:)*100, rz.Z(:,1)*100, rz.J);
        colormap(ax,jet(256)); caxis(zClims);
        hold on; plot([p.pad1L p.pad1R]*100,[0 0],'w-','LineWidth',2); hold off;
        xlabel('x (cm)'); ylabel('z (cm)'); title(lbl);
        xlim([p.pad1L-0.005,p.pad1R+0.005]*100); ylim([0 2.5]);
        set(ax,'YDir','normal'); setAx(ax);
    end
    cb=colorbar('Position',[0.925 0.17 0.018 0.70]); ylabel(cb,'|J| (A/m^2)');
end

% Fig 5: Solid vs 0.4mm overlay
figure('Units','cm','Position',[2 0 12 8],'Color','w','PaperPositionMode','auto');
plot(distSolid*100, JlineSolid, 'k-', 'LineWidth',1.8); hold on;
if ~isempty(bIdx)
    plot(results(bIdx).dist*100, results(bIdx).Jline, '-','Color',[.8 .15 .15],'LineWidth',1.8);
end
for rr = 1:nRings
    xline(-ringRadii_m(rr)*100,':','Color',[.6 .6 .6],'LineWidth',0.5);
    xline( ringRadii_m(rr)*100,':','Color',[.6 .6 .6],'LineWidth',0.5);
end
xline(-p.R_pad*100,'-','Color',[.3 .3 .3],'LineWidth',0.8);
xline( p.R_pad*100,'-','Color',[.3 .3 .3],'LineWidth',0.8); hold off;
xlabel('Distance from pad centre (cm)'); ylabel('J_S (A/m^2)');
title('Solid vs Ring (d=0.4mm)');
legend('Solid','Ring (0.4mm)','Location','south','FontSize',8,'Box','off');
xlim([-6 6]); setAx(gca); grid on;

% Fig 6: % change
figure('Units','cm','Position',[2 0 12 8],'Color','w','PaperPositionMode','auto');
dK_pct = ([results.K]-Ks)/Ks*100;
dJ_pct = ([results.peakJ]-pJs)/pJs*100;
plot(diams,dK_pct,'o-','Color',[.2 .4 .8],'MarkerFaceColor',[.2 .4 .8],'MarkerSize',5,'LineWidth',1.2); hold on;
plot(diams,dJ_pct,'s-','Color',[.8 .3 .2],'MarkerFaceColor',[.8 .3 .2],'MarkerSize',5,'LineWidth',1.2);
yline(0,'k:','LineWidth',0.7);
if ~isempty(bIdx), xline(0.4,'--','Color',[.4 .4 .4],'LineWidth',0.8); end
hold off;
xlabel('Hole diameter (mm)'); ylabel('Change from solid (%)');
title('Performance vs Hole Diameter');
legend('\DeltaK','\Deltapeak|J|','Location','best','FontSize',8,'Box','off');
setAx(gca); grid on;

%% 9. Save
save('ring_electrode_results.mat', ...
    'results','resSolid','rzSolid','Ks','pJs', ...
    'ringRadii_m','nHolesPerRing','ringDiameters_cm', ...
    'holeDiameters_mm','distSolid','JlineSolid','p');
fprintf('Results saved. Figures generated.\n');
end


%% Forward solver: discrete insulating holes
function [res, resZoom, K, JAtProbe, peakJtissue] = solveWithHoles(holeLArr, holeRArr, p, highRes)
    nHoles = length(holeLArr);
    useHoles = nHoles > 0;

    function tf = inHole(xi)
        if ~useHoles, tf = false;
        else, tf = any((xi >= holeLArr) & (xi <= holeRArr)); end
    end

    mdl = createpde();
    geometryFromEdges(mdl, p.dl);
    if highRes
        generateMesh(mdl,'Hmax',0.0012,'Hmin',0.0003,'GeometricOrder','linear');
    else
        generateMesh(mdl,'Hmax',0.002,'Hmin',0.0005,'GeometricOrder','linear');
    end

    function sig = condCoeff(region, ~)
        x = region.x;  z = region.y;
        sig = p.sigma_tissue * ones(size(x));
        for ii = 1:numel(x)
            xi = x(ii);  zi = z(ii);
            if (xi-p.heartX)^2+(zi-p.heartZ)^2 <= p.heartRadius^2
                sig(ii) = p.sigma_heart;
            elseif zi<=p.d && xi>=p.pad1L && xi<=p.pad1R
                if inHole(xi), sig(ii)=p.sigma_hole; else, sig(ii)=p.sigma_metal; end
            elseif zi>=p.domainHeight-p.d && xi>=p.pad2L && xi<=p.pad2R
                if inHole(xi), sig(ii)=p.sigma_hole; else, sig(ii)=p.sigma_metal; end
            end
        end
    end

    specifyCoefficients(mdl,'m',0,'d',0,'c',@condCoeff,'a',0,'f',0);
    applyBoundaryCondition(mdl,'dirichlet','Edge',p.elec1Edge,'u',0);
    applyBoundaryCondition(mdl,'dirichlet','Edge',p.elec2Edge,'u',p.Vapp);
    for eID = setdiff(p.externalEdges,[p.elec1Edge,p.elec2Edge])
        applyBoundaryCondition(mdl,'neumann','Edge',eID,'g',0,'q',0);
    end

    result = solvepde(mdl);
    nodeX = mdl.Mesh.Nodes(1,:)';  nodeZ = mdl.Mesh.Nodes(2,:)';
    [gradX,gradZ] = evaluateGradient(result,nodeX,nodeZ);

    sigN = p.sigma_tissue*ones(size(nodeX));
    for ii = 1:numel(nodeX)
        xi = nodeX(ii);  zi = nodeZ(ii);
        if (xi-p.heartX)^2+(zi-p.heartZ)^2 <= p.heartRadius^2
            sigN(ii) = p.sigma_heart;
        elseif zi<=p.d && xi>=p.pad1L && xi<=p.pad1R
            if inHole(xi), sigN(ii)=p.sigma_hole; else, sigN(ii)=p.sigma_metal; end
        elseif zi>=p.domainHeight-p.d && xi>=p.pad2L && xi<=p.pad2R
            if inHole(xi), sigN(ii)=p.sigma_hole; else, sigN(ii)=p.sigma_metal; end
        end
    end
    Jmag = sqrt((-sigN.*gradX).^2 + (-sigN.*gradZ).^2);

    % Interpolate onto grids
    [gX,gZ] = meshgrid(linspace(0,p.domainWidth,250), linspace(0,p.domainHeight,125));
    JG = griddata(nodeX,nodeZ,Jmag,gX,gZ);
    VG = griddata(nodeX,nodeZ,result.NodalSolution,gX,gZ);
    res = struct('J',JG, 'V',VG, 'X',gX, 'Z',gZ);

    [zX,zZ] = meshgrid(linspace(p.pad1L-0.005,p.pad1R+0.005,500), linspace(0,0.025,200));
    resZoom = struct('J',griddata(nodeX,nodeZ,Jmag,zX,zZ), ...
        'V',griddata(nodeX,nodeZ,result.NodalSolution,zX,zZ), 'X',zX, 'Z',zZ);

    % K = Jmax/Jmin at 0.5mm below interface
    evalZ = p.d + 0.0005;
    [~,zi] = min(abs(gZ(:,1)-evalZ));
    xm = (gX(1,:)>=p.pad1L) & (gX(1,:)<=p.pad1R);
    Js = JG(zi,xm);  Js = Js(~isnan(Js) & Js>0);
    if ~isempty(Js), K = max(Js)/min(Js); else, K = NaN; end

    JAtProbe = interp2(gX,gZ,JG,p.probeX,p.probeZ,'linear');
    if isnan(JAtProbe)
        dd = (gX-p.probeX).^2+(gZ-p.probeZ).^2; dd(isnan(JG))=inf;
        [~,idx] = min(dd(:)); JAtProbe = JG(idx);
    end

    heartM = ((gX-p.heartX).^2+(gZ-p.heartZ).^2) <= p.heartRadius^2;
    intfM = (gZ<=p.d) | (gZ>=p.domainHeight-p.d);
    Jt = JG(~heartM & ~intfM & ~isnan(JG));
    if ~isempty(Jt), peakJtissue = max(Jt); else, peakJtissue = NaN; end
end


%% Extract 1D profile
function [dist, Jline] = extractProfile(resZoom, p, evalZ)
    [~,zIdx] = min(abs(resZoom.Z(:,1)-evalZ));
    distFull = resZoom.X(1,:) - p.pad1X;
    Jrow = resZoom.J(zIdx,:);
    mask = abs(distFull) <= p.R_pad+0.005;
    dist = distFull(mask);  Jline = Jrow(mask);
end


%% Build geometry
function [dl, externalEdges, elec1Edge, elec2Edge] = buildRingGeometry(p)
    vertices = [0,0; p.pad1L,0; p.pad1R,0; p.domainWidth,0;
                p.domainWidth,p.domainHeight; p.pad2R,p.domainHeight;
                p.pad2L,p.domainHeight; 0,p.domainHeight];
    nv = size(vertices,1);
    pgon = [2;nv;vertices(:,1);vertices(:,2)];
    circ = [1;p.heartX;p.heartZ;p.heartRadius;zeros(6,1)];
    maxLen = max(length(pgon),length(circ));
    pgon(end+1:maxLen)=0; circ(end+1:maxLen)=0;

    [dl,~] = decsg([pgon,circ], 'P1+C1', char('P1','C1')');

    numEdges = size(dl,2);  externalEdges = [];
    for eID = 1:numEdges
        if dl(6,eID)==0 || dl(7,eID)==0, externalEdges(end+1)=eID; end %#ok<AGROW>
    end

    elec1Edge=2; elec2Edge=6; tol=1e-6;
    for eID = externalEdges
        mx=(dl(2,eID)+dl(3,eID))/2; my=(dl(4,eID)+dl(5,eID))/2;
        if abs(my)<tol && mx>p.pad1L-tol && mx<p.pad1R+tol, elec1Edge=eID; end
        if abs(my-p.domainHeight)<tol && mx>p.pad2L-tol && mx<p.pad2R+tol, elec2Edge=eID; end
    end
    fprintf('Geometry: %d edges, elec1=%d, elec2=%d\n\n', numEdges, elec1Edge, elec2Edge);
end

\end{lstlisting}