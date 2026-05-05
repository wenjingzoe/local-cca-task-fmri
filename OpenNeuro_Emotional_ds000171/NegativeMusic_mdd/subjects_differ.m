%% ================== Method–CCA Functional Difference per Component Group ==================
clear; clc;
root_dir = '/Users/wenjingy/Desktop/OpenNeuro_Emotional_ds000171/NegativeMusic_mdd';
methods = ["G-FWCCA","LG-FWCCA","GL-FWCCA"]; % 与 CCA 比较
Kcomp = 8;
groups = {[1 2],[3 4],[5 6],[7 8]};
groupNames = {'Comp1-2','Comp3-4','Comp5-6','Comp7-8'};






sub_dirs = dir(fullfile(root_dir, 'sub-mdd*'));
AllDiffSys = {};

for s = 1:numel(sub_dirs)
    subj = sub_dirs(s).name;
    subj_dir = fullfile(root_dir, subj);
    if ~sub_dirs(s).isdir, continue; end

    xl = dir(fullfile(subj_dir, '*TopRegions*.xlsx'));
    if isempty(xl)
        fprintf('[Skip] No TopRegions file in %s\n', subj_dir);
        continue;
    end
    T = readtable(fullfile(xl(1).folder, xl(1).name));

    % --- 准备 CCA 的 component sets ---
    RegionSets_CCA = cell(1, Kcomp);
    for k = 1:Kcomp
        rows = contains(string(T.Method),'CCA') & (T.Comp == k);
        all_str = strjoin(T.TopRegions(rows & ~cellfun(@isempty,T.TopRegions)), ',');
        toks = lower(strtrim(strsplit(all_str, ',')));
        toks = toks(~cellfun(@isempty,toks));
        toks = toks(~contains(toks,'unknown'));
        RegionSets_CCA{k} = unique(toks);
    end

    % --- 分组后再比较 ---
    for g = 1:numel(groups)
        idx_group = groups{g};
        CCA_set = unique([RegionSets_CCA{idx_group}]); % 合并组内两个 comp 的脑区

        for m = 1:numel(methods)
            method = methods(m);

            % 合并该方法组内的所有脑区
            RegionSets_m = {};
            for k = idx_group
                rows = contains(string(T.Method), method) & (T.Comp == k);
                all_str = strjoin(T.TopRegions(rows & ~cellfun(@isempty,T.TopRegions)), ',');
                toks = lower(strtrim(strsplit(all_str, ',')));
                toks = toks(~cellfun(@isempty,toks));
                toks = toks(~contains(toks,'unknown'));
                RegionSets_m = [RegionSets_m; unique(toks(:))];
            end
            RegionSets_m = unique(RegionSets_m);

            % --- 差集（method 独有的部分） ---
            diffRegions = setdiff(RegionSets_m, CCA_set);

            % --- 分类统计 ---
            sysCount = zeros(1, numel(sysNames));
            for si = 1:numel(sysNames)
                kw = systems.(sysNames{si});
                sysCount(si) = sum(contains(diffRegions, kw));
            end

            % --- 保证 prop 合法 ---
            total = sum(sysCount);
            if total == 0 || ~isnumeric(total)
                prop = zeros(1, numel(sysNames));
            else
                prop = sysCount ./ total;
            end
            prop = double(prop(:))';  % 转为行向量

            AllDiffSys = [AllDiffSys; {subj, char(method), groupNames{g}, prop}];
        end
    end
end

%% === 转表并修复非法行 ===
Tbl = cell2table(AllDiffSys, 'VariableNames', {'Subject','Method','Group','PropVec'});

% ---- 统一 PropVec 列为 cell 数组 ----
if isnumeric(Tbl.PropVec)
    Tbl.PropVec = num2cell(Tbl.PropVec, 2);
elseif ~iscell(Tbl.PropVec)
    Tbl.PropVec = arrayfun(@(x) {x}, Tbl.PropVec);
end

% ---- 过滤非法行 ----
validRows = cellfun(@(x) isnumeric(x) && numel(x)==numel(sysNames), Tbl.PropVec);
Tbl = Tbl(validRows, :);

% ---- 拼接成矩阵 ----
PropMat = cell2mat(Tbl.PropVec);

% ---- 展开到列 ----
for i = 1:numel(sysNames)
    Tbl.(sysNames{i}) = PropMat(:,i);
end
Tbl.PropVec = [];

% === 保存结果 ===
out_file = fullfile(root_dir, 'FunctionalDiff_vsCCA_MDD.xlsx');
writetable(Tbl, out_file);
disp(['✅ Saved subject-level table: ' out_file]);

%% === Group-level mean summary ===
numVars = Tbl.Properties.VariableNames(4:end);  % 从第4列开始都是数值
MeanTbl = groupsummary(Tbl, {'Method','Group'}, 'mean', numVars);
out_mean = fullfile(root_dir, 'FunctionalDiff_vsCCA_MDD_mean.xlsx');
writetable(MeanTbl, out_mean);
disp(['✅ Saved mean table: ' out_mean]);


%% === 可视化 ===
figure('Position',[300 300 950 450]);
for m = 1:numel(methods)
    subplot(1,numel(methods),m);
    Tm = MeanTbl(strcmp(MeanTbl.Method,methods(m)),:);
    bar(categorical(Tm.Group), table2array(Tm(:,4:end)),'stacked');
    ylim([0 1]);
    legend(sysNames,'Location','eastoutside');
    title(sprintf('%s vs CCA', methods(m)),'FontWeight','bold');
    ylabel('Proportion of Diff Regions');
    xlabel('Component Group');
end

sgtitle('Functional System Differences (relative to CCA, MDD)');
%saveas(gcf, fullfile(root_dir,'FunctionalDiff_vsCCA_MDD.png'));

%% === 自动摘要输出 ===
fprintf('\n================ Summary of Functional Shifts ================\n');
for m = 1:numel(methods)
    Tm = MeanTbl(strcmp(MeanTbl.Method,methods(m)),:);
    [~, maxIdx] = max(table2array(Tm(:,4:end)),[],2); % 每个 group 哪个系统最高
    for g = 1:height(Tm)
        fprintf('%-10s | %-8s → strongest in %-20s\n', ...
            methods(m), string(Tm.Group(g)), sysNames{maxIdx(g)});
    end
end
fprintf('==============================================================\n\n');
