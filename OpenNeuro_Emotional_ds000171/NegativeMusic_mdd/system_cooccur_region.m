function C = system_cooccur_region(xlsx_path, sheet_names, systems_interest)
% SYSTEM_COCOCCUR_REGION
%   统计 RegionFreq 表中 "Systems" 多标签的两两共现次数（按 Count 加权）
%
% USAGE:
%   sheets = sheetnames(xlsx_path);     % 例如 {'sub_mdd02','sub_mdd03',...}
%   C = system_cooccur_region(xlsx_path, sheets(startsWith(lower(sheets),'sub_')));
%
% INPUTS
%   xlsx_path        - RegionFreq.xlsx 的完整路径
%   sheet_names      - 要统计的 sheet 名称 cell 数组（每个 sheet 一位被试）
%   systems_interest -（可选）关心的系统名（与表中一致的拼写）
%                      默认：{'Attention','EmotionControl','Visual','Auditory',...
%                             'LanguageCognition','Memory','Motor'}
%
% OUTPUT
%   C.pairs   - N×2 cell，系统对（无向对，已排序）
%   C.count   - N×1 double，对应每个系统对的加权共现次数
%   C.table   - table，包含列 {'SysA','SysB','CoCount'}（按 CoCount 降序）
%   C.mat     - S×S 矩阵（S=systems_interest 数量），对称邻接矩阵
%   C.systems - 记录使用的 systems_interest（作为 C.mat 的行/列标签）
%   C.meta    - 一些元信息（处理的 sheet 数、有效行数等）
%
% 说明
%   - 每行的共现计数按该行的 Count 加权（若缺失/<=0，则记为 1）
%   - Systems 列支持 char/string/categorical；多标签用 , 或 ; 分隔
%   - 仅统计 systems_interest 中的标签；一行若少于 2 个有效标签则跳过
%
% 作者：你可以随意修改/扩展

    if nargin < 3 || isempty(systems_interest)
        systems_interest = {'Attention','EmotionControl','Visual','Auditory', ...
                            'LanguageCognition','Memory','Motor'};
    end
    % 统一成 cellstr，保持与 Excel 中标签完全一致的拼写
    if isstring(systems_interest); systems_interest = cellstr(systems_interest); end

    % 用 Map 累计无向对的频次
    pairMap = containers.Map('KeyType','char','ValueType','double');

    total_rows   = 0;
    valid_rows   = 0;
    processed_sh = 0;

    for s = 1:numel(sheet_names)
        try
            T = readtable(xlsx_path, 'Sheet', sheet_names{s});
        catch
            warning('读取 sheet 失败：%s（已跳过）', sheet_names{s});
            continue;
        end
        processed_sh = processed_sh + 1;
        total_rows   = total_rows + height(T);

        % 找列索引（更鲁棒）
        vnames = string(T.Properties.VariableNames);
        iSys = find(strcmpi(vnames,'Systems') | contains(lower(vnames),'system'), 1, 'first');
        iCnt = find(strcmpi(vnames,'Count')  | contains(lower(vnames),'count'),  1, 'first');
        if isempty(iSys)
            warning('未找到 Systems 列：sheet %s（已跳过）', sheet_names{s});
            continue;
        end
        if isempty(iCnt)
            warning('未找到 Count 列：sheet %s，默认权重=1', sheet_names{s});
        end

        % 取出列并强制转为 cellstr
        rawSys = T{:, iSys};
        if iscell(rawSys)
            sysCol = rawSys;
        elseif isstring(rawSys)
            sysCol = cellstr(rawSys);
        elseif iscategorical(rawSys)
            sysCol = cellstr(string(rawSys));
        else
            sysCol = cellstr(string(rawSys));
        end

        if ~isempty(iCnt)
            cntCol = T{:, iCnt};
        else
            cntCol = ones(height(T),1);
        end

        % 遍历每一行
        for r = 1:height(T)
            val = sysCol{r};
            if isempty(val) || (isstring(val) && strlength(val)==0)
                continue;
            end
            % 统一为 string 再拆分，支持逗号/分号
            sval  = string(val);
            parts = regexp(sval, '\s*[;,]\s*', 'split');   % cellstr
            tags  = strtrim(cellstr(string(parts)));
            % 去空、去重复、只保留关心的系统
            tags = tags(~cellfun('isempty', tags));
            tags = tags(ismember(tags, systems_interest));
            tags = unique(tags);
            if numel(tags) < 2
                continue;
            end

            % 行权重（优先用 Count；若缺失/<=0 -> 1）
            w = 1;
            if r <= numel(cntCol) && ~isempty(cntCol(r)) && ~isnan(cntCol(r)) && cntCol(r) > 0
                w = cntCol(r);
            end
            valid_rows = valid_rows + 1;

            % 为这一行的所有两两组合累计（无向对）
            for i = 1:numel(tags)-1
                for j = i+1:numel(tags)
                    a = tags{i}; b = tags{j};
                    key = strjoin(sort({a,b}), '|');  % 无向对的唯一 key
                    if isKey(pairMap, key)
                        pairMap(key) = pairMap(key) + w;
                    else
                        pairMap(key) = w;
                    end
                end
            end
        end
    end

    % ------- 输出结果整理 -------
    keys = pairMap.keys;
    vals = pairMap.values;
    if isempty(keys)
        % 没有共现对
        C.pairs   = cell(0,2);
        C.count   = zeros(0,1);
        C.table   = cell2table(cell(0,3), 'VariableNames', {'SysA','SysB','CoCount'});
        C.mat     = zeros(numel(systems_interest));
        C.systems = systems_interest;
        C.meta    = struct('sheets', {sheet_names}, ...
                           'processedSheets', processed_sh, ...
                           'totalRows', total_rows, ...
                           'validRows', valid_rows);
        return;
    end

    % 拆 key -> (SysA, SysB)
    pairAB = cellfun(@(k) strsplit(k, '|'), keys, 'UniformOutput', false);
    pairAB = vertcat(pairAB{:});
    counts = cell2mat(vals(:));

    % 按频次降序
    [counts, order] = sort(counts, 'descend');
    pairAB = pairAB(order, :);

    % 生成 table
    C.table = table(pairAB(:,1), pairAB(:,2), counts, ...
                    'VariableNames', {'SysA','SysB','CoCount'});
    C.pairs   = pairAB;
    C.count   = counts;

    % 生成对称邻接矩阵（便于热力图/网络图）
    S  = systems_interest;
    ns = numel(S);
    M  = zeros(ns);
    for k = 1:numel(counts)
        i = find(strcmp(S, pairAB{k,1}));
        j = find(strcmp(S, pairAB{k,2}));
        if ~isempty(i) && ~isempty(j)
            M(i,j) = M(i,j) + counts(k);
            M(j,i) = M(i,j);
        end
    end
    C.mat     = M;
    C.systems = S;

    % 元信息
    C.meta = struct('sheets', {sheet_names}, ...
                    'processedSheets', processed_sh, ...
                    'totalRows', total_rows, ...
                    'validRows', valid_rows);
end