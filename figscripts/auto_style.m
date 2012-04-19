function [ style_array ] = auto_style( style_num )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%   style_num, number of styles to be generated
%   
% colrs = {'b','g','r','c','m','k','y'};
colrs = {'b', 'g', 'r', 'k'};
%stys = {'.', 'o', 'x', '+', '*', 's', 'd', 'v', '^', '<', '>', 'p', 'h'};
stys = {'o','+', '*', 'd', 's', 'v', '^', '<', '>'};
% lines = {'-', ':', '-.', '--'};
lines = {'-'};

styles = {1, numel(colrs)*numel(stys)*numel(lines)};
indx = 1;


for j = 1:numel(stys)
    for i = 1:numel(colrs)
        for k = 1:numel(lines)
            
            styles{1, indx} = strcat(colrs{i}, stys{j}, lines{k});
            indx = indx +1;
            
        end
    end
end
sidx = rand_int(1,numel(styles), style_num, 1,1,0);
style_cell = styles(sidx);
style_array = style_cell;
%style_array = cell2mat(style_cell);

end

