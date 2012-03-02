function tsavefigures(path, num_figs, extra_info)
% this script is used to save the figure[1-5] into pdf files
% path, figure storage path
% numfigs, number of figures to be stored
% extra_info, extra infomation like xlabel, title etc needed

    %%
if  exist(path, 'dir') || mkdir(path)
    
    %num_figs = 8;
    for i=1:num_figs
        if ishandle(i) %if current figure i exists
            figure(i);
            set(gcf, 'PaperUnits', 'inches');
            set(gcf,'PaperSize', [8, 6]);
            set(gcf, 'PaperPositionMode', 'manual');
            %left_margin always equals bottom_margin, and set to 0.05 inches, thus the
            %fig_width = pwidth - 0.1 and fig_height = pheight - 0.1

            %position = [left_margin, bottom_margin, fig_width, fig_height];
            position = [0.05, 0.05, 7.9, 5.9];
            set(gcf, 'PaperPosition', position);
            xxtitle = get(get(gca, 'title'), 'string');
            xxtitle = sprintf('%s',xxtitle);
            
            %outfile = sprintf('%s/figure%s%d', path, xxtitle, i);
            %disp(outfile);
            
            outfile = sprintf('%s/figure_result_%d', path, i);
            disp(outfile);
            if extra_info
                if ismember(i, [3,5])
                    xlabel('time');
                    ylabel('probe');
                elseif i == 1
                    xlabel('atoms');
                    ylabel('probe');
                    colorbar;
                elseif i == 2
                    xlabel('time');
                    ylabel('atoms');
                    colorbar;
                elseif i == 4
                    xlabel('time');
                    ylabel('probe');
                    colorbar;
                end
            end

            print(gcf, '-dpdf', outfile);
            saveas(gcf, outfile, 'fig');
            %clf(gcf);
        end

    end
    
end

end