function dict_draw(res_maps, path)
% res_maps, map with the learned dictionary and coefficients matrix inside
% path, path to save the figures using tsavefigures() below, default is
% commented

fig_indx =0;
mkeys = res_maps.keys();
for imp = 1:size(mkeys,2)
    mdat = res_maps(mkeys{imp});
    outD = mdat{1};
    outX = mdat{2};

    dict_title = sprintf('dict_%s', mkeys{imp});
    cofs_title = sprintf('cofs_%s', mkeys{imp});
    
    fig_indx = fig_indx +1;
    figure(fig_indx);
    clf(gcf);
    pcolor(outD);
    shading flat;
    title(dict_title);
    ylabel('probes');
    xlabel('atoms');
    
    fig_indx = fig_indx +1;
    figure(fig_indx);
    clf(gcf);
    pcolor(outX);
    shading flat;
    title(cofs_title);
    xlabel('time');
    ylabel('atoms');
   
    break;

end

%tsavefigures(path, size(mkeys,2)*2, 0);

end