function save_pdf_and_fig(fig,filename)
    h=fig;
    set(h,'Units','Inches','PaperPositionMode','Auto',...
        'PaperUnits','Inches');
    pos = get(h,'Position');
    set(h, 'PaperSize',[pos(3), pos(4)])
    full_name=sprintf("gen/%s.pdf",filename);
    print(h,full_name,'-dpdf','-r0');
    fprintf("Saved %s\n",full_name);
    full_name=sprintf("gen/%s.fig",filename);
    savefig(h,full_name);
    fprintf("Saved %s\n",full_name);
end