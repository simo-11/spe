function save_pdf(fig,filename)
    h=fig;
    set(h,'Units','Inches','PaperPositionMode','Auto',...
        'PaperUnits','Inches');
    pos = get(h,'Position');
    set(h, 'PaperSize',[pos(3), pos(4)])
    pdf_name=sprintf("gen/%s.pdf",filename);
    print(h,pdf_name,'-dpdf','-r0');
    fprintf("Saved %s\n",pdf_name);
end