function savefigure(name,number,printfigs)
%
ffname = sprintf('./figs/%s_fig%i', name, number);

h=figure(number);
set(h,'Units','centimeters');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

if(printfigs==true)
    %savefig(ffname);
    print(sprintf('%s.pdf',ffname),'-dpdf','-r0');
end

end
