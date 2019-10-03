function f = plotOAM(esw1,esw2,det1,det2,diff,def)
if nargin<3
    det1 = abs(fftshift(fft2(ifftshift(esw1)))).^2;
end
if nargin<4
    det2 = abs(fftshift(fft2(ifftshift(esw2)))).^2;
end
if nargin<5
%     diff = det2-det1;
    diff = fftshift(ifft2(ifftshift(det2))) - fftshift(ifft2(ifftshift(det1)));
end
if nargin<6
    def = fftshift(ifft2(ifftshift(diff)));
end


sz1 = 4;
sz2 = 12;
v = [1:2 [1:2]+sz2];
f=figure(1);
clf('reset')
set(gcf,'color','w')
ax(1) = subplot_er(sz1,sz2,v);imagesc(abs(esw1));axis image xy;title('Defect','FontSize',18);ylabel('L: +1','FontSize',16)
ax(2) = subplot_er(sz1,sz2,v+2*sz2);imagesc(abs(esw2));axis image xy;ylabel('L: -1','FontSize',16)

ax(3) = subplot_er(sz1,sz2,v+3);imagesc(log(abs(det1)));axis image xy off;title(['Diffraction Patterns, max count: ', num2str(max(abs(det1(:))))],'FontSize',18)
ax(4) = subplot_er(sz1,sz2,v+3+2*sz2);imagesc(log(abs(det2)));axis image xy off

ax(5) = subplot_er(sz1,sz2,v+7+sz2);imagesc(abs(diff));axis image xy off;title(['Difference, max count: ', num2str(max(abs(diff(:))))],'FontSize',18)

ax(6) = subplot_er(sz1,sz2,v+10+ 0);imagesc(abs(def));axis image xy ;title('Defect Location','FontSize',18);ylabel('Amplitude','FontSize',16)
ax(7) = subplot_er(sz1,sz2,v+10+2*sz2);imagesc(angle(def));axis image xy ;ylabel('Phase','FontSize',16)

annotation('textbox',[0.185 0.755 0.5 0.5],'String','Propagate','FitBoxToText','on','FontSize',18,'color','r','EdgeColor','none','Interpreter','Latex','VerticalAlignment','bottom')
annotation('arrow',[0.2 0.25],[0.75 0.75],'Color',[1 0 0],'LineWidth',5,'HeadLength',18,'HeadWidth',18)
annotation('arrow',[0.2 0.25],[0.25 0.25],'Color',[1 0 0],'LineWidth',5,'HeadLength',18,'HeadWidth',18)

% a = listfonts;
% for i = 1:length(a)
    annotation('textbox',[0.35 0.02 0.5 0.5],'String','\}','FitBoxToText','on','FontSize',525,'FontName','DejaVu Sans Light','color','r','EdgeColor','none','VerticalAlignment','bottom')
%     pause(0.01)
% end
annotation('line',[0.5 0.525]+.05,[0.5 0.5]+.019,'Color',[1 0 0],'LineWidth',5)

annotation('arrow',[0.75 0.8]+.025,[0.5 0.5]+.007,'Color',[1 0 0],'LineWidth',5,'HeadLength',18,'HeadWidth',18)
annotation('textbox',[0.78 0.505 0.5 0.5],'String','$\mathcal{F}^{-1}$','FitBoxToText','on','FontSize',18,'color','r','EdgeColor','none','Interpreter','Latex','VerticalAlignment','bottom')


% s = listfonts;
% for i = 1:length(s)
%     set(x, 'FontName', s{i});
%     title(s{i});
%     pause
% end


