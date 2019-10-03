%% make two OAM beams l = +/-1
load('OAM beam l=p=1.mat');

% number of OAM charge, if set to 0, it generates Gaussian beams
N_OAM = 1;

prb1 = zeros(257);prb1(1:end-1, 1:end-1) = LGBeams;

phase = angle(prb1) ;
phase = mod(phase * N_OAM, 2*pi);

prb1 = abs(prb1) .* exp(1i * phase);

prb2 = conj(prb1);

% prb1 = abs(prb1);
% prb2 = abs(prb2);

figure(17);
subplot(221);imagesc(abs(prb1));axis image;title('OAM +1 Amplitude');
subplot(223);imagesc(abs(prb2));axis image;title('OAM -1Amplitude');
subplot(222);imagesc(angle(prb1));axis image;title('OAM +1 Phase');
subplot(224);imagesc(angle(prb2));axis image;title('OAM -1 Phase');
pause
close(17)

%% a point defect on flat substrate

temp1 = ones(257) * (1 + 1i);
figure(17);imagesc(abs(temp1));axis image;pause

defect = 0.5 * (1 + 1i);

temp2 = temp1;
temp2(149:151,149:151) = temp2(149:151,149:151) - defect;   %0.5 * (1 + 1i);
% temp2((149:151),(149:151)+5) = temp2((149:151),(149:151)+5) + defect/1.1;

temp3 = temp1;
temp3(149:151,149:151) = temp3(149:151,149:151) + defect;
% temp3((149:151),(149:151)+5) = temp3((149:151),(149:151)+5) + defect/100;
% temp3(145:155,145:155) = 0.5 * (1 + 1i);

maskClean = temp1;

maskDefect1 = temp2;

maskDefect2 = temp3;

figure(17);
subplot(131);imagesc(abs(maskClean));
axis image;caxis([0 2]);title('clean flat substrate');
subplot(132);imagesc(abs(maskDefect1));
axis image;caxis([0 2]);title('a negative point defect on flt substrate');
subplot(133);imagesc(abs(maskDefect2));
axis image;caxis([0 2]);title('a positive point defect on flt substrate');

pause
close(17)

%% scan the probe across the flat substrate:
%   1. clean mask (a 2D contact array)
%   2. defect mask #1 (a missing contact in a 2D contact array)
%   3. defect mask #2 (a redundant contact in a a 2D contact array)

% xshifts = 0;%   -40:10;
xshifts = -75:-10;
Nshifts = length(xshifts);

% can shift all three objects the same amount in y
yshift =2;

% make a video
videoFlag = input('>> save a video??? [0/1]');
if videoFlag == 1
    videoName = input('>> Please name the video (AS A STRING, with .avi extension,): ');
    v = VideoWriter(videoName);
    v.FrameRate = 5;
    open(v);
end

% h1 = figure(18);
for i = 1:Nshifts
    
    % ================= OAM = +1 ========================    
    % a small point defect on a flat substrate
    maskDefectShifted1 = circshift(maskDefect1, [xshifts(i) yshift]);    
    ESW1_d1 = maskDefectShifted1 .* (prb1);
    Fesw1_d1 = fftshift(fft2(ifftshift(ESW1_d1)));
    
    % ================= OAM = -1 ========================
    ESW2_d1 = maskDefectShifted1 .* (prb2);
    Fesw2_d1 = fftshift(fft2(ifftshift(ESW2_d1)));
    
    % ========== difference in diffraction patterns: a small point defect on a flat substrate ===========
    temp2 = abs((Fesw2_d1)).^2 - abs(Fesw1_d1).^2;
    
    nrow = 2;
    ncol = 6;
    
%     xrange = 80:178;
%     yrange = 80:178;
    
    xrange = 1:257;
    yrange = 1:257;

    f = plotOAM(ESW1_d1,ESW2_d1);
    
%     % ========== OAM BEAM +1 ===========
%     subplot(nrow, ncol, 1);imagesc(angle(prb1(xrange, yrange)));axis image off;
%     title({[num2str(i), ' out of ', num2str(Nshifts)], 'phase of OAM beam +1'});
% %     caxis([0 0.4])
%     
%     % ========== OAM BEAM -1 ===========
%     subplot(nrow, ncol, ncol+1);imagesc(angle(prb2(xrange, yrange)));axis image off;
%     title('phase of OAM beam -1');
% %     caxis([0.5 0.9])
%     
%     % ========== object: a small point defect on a flat substrate ===========
%     subplot(nrow, ncol, [2 ncol+2]);imagesc(abs(maskDefectShifted1(xrange, yrange)));axis image off;
%     title('ESW from OAM beam -1');
% %     caxis([0.5 0.9])
% 
%     % ========== ESW: +1 OAM ===========
%     subplot(nrow, ncol, 3);imagesc(abs(ESW1_d1(xrange, yrange)));axis image off;
%     title('ESW from OAM beam +1');
% %     caxis([0 0.4])
%     
%     % ========== ESW: -1 OAM ===========
%     subplot(nrow, ncol, ncol+3);imagesc(abs(ESW2_d1(xrange, yrange)));axis image off;
%     title('ESW from OAM beam -1');
% %     caxis([0.5 0.9])
% 
%     % ========== Diff Patt: +1 OAM ===========
%     subplot(nrow, ncol, 4);imagesc(abs(Fesw1_d1(xrange, yrange)));axis image off;
%     title('Diff patt from OAM beam +1');
% %     caxis([0 0.4])
%     
%     % ========== Diff Patt: -1 OAM ===========
%     subplot(nrow, ncol, ncol+4);imagesc(abs(Fesw2_d1(xrange, yrange)));axis image off;
%     title('Diff patt from OAM beam -1');
% %     caxis([0.5 0.9])
%         
%     subplot(nrow, ncol, [5, ncol+5]);imagesc(temp2(xrange, yrange));axis image off;
%     title({'difference between OAM +1 and -1',['max value ', num2str(max(abs(temp2(:))))]});
%         
%     % ========== difference in diffraction patterns: a small point defect on a flat substrate ===========
%     temp22 = fftshift(ifft2(ifftshift(temp2)));
% %     subplot(nrow, ncol, ncol+3);imagesc( abs(temp22(80:178, 80:178)).*(mod( angle(temp22(80:178, 80:178)) - pi/2, pi)-pi/2<0));axis image;
%     
%     % amplitude of FT{differential measurement}
%     subplot(nrow, ncol, 6);imagesc( abs(temp22(xrange, yrange)) );axis image off;
%     title({'amplitude of the FT',['max value ', num2str(max(abs(temp22(:))))]});
%     
%     % phase of FT{differential measurement}
%     subplot(nrow, ncol, ncol+6);imagesc( angle(temp22(xrange, yrange)) );axis image off;
%     title('Phase of the FT');
    
    pause(0.01)
    
    if videoFlag == 1
        frame = getframe(f);
        writeVideo(v, frame);
    end
end

if videoFlag == 1
    close(v);
end


%% make 1D line gratings

% clean grating
maskClean = ones(255);
maskClean(2:5:end,:) = 1 + 1i;
maskClean(3:5:end,:) = 1 + 1i;

% defect
defect = (1+1i)/2;

% defect grating #1
maskDefect1 = maskClean;
maskDefect1(127:128, 117:119) = maskDefect1(127:128, 117:119) - defect;%(1 + 1i) / 2;% + 1i/2;
% maskDefect1(:, 117:119) = (1 + 1i) / 2;% + 1i/2;

% defect grating #2
maskDefect2 = maskClean;
maskDefect2(127:128, 117:119) = maskDefect2(127:128, 117:119) + defect;%(1 + 1i) * 1.5;% + 1i/2;
% maskDefect2(:, 117:119) = (1 + 1i) * 2;% + 1i/2;
% maskDefect2(128, 118) = (1 + 1i)/2;% + 1i/2;

figure(4);
ax(1) = subplot(131);imagesc(abs(maskClean));axis image;caxis([0.7 2.2]);
ax(2) = subplot(132);imagesc(abs(maskDefect1));axis image;caxis([0.7 2.2]);
ax(3) = subplot(133);imagesc(abs(maskDefect2));axis image;caxis([0.7 2.2]);

linkaxes(ax)    

%% scan the probe across the flat substrate:
%   1. clean mask (a 2D contact array)
%   2. defect mask #1 (a missing contact in a 2D contact array)
%   3. defect mask #2 (a redundant contact in a a 2D contact array)

% xshifts = 0;%   -40:10;
xshifts = -55:10;
Nshifts = length(xshifts);

% can shift all three objects the same amount in y
yshift = -10;

% make a video
videoFlag = input('>> save a video??? [0/1]');
if videoFlag == 1
    videoName = input('>> Please name the video (AS A STRING, with .avi extension,): ');
    v = VideoWriter(videoName);
    v.FrameRate = 5;
    open(v);
end

% h1 = figure(18);
for i = 1:Nshifts
    
    % ================= OAM = +1 ========================    
    % a small point defect on a flat substrate
    maskDefectShifted1 = circshift(maskDefect1, [xshifts(i) yshift]);    
    ESW1_d1 = maskDefectShifted1 .* (prb1(2:end-1,2:end-1));
    Fesw1_d1 = fftshift(fft2(ifftshift(ESW1_d1)));
    
    % ================= OAM = -1 ========================
    ESW2_d1 = maskDefectShifted1 .* (prb2(2:end-1,2:end-1));
    Fesw2_d1 = fftshift(fft2(ifftshift(ESW2_d1)));
    
    % ========== difference in diffraction patterns: a small point defect on a flat substrate ===========
    temp2 = abs((Fesw2_d1)).^2 - abs(Fesw1_d1).^2;
    
    nrow = 2;
    ncol = 6;
    
%     xrange = 80:178;
%     yrange = 80:178;
    
    xrange = 1:257;
    yrange = 1:257;

    f = plotOAM(ESW1_d1,ESW2_d1);
    
%     % ========== OAM BEAM +1 ===========
%     subplot(nrow, ncol, 1);imagesc(angle(prb1(xrange, yrange)));axis image off;
%     title({[num2str(i), ' out of ', num2str(Nshifts)], 'phase of OAM beam +1'});
% %     caxis([0 0.4])
%     
%     % ========== OAM BEAM -1 ===========
%     subplot(nrow, ncol, ncol+1);imagesc(angle(prb2(xrange, yrange)));axis image off;
%     title('phase of OAM beam -1');
% %     caxis([0.5 0.9])
%     
%     % ========== object: a small point defect on a flat substrate ===========
%     subplot(nrow, ncol, [2 ncol+2]);imagesc(abs(maskDefectShifted1(xrange, yrange)));axis image off;
%     title('ESW from OAM beam -1');
% %     caxis([0.5 0.9])
% 
%     % ========== ESW: +1 OAM ===========
%     subplot(nrow, ncol, 3);imagesc(abs(ESW1_d1(xrange, yrange)));axis image off;
%     title('ESW from OAM beam +1');
% %     caxis([0 0.4])
%     
%     % ========== ESW: -1 OAM ===========
%     subplot(nrow, ncol, ncol+3);imagesc(abs(ESW2_d1(xrange, yrange)));axis image off;
%     title('ESW from OAM beam -1');
% %     caxis([0.5 0.9])
% 
%     % ========== Diff Patt: +1 OAM ===========
%     subplot(nrow, ncol, 4);imagesc(abs(Fesw1_d1(xrange, yrange)));axis image off;
%     title('Diff patt from OAM beam +1');
% %     caxis([0 0.4])
%     
%     % ========== Diff Patt: -1 OAM ===========
%     subplot(nrow, ncol, ncol+4);imagesc(abs(Fesw2_d1(xrange, yrange)));axis image off;
%     title('Diff patt from OAM beam -1');
% %     caxis([0.5 0.9])
%         
%     subplot(nrow, ncol, [5, ncol+5]);imagesc(temp2(xrange, yrange));axis image off;
%     title({'difference between OAM +1 and -1',['max value ', num2str(max(abs(temp2(:))))]});
%         
%     % ========== difference in diffraction patterns: a small point defect on a flat substrate ===========
%     temp22 = fftshift(ifft2(ifftshift(temp2)));
% %     subplot(nrow, ncol, ncol+3);imagesc( abs(temp22(80:178, 80:178)).*(mod( angle(temp22(80:178, 80:178)) - pi/2, pi)-pi/2<0));axis image;
%     
%     % amplitude of FT{differential measurement}
%     subplot(nrow, ncol, 6);imagesc( abs(temp22(xrange, yrange)) );axis image off;
%     title({'amplitude of the FT',['max value ', num2str(max(abs(temp22(:))))]});
%     
%     % phase of FT{differential measurement}
%     subplot(nrow, ncol, ncol+6);imagesc( angle(temp22(xrange, yrange)) );axis image off;
%     title('Phase of the FT');
    
    pause(0.01)
    
    if videoFlag == 1
        frame = getframe(f);
        writeVideo(v, frame);
    end
end

if videoFlag == 1
    close(v);
end

%% make 2 1D line gratings with gap in between

% clean grating
maskClean = ones(255);
maskClean(2:5:end,:) = 1 + 1i;
maskClean(3:5:end,:) = 1 + 1i;
maskClean(:,117:119) = (1 + 1i) / 2;
% defect
defect = (1+1i) / 2;

% defect grating #1
maskDefect1 = maskClean;
maskDefect1(127:128, 117:119) = maskDefect1(127:128, 117:119) + defect;%(1 + 1i) / 2;% + 1i/2;
% maskDefect1(:, 117:119) = (1 + 1i) / 2;% + 1i/2;

% defect grating #2
maskDefect2 = maskClean;
maskDefect2(127:128, 117:119) = maskDefect2(127:128, 117:119) + 2 * defect;%(1 + 1i) * 1.5;% + 1i/2;
% maskDefect2(:, 117:119) = (1 + 1i) * 2;% + 1i/2;
% maskDefect2(128, 118) = (1 + 1i)/2;% + 1i/2;

figure(4);
ax(1) = subplot(131);imagesc(abs(maskClean));axis image;caxis([0.7 2.2]);
ax(2) = subplot(132);imagesc(abs(maskDefect1));axis image;caxis([0.7 2.2]);
ax(3) = subplot(133);imagesc(abs(maskDefect2));axis image;caxis([0.7 2.2]);

linkaxes(ax)

%% scan the probe across the flat substrate:
%   1. clean mask (a 2D contact array)
%   2. defect mask #1 (a missing contact in a 2D contact array)
%   3. defect mask #2 (a redundant contact in a a 2D contact array)

% xshifts = 0;%   -40:10;
xshifts = 10;
Nshifts = length(xshifts);

% can shift all three objects the same amount in y
yshift =2;

% make a video
videoFlag = input('>> save a video??? [0/1]');
if videoFlag == 1
    videoName = input('>> Please name the video (AS A STRING, with .avi extension,): ');
    v = VideoWriter(videoName);
    v.FrameRate = 5;
    open(v);
end

% h1 = figure(18);
for i = 1:Nshifts
    
    % ================= OAM = +1 ========================    
    % a small point defect on a flat substrate
    maskDefectShifted1 = circshift(maskDefect1, [xshifts(i) yshift]);    
    ESW1_d1 = maskDefectShifted1 .* (prb1(2:end-1,2:end-1));
    Fesw1_d1 = fftshift(fft2(ifftshift(ESW1_d1)));
    
    % ================= OAM = -1 ========================
    ESW2_d1 = maskDefectShifted1 .* (prb2(2:end-1,2:end-1));
    Fesw2_d1 = fftshift(fft2(ifftshift(ESW2_d1)));
    
    % ========== difference in diffraction patterns: a small point defect on a flat substrate ===========
    temp2 = abs((Fesw2_d1)).^2 - abs(Fesw1_d1).^2;
    temp3 = temp2;
%     temp3(110:146,:) = 0;
    
    nrow = 2;
    ncol = 6;
    
%     xrange = 80:178;
%     yrange = 80:178;
    
    xrange = 1:257;
    yrange = 1:257;

    f = plotOAM(ESW1_d1,ESW2_d1, Fesw1_d1, Fesw2_d1, temp3);
    
%     % ========== OAM BEAM +1 ===========
%     subplot(nrow, ncol, 1);imagesc(angle(prb1(xrange, yrange)));axis image off;
%     title({[num2str(i), ' out of ', num2str(Nshifts)], 'phase of OAM beam +1'});
% %     caxis([0 0.4])
%     
%     % ========== OAM BEAM -1 ===========
%     subplot(nrow, ncol, ncol+1);imagesc(angle(prb2(xrange, yrange)));axis image off;
%     title('phase of OAM beam -1');
% %     caxis([0.5 0.9])
%     
%     % ========== object: a small point defect on a flat substrate ===========
%     subplot(nrow, ncol, [2 ncol+2]);imagesc(abs(maskDefectShifted1(xrange, yrange)));axis image off;
%     title('ESW from OAM beam -1');
% %     caxis([0.5 0.9])
% 
%     % ========== ESW: +1 OAM ===========
%     subplot(nrow, ncol, 3);imagesc(abs(ESW1_d1(xrange, yrange)));axis image off;
%     title('ESW from OAM beam +1');
% %     caxis([0 0.4])
%     
%     % ========== ESW: -1 OAM ===========
%     subplot(nrow, ncol, ncol+3);imagesc(abs(ESW2_d1(xrange, yrange)));axis image off;
%     title('ESW from OAM beam -1');
% %     caxis([0.5 0.9])
% 
%     % ========== Diff Patt: +1 OAM ===========
%     subplot(nrow, ncol, 4);imagesc(abs(Fesw1_d1(xrange, yrange)));axis image off;
%     title('Diff patt from OAM beam +1');
% %     caxis([0 0.4])
%     
%     % ========== Diff Patt: -1 OAM ===========
%     subplot(nrow, ncol, ncol+4);imagesc(abs(Fesw2_d1(xrange, yrange)));axis image off;
%     title('Diff patt from OAM beam -1');
% %     caxis([0.5 0.9])
%         
%     subplot(nrow, ncol, [5, ncol+5]);imagesc(temp2(xrange, yrange));axis image off;
%     title({'difference between OAM +1 and -1',['max value ', num2str(max(abs(temp2(:))))]});
%         
%     % ========== difference in diffraction patterns: a small point defect on a flat substrate ===========
%     temp22 = fftshift(ifft2(ifftshift(temp2)));
% %     subplot(nrow, ncol, ncol+3);imagesc( abs(temp22(80:178, 80:178)).*(mod( angle(temp22(80:178, 80:178)) - pi/2, pi)-pi/2<0));axis image;
%     
%     % amplitude of FT{differential measurement}
%     subplot(nrow, ncol, 6);imagesc( abs(temp22(xrange, yrange)) );axis image off;
%     title({'amplitude of the FT',['max value ', num2str(max(abs(temp22(:))))]});
%     
%     % phase of FT{differential measurement}
%     subplot(nrow, ncol, ncol+6);imagesc( angle(temp22(xrange, yrange)) );axis image off;
%     title('Phase of the FT');
    
    pause(0.01)
    
    if videoFlag == 1
        frame = getframe(f);
        writeVideo(v, frame);
    end
end

if videoFlag == 1
    close(v);
end

%% make 2D contact array

temp1 = ones(8);
figure(17);imagesc(abs(temp1));axis image;pause
temp1(3:5,3:5) = 1+1i;
figure(17);imagesc(abs(temp1));axis image;pause
temp1 = repmat(temp1, 32, 32);
figure(17);imagesc(abs(temp1));axis image;pause
temp2 = ones(257);
temp2(2:end,2:end) = temp1;
figure(17);imagesc(abs(temp2));axis image;

maskClean = temp2;

% defect position
rangex = 148:150;   % centered on a contact
rangey = 148:150;

% rangex = 144:146;   % centered in between contacts
% rangey = 144:146

defect = 0.5 * (1 + 1i);

maskDefect1 = temp2;
maskDefect1(rangex,rangey) = maskDefect1(rangex,rangey) - defect;   %exp(1i*pi/4) * (1 + 1i) * 0.5;

maskDefect2 = temp2;
% maskDefect2(144:146, 144:146) = 1 + 1i;
maskDefect2(rangex,rangey) = maskDefect2(rangex,rangey) + defect;

figure(17);
subplot(131);imagesc(abs(maskClean));axis image;caxis([0.7 2.2]);
subplot(132);imagesc(abs(maskDefect1));axis image;caxis([0.7 2.2]);
subplot(133);imagesc(abs(maskDefect2));axis image;caxis([0.7 2.2]);

%% scan the probe across the flat substrate:
%   1. clean mask (a 2D contact array)
%   2. defect mask #1 (a missing contact in a 2D contact array)
%   3. defect mask #2 (a redundant contact in a a 2D contact array)

% xshifts = 0;%   -40:10;
xshifts = -70:-8;
Nshifts = length(xshifts);

% can shift all three objects the same amount in y
yshift = 4;

% make a video
videoFlag = input('>> save a video??? [0/1]');
if videoFlag == 1
    videoName = input('>> Please name the video (AS A STRING, with .avi extension,): ');
    v = VideoWriter(videoName);
    v.FrameRate = 5;
    open(v);
end

% h1 = figure(18);
for i = 1:Nshifts
    
    % ================= OAM = +1 ========================    
    % a small point defect on a flat substrate
    maskDefectShifted1 = circshift(maskDefect1, [xshifts(i) yshift]);    
%     maskDefectShifted1 = circshift(maskClean, [xshifts(i) yshift]);    
    ESW1_d1 = maskDefectShifted1 .* (prb1);
    Fesw1_d1 = fftshift(fft2(ifftshift(ESW1_d1)));
    
    % ================= OAM = -1 ========================
    ESW2_d1 = maskDefectShifted1 .* (prb2);
    Fesw2_d1 = fftshift(fft2(ifftshift(ESW2_d1)));
    
    % ========== difference in diffraction patterns: a small point defect on a flat substrate ===========
    temp2 = abs(fliplr(Fesw2_d1)).^2 - abs(Fesw1_d1).^2;
    temp3 = abs((Fesw2_d1)).^2 - abs(Fesw1_d1).^2;
%     temp3(110:146,:) = 0;
    
    nrow = 2;
    ncol = 8;
    
%     xrange = 80:178;
%     yrange = 80:178;
    
    xrange = 1:257;
    yrange = 1:257;

%     f = plotOAM(ESW1_d1,ESW2_d1, Fesw1_d1, Fesw2_d1, temp3);
    
    % ========== OAM BEAM +1 ===========
    subplot(nrow, ncol, 1);imagesc(angle(prb1(xrange, yrange)));axis image off;
    title({[num2str(i), ' out of ', num2str(Nshifts)], 'phase of OAM beam +1'});
%     caxis([0 0.4])
    
    % ========== OAM BEAM -1 ===========
    subplot(nrow, ncol, ncol+1);imagesc(angle(prb2(xrange, yrange)));axis image off;
    title('phase of OAM beam -1');
%     caxis([0.5 0.9])
    
    % ========== object: a small point defect on a flat substrate ===========
    subplot(nrow, ncol, [2 ncol+2]);imagesc(abs(maskDefectShifted1(xrange, yrange)));axis image off;
    title('ESW from OAM beam -1');
%     caxis([0.5 0.9])

    % ========== ESW: +1 OAM ===========
    subplot(nrow, ncol, 3);imagesc(abs(ESW1_d1(xrange, yrange)));axis image off;
    title('ESW from OAM beam +1');
%     caxis([0 0.4])
    
    % ========== ESW: -1 OAM ===========
    subplot(nrow, ncol, ncol+3);imagesc(abs(ESW2_d1(xrange, yrange)));axis image off;
    title('ESW from OAM beam -1');
%     caxis([0.5 0.9])

    % ========== Diff Patt: +1 OAM ===========
    subplot(nrow, ncol, 4);imagesc(abs(Fesw1_d1(xrange, yrange)));axis image off;
    title('Diff patt from OAM beam +1');
%     caxis([0 0.4])
    
    % ========== Diff Patt: -1 OAM ===========
    subplot(nrow, ncol, ncol+4);imagesc(abs(Fesw2_d1(xrange, yrange)));axis image off;
    title({'Diff patt from OAM beam -1', ['max value: ', num2str(max(abs(Fesw2_d1(:))).^2);]});
%     caxis([0.5 0.9])
        
    subplot(nrow, ncol, [5, ncol+5]);imagesc(temp2(xrange, yrange));axis image off;
    title({'difference between OAM +1 and -1',['max value ', num2str(max(abs(temp2(:))))]});
        
    % ========== difference in diffraction patterns: a small point defect on a flat substrate ===========
    temp22 = fftshift(ifft2(ifftshift(temp2)));
    temp33 = fftshift(ifft2(ifftshift(temp3)));
%     subplot(nrow, ncol, ncol+3);imagesc( abs(temp22(80:178, 80:178)).*(mod( angle(temp22(80:178, 80:178)) - pi/2, pi)-pi/2<0));axis image;
    
    % amplitude of FT{differential measurement}
    subplot(nrow, ncol, 6);imagesc( abs(temp22(xrange, yrange)) );axis image off;
    title({'amplitude of the FT',['max value ', num2str(max(abs(temp22(:))))]});
    
    % phase of FT{differential measurement}
    subplot(nrow, ncol, ncol+6);imagesc( angle(temp22(xrange, yrange)) );axis image off;
    title('Phase of the FT');
    
    subplot(nrow, ncol, [7, ncol+7]);imagesc(temp3(xrange, yrange));axis image off;
    title({'difference between OAM +1 and -1',['max value ', num2str(max(abs(temp3(:))))]});
    
    % amplitude of FT{differential measurement}
    subplot(nrow, ncol, 8);imagesc( abs(temp33(xrange, yrange)) );axis image off;
    title({'amplitude of the FT',['max value ', num2str(max(abs(temp33(:))))]});
    
    % phase of FT{differential measurement}
    subplot(nrow, ncol, ncol+8);imagesc( angle(temp33(xrange, yrange)) );axis image off;
    title('Phase of the FT');
    
    pause(0.01)
    
    if videoFlag == 1
        frame = getframe(f);
        writeVideo(v, frame);
    end
end

if videoFlag == 1
    close(v);
end