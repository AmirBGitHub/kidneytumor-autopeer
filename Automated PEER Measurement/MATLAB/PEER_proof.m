clear all; close all; clc;

fid = fopen('PEERs.txt');
PEER_man = textscan(fid,'%s%f');
fclose(fid);

for i = 20:20
    img{i,:} = imread(['subj',num2str(i),'.jpg']);
    img_K{i,:} = imread(['subj',num2str(i),'-K.jpg']);
    img_T{i,:} = imread(['subj',num2str(i),'-T.jpg']);
    imgC{i,:} = imread(['subj',num2str(i),'Ints.jpg']); % with contrast
    imgC_K{i,:} = imread(['subj',num2str(i),'Ints-K.jpg']);
    imgC_T{i,:} = imread(['subj',num2str(i),'Ints-T.jpg']);
    img_TK_NC{i,1} = img_T{i,:}-mode(mode(img_T{i,:}));
    img_TK_NC{i,2} = img_K{i,:}-mode(mode(img_K{i,:})); 
    %- img_TK_NC{i,1};
    img_TK_C{i,1} = imgC_T{i,:}-mode(mode(imgC_T{i,:}));
    img_TK_C{i,2} = imgC_K{i,:}-mode(mode(imgC_K{i,:})); 
    %- img_TK_C{i,1};
    
    h(i) = figure;
    imshow(img{i,1})
    hold on
    saveas(gcf, 'rawNC_.jpg');
    h(i+5) = figure;
    imshow(imgC{i,1})
    hold on
    saveas(gcf, 'rawC_.jpg');
    

    for k = 1:2
        %% edge detection
        [BWoutline, BWoutlineShrinked] = edgeFind(img_TK_NC{i,k},h(i),k);
        saveas(gcf, ['edgeNC_',num2str(k),'.jpg']);
        [BWoutlineC, BWoutlineShrinkedC] = edgeFind(img_TK_C{i,k},h(i+5),k);
        saveas(gcf, ['edgeC_',num2str(k),'.jpg']);
        
        %% finding the size of the tumor
        if k == 1 % if tumor
           tumorArea(i,1) = areaFind(BWoutline);
           tumorAreaC(i,1) = areaFind(BWoutlineC);
        end
        
        %% finding the intensity and location of the brightest spot
        [centy, centx, closestDist, closestPnt] = brightSpotFind(img_TK_NC{i,k}, BWoutline, BWoutlineShrinked, h(i), k);
        [centyC, centxC, closestDistC, closestPntC] = brightSpotFind(img_TK_C{i,k}, BWoutlineC, BWoutlineShrinkedC, h(i+5), k);
        
        %% finding the direction of ellipse
        rotAngle = ellipAngleFind(BWoutline,closestPnt,h(i)); 
        rotAngleC = ellipAngleFind(BWoutlineC,closestPntC,h(i+5)); 
        
        %% drawing the shape finding the ROI intensity
        [maxBrtVal(i,k), shape] = drawShape(img_TK_NC{i,k}, tumorArea(i,1), closestDist, rotAngle, centx, centy, 'b', k, h(i));
        h1 = gcf;
        saveas(gcf, ['roiNC_',num2str(k),'.jpg']);
        [maxBrtValC(i,k), shapeC] = drawShape(img_TK_C{i,k}, tumorAreaC(i,1), closestDistC, rotAngleC, centxC, centyC, 'b', k, h(i+5));
        h2 = gcf;
        saveas(gcf, ['roiC_',num2str(k),'.jpg']);
    end
    %% salculating PEER
    PEER(i,1) = (maxBrtValC(i,1)-maxBrtVal(i,1))/(maxBrtValC(i,2)-maxBrtVal(i,2));
    
    %% save image results
    saveas(h(i),['subj',num2str(i),'_post.jpg'])
    saveas(h(i+5),['subj',num2str(i),'C_post.jpg'])
end
PEER_RMSE = sqrt(mean((PEER-PEER_man{1,2}(2:end)).^2))

BWoutline = edgeFind(img_TK_NC{i,1},h(i),1);
BWoutlineC = edgeFind(img_TK_C{i,1},h(i+5),1);

function [BWoutline, BWoutlineShrinked] = edgeFind(img,h,K)
    BWin = imbinarize(img);
    BWc = imclearborder(BWin,4);
    windowSize = 10;
    kernel = ones(windowSize) / windowSize ^ 2;
    blurryImage = conv2(single(BWc), kernel, 'same');
    BW = blurryImage > 0.3; % Rethreshold
    [B,~] = bwboundaries(BW,'noholes');

    [r, c] = find(BW);
    smallMask = imerode(BW, true(round(0.2*sqrt(range(r)*range(c)))));
    [Bsmall, ~] = bwboundaries(smallMask,'noholes');

    BWt{1,1} = Bsmall;
    BWt{2,1} = B;
    
    Bmat = [];
    for m = 1:2
        Bmat = cell2mat(BWt{m,1});
        if ~isempty(Bmat)
            x = Bmat(:,1);
            y = Bmat(:,2);
            z = ones(size(Bmat,1),1);
            BWout{m,1} = logical(accumarray([x(:),y(:)],z(:),[size(img,1) size(img,2)]));
        end
    end
    
    BWoutlineShrinked = BWout{1,1};
    BWoutline = BWout{2,1};
    
    if K == 1
        color = 'r';
    elseif K == 2
        color = 'g';
    end
    for n = 1:size(B,1)
        figure(h);
        plot(B{n,1}(:,2), B{n,1}(:,1), color, 'LineWidth', 2.5)
        hold on
    end
end

function tumorArea = areaFind(border)
    [idx_y, idx_x] = find(border); 
    maxX_r = max(idx_x);
    maxX_l = min(idx_x);
    maxY_d = max(idx_y);
    maxY_u = min(idx_y);
    [xRange, yRange] = meshgrid(maxX_l:maxX_r, maxY_u:maxY_d);
    [inTum,~] = inpolygon(xRange,yRange,idx_x,idx_y);
    tumorArea = 0.0036*numel(xRange(inTum));
end

function [centy, centx, closestDist, closestPnt] = brightSpotFind(img, border, borderShrnk, h, K)
    
    % focusing the image area inside the tumor and on the renal cortex
    borderFill = imfill(border,'holes');
    borderShrnkFill = imfill(borderShrnk,'holes');

    if K == 1 && ~isempty(borderFill)
        img(~borderShrnkFill) = 0;
    elseif K == 2 && ~isempty(borderFill) && ~isempty(borderShrnkFill)
        img(~(borderFill-borderShrnkFill)) = 0;
    end

    % get the gaussian filtered image
    v = imgaussfilt(img);
    %find the max pixel value
    max_v = max(max(v));
    % find the position of pixels having this value.
    possibles = v > 0.8*max_v;
    v_possibles = double(possibles).*double(v);
    [r, c] = find(v_possibles == min(v_possibles(v_possibles(:) > 0)));
    %[r, c] = find(v == 0.8*max_v);

    radius = [6 6];
    [X, Y] = meshgrid(-radius(1):radius(1), -radius(2):radius(2));
    theta = 0 : (2 * pi / 10000) : (2 * pi);
    pline_x = radius(1) * cos(theta);
    pline_y = radius(2) * sin(theta);
    idx_pix = inpolygon(X,Y,pline_x,pline_y);
    brtId = [];

    for j = 1:length(r)
        img_crop = img(r(j)-radius(2):r(j)+radius(2), c(j)-radius(1):c(j)+radius(1)); 
        brtId(j) = mean(double(img_crop(idx_pix)));
    end
    [~, maxBrtIdx] = max(brtId);
    [yEg, xEg] = find(border);    % location of edges
    % here, r are rows, c are columns.
    % draw shapes
    centx = c(maxBrtIdx);
    centy = r(maxBrtIdx);
    [closestDist, closestPnt] = min(sqrt(sum(([xEg yEg]'-[centx centy]').^2)));
    
%     figure(h)
%     %plot(xEg, yEg, 'y', 'LineWidth', 1.5)
%     plot(pline_x+centx,pline_y+centy,'g')
%     hold on
%     plot(xEg(closestPnt),yEg(closestPnt),'g*')
%     hold on
end

function rotAngle = ellipAngleFind(BWoutline, clsPnt, h)
    Orientations = skeletonOrientation(BWoutline,5);    %5x5 box
    Onormal = Orientations+90;    % easier to view normals
    Onr = sind(Onormal);   %vv
    Onc = cosd(Onormal);   %uu
    [r,c] = find(BWoutline);   %row/cols
    idx_ang = find(BWoutline); % Linear indices into Onr/Onc
    % Overlay normals to verify
%     figure(h);
%     quiver(c,r,-Onc(idx_ang),Onr(idx_ang));
%     hold on
    
    angle_all = atand(Onr(idx_ang)./(-Onc(idx_ang))) - 90;
    %[Onr(idx_ang(clsPnt)) -Onc(idx_ang(clsPnt))]
    st = clsPnt-1;
    nd = clsPnt+1;
    if st < 1
        angle = mean([angle_all(end+st:end); angle_all(1:nd)]);
    elseif nd > length(angle_all)
        angle = mean([angle_all(st:end); angle_all(1:nd-length(angle_all))]);
    else 
        angle = mean(angle_all(st:nd));
    end
    rotAngle = deg2rad(angle);
end

function [maxBrtValAct, shape] = drawShape(img, tumorArea, closestDist, rotAngle, centx, centy, color, K, h)
    % determining the ellipse/circle dimensions
    if K == 1
        if tumorArea >= 3
            radius = [7 7];
        elseif tumorArea < 3
            radius = [7 3.5];
        end
    elseif K == 2
        radius = round(10.*[0.5+0.5*rand(1), 0.2+0.4*rand(1)]);
    end
    xa = -radius(1): radius(1);
    xb = -radius(2): radius(2);
    theta = 0 : (2 * pi / 10000) : (2 * pi);
    pline_x = radius(1) * cos(theta);
    pline_y = radius(2) * sin(theta);
    % adjusting the location of the shape adjacent to the edges
    if closestDist < max(radius)
        centx = max(radius)*cos(rotAngle)+centx;
        centy = max(radius)*sin(rotAngle)+centy;
    end
    pline_xRot = pline_x*cos(rotAngle) - pline_y*sin(rotAngle) + centx;
    pline_yRot = pline_x*sin(rotAngle) + pline_y*cos(rotAngle) + centy;
    
    [X, Y] = meshgrid(xa, xb);
    idx_pix = inpolygon(X,Y,pline_x,pline_y);
%     [X, Y] = meshgrid(-(range(pline_xRot)/2):(range(pline_xRot)/2), -(range(pline_yRot)/2):(range(pline_yRot)/2));
%     idx_pix = inpolygon(X,Y,pline_xRot-centx,pline_yRot-centy);
    img_crop = img(round(centy-radius(2)):round(centy+radius(2)), round(centx-radius(1)):round(centx+radius(1))); 
    maxBrtValAct = mean(double(img_crop(idx_pix)));
    
    figure(h)
    shape = plot(pline_xRot, pline_yRot, [color,'-'],'LineWidth',2);
    hold on
end

function [Orientations] = skeletonOrientation(skel,blksz)
% SKELETONORIENTATION Calculate the local orientation of a skeleton
%
%Inputs: 
%  skel:  MxN binary skeleton image.
%
%  blksz: Size of block to look around for local orientation
%         Both elements must be odd and greater than or equal to three
%         blksz can be a 1x2 vector meaning [row x col] block size or 1x1
%               for square block 
%         optional, defaults to [5 5]
%
%Outputs:
%  Orientation: image that is the size of skel with zeros outside of skel
%               and orientations elsewhere
%
%

%
% Copyright 2013 The MathWorks, Inc. 
% SCd 9/25/2013
%

    %Error checking:
    assert(nargin==1||nargin==2,'One or two inputs expected');
    assert(islogical(skel),'skel should be logical');
    assert(ismatrix(skel),'skel should be a matrix');    
    sz = size(skel);
    assert(isequal(sz,size(skel)),'Sizes of M and C expected to be equal');
    
    %Handling blksz
    if nargin == 1
        %Default
        blksz = [5 5];
    else
        %assertions
        assert(all(blksz>=3),'blksz elements must be greater than or equal to three');
        assert(all(mod(blksz,2)==1),'blksz elements expected to be odd');
        if isscalar(blksz)
            blksz = blksz([1,1]);
        else
            assert(isequal(size(blksz),[1 2]),'blksz is expected to be a variable of size [1x1] or [1x2]');
        end
    end
    
    %Find the skeleton pixels' index
    [row,col] = find(skel);
    npts      = numel(row);
    
    %Pad the array and offset the rows/cols so every local block fits
    padAmount = floor(blksz./2); %distance from center to edge of block
    skelPad   = padarray(skel,padAmount); %We need to pad so that image boundary pixels are contained in a block
    
    %Preallocate Orientations
    Orientations = zeros(sz); 
    
    %Some parameters
    %-Bottom of block will be the same as center before pad
    %-Top will be bottom + block size - 1 (inclusive    
    rowHigh = row+blksz(1)-1; 
    colHigh = col+blksz(2)-1;
    center  = padAmount+1; %Center of small block
    
    %Now the engine, we will loop over the image, create each local block
    %of pixels that are touching the index of interest. Do a connected
    %components analysis, remove pixels not connected to the center pixel
    %and then calculate orientation.  
    
    %Start the engine!
    for ii = 1:npts
        %Extract small block
        block = skelPad(row(ii):rowHigh(ii),col(ii):colHigh(ii));
        
        %Label and calculate orientation
        Label = bwlabel(block);
        center_label = Label==Label(center(1),center(2)); %only label of center pixel
        rp    = regionprops(center_label,'Orientation');
        
        %Set orientation of the center pixel equal to the calculated one
        Orientations(row(ii),col(ii)) = rp.Orientation;    
    end
    
end