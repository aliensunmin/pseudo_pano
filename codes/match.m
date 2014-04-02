classdef match
	properties (Constant)
		pborderW = 0.10;
		plist = linspace(0,255,8);
		pnrow = 4;
		miter = 1e5;
        distFunc = str2func('match.colorHist'); % define your own feature
	end
	methods(Static)
		function pano = PanoV1(order, flipFlags, imFs, Height, folder)
			N = length(imFs);
			pano = [];
			for n = 1:N
				im = imread(imFs{order(n)});
				h = size(im,1);
				r = Height/h;
                w = round(size(im,2)*r);
				im = imresize(im, [Height w]);
				if flipFlags(n) == 1
					im = im(:,end:-1:1,:);
                end
                fs{n} = fullfile(folder,sprintf('%03d.jpg',n));
				imwrite(im,fs{n});
				if Height <= 200
					pano = [pano im];
				end
            end
            fs_str = sprintf('%s ',fs{:});
            pano_str = fullfile(folder,'pano.jpg');
            cmd = sprintf('convert %s +append %s',fs_str,pano_str);
            [status msg] = unix(cmd);
		end
        function [order, flipFlags, imFs] = getOrder(folder,distFunc) % start from this function
            if nargin < 2
                distFunc = match.distFunc;
            end
			imFs = match.getImFs(folder);
            [optRoute,minDist] = match.tspath(imFs,distFunc);
            tmp = optRoute(1:2:end)/2;
            order = round(tmp);
            flipFlags = order==tmp;
        end
		function [optRoute,minDist] = tspath(imFs,distFunc)
			D = match.getD(imFs,distFunc);
			L = size(D,1);
            xy = rand(L,2);
            [optRoute,minDist] = tsp_ga(xy,D);
		end
		function [order, imFs] = SampleOrderWrap(folder,distFunc)
            if nargin < 2
                distFunc = match.distFunc;
            end
			imFs = match.getImFs(folder);
			D = match.getD(imFs,distFunc);
			[order] = match.SampleOrder(D,match.miter);
		end
		function [order] = SampleOrder(D,miter)
			N = size(D,1)/2;
			order = 1:N;
			dist = match.distOrder(order,D);
            fprintf('%02f.',dist);
			for iter = 1:miter
				list = randperm(N);
				list = list(1:2);
				order_tmp = order;
				order_tmp(list) = order_tmp(list(end:-1:1));
				dist_tmp = match.distOrder(order_tmp,D);
				if dist_tmp < dist
					dist = dist_tmp;
					order = order_tmp;
                    fprintf('%02f.',dist);
				end
			end
		end
		function D = getD(imFs,distFunc)
			feats = match.getFeatsWrap(imFs,distFunc);
            feats = [feats{:}];
            dist = pdist(feats','cosine');
            L = size(feats,2);
            U = triu(ones(L),1);
            inds = find(U==1);
            D = zeros(L);
            D(inds) = dist;
            D = D'+D;
			inds = sub2ind(size(D),1:2:L,2:2:L);
            D(inds) = -100;
            inds = sub2ind(size(D),2:2:L,1:2:L);
            D(inds) = -100;
            inds = sub2ind(size(D),1:L,1:L);
            D(inds) = inf;
		end
		function dist = distOrder(ord,D)
			right = ord*2-1;
			left = ord*2;
			inds = sub2ind(size(D),right(1:end-1),left(2:end));
			dist = sum(D(inds));
		end
		function imFs = getImFs(folder)
			files = dir(fullfile(folder,'*.jpg'));
			N = length(files);
			for n = 1:N
				imFs{n} = fullfile(folder,files(n).name);
			end
		end
		function feats = getFeatsWrap(imFs,distFunc) 
			N = length(imFs);
			feats = {};
			% extract features
			for n = 1:N
				[pathstr,name] = fileparts(imFs{n});
				[pathstr] = fileparts(pathstr);
				file = fullfile(pathstr,'feats');
				if n == 1
					[status msg] = unix(['mkdir -p ' file]);					
				end
				file = fullfile(file,[name '.mat']);
 				try
 					tmp = load(file,'tmpFeats');
 					feats = [feats tmp.tmpFeats];
 				catch
					im = imread(imFs{n});
					tmpFeats = match.leftRightFeats(im,distFunc);
					save(file,'tmpFeats');
					feats = [feats match.leftRightFeats(im,distFunc)];
 				end
			end
		end
		function feats = leftRightFeats(im,distFunc)
% 			feats{1} = match.colorHist(im);% right
% 			feats{2} = match.colorHist(im(:,end:-1:1,:));% left
            feats{1} = distFunc(im);% right
			feats{2} = distFunc(im(:,end:-1:1,:));% left
		end
		function afeat = colorHist(im,p)
			if nargin < 2
				p.list = match.plist;
				p.bw = match.pborderW;
				p.nrow = match.pnrow;
			end
			afeat = [];
			H = size(im,1);
			hstep = round(linspace(1,H,p.nrow+1));
			W = size(im,2);
			bW = round(W*p.bw);
			for r = 1:p.nrow
				subim = im(hstep(r):hstep(r+1),(end-bW):end,:);
				nc = size(subim,3);
                feat = [];
				for c = 1:nc
					tmp = subim(:,:,c);
					feat(c,:) = histc(tmp(:),p.list);
					feat(c,:) = feat(c,:)./prod(size(subim(:,:,c)));
				end
				feat = feat(:);
				afeat = [afeat; feat];
			end
		end
	end
end
