classdef match
	properties (Constant)
		pborderW = 0.10;
		plist = linspace(0,255,8);
	end
	methods(Static)
        function getOrder(imFs)
            [optRoute,minDist] = tspath(imFs);
            tmp = optRoute(1:2:end)/2;
            order = round(tmp);
            flipFlags = order==tmp;
        end
		function [optRoute,minDist] = tspath(imFs)
			feats = match.getFeatsWrap(imFs);
            feats = [feats{:}];
            D = pdist(feats,'cosine');
            L = size(feats,2);
			inds = sub2ind(size(D),1:2:L,2:2:L);
            D(inds) = 0;
            inds = sub2ind(size(D),2:2:L,1:2:L);
            D(inds) = 0;
            inds = sub2ind(size(D),1:L,1:L);
            D(inds) = inf;
            xy = rand(L,2);
            [optRoute,minDist] = tsp_ga(xy,D);
		end
		function imFs = getImFs(folder)
			files = dir(fullfile(folder,'*.jpg'));
			N = length(files);
			for n = 1:N
				imFs{n} = fullfile(folder,files(n).name);
			end
		end
		function feats = getFeatsWrap(imFs)
			N = length(imFs);
			feats = {};
			% extract features
			for n = 1:N
				[pathstr,name] = fileparts(imFs{n});
				file = fullfile(pathstr,'feats');
				if n == 1
					[status msg] = unix(['mkdir -p ' file]);					
				end
				file = fullfile(file,[name '.mat']);
% 				try
% 					tmp = load(file,'feats');
% 					feats = [feats tmp.feats];
% 				catch
					im = imread(imFs{n});
					feats = [feats match.leftRightFeats(im)];
% 				end
			end
		end
		function feats = leftRightFeats(im)
			feats{1} = match.colorHist(im);% right
			feats{2} = match.colorHist(im(:,end:-1:1,:));% left
		end
		function feat = colorHist(im,p)
			if nargin < 2
				p.list = match.plist;
				p.bw = match.pborderW;
			end
			W = size(im,2);
			bW = round(W*p.bw);
			subim = im(:,(end-bW):end,:);
			nc = size(subim,3);
			for c = 1:nc
                tmp = subim(:,:,c);
				feat(c,:) = histc(tmp(:),p.list);
                feat(c,:) = feat(c,:)./prod(size(subim(:,:,c)));
            end
            feat = feat(:);
		end
	end
end
