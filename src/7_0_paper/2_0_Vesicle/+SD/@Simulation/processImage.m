function processImage(obj, videoTime)

	% create video from imges
	
		imageNames = dir(fullfile(obj.JPG,'*.jpg'));
		imageNames = {imageNames.name}';
	
		videoOutput = fullfile(obj.GIF,[obj.simulationName,'.gif']);
	
		ImageNum = length(imageNames);
		DelayTime = videoTime / ImageNum;
	
		for ii = 1:ImageNum
			img = imread(fullfile(obj.JPG,imageNames{ii}));
	
			% save image with transparent background
			%alphaChannel = all(img>150,3);
			%imwrite(img, fullfile(obj.PNG, [sprintf('%05d',ii), '.png']), ...
			%		'Alpha', double(~alphaChannel));
			imwrite(img, fullfile(obj.PNG, [sprintf('%05d',ii), '.png'] ));
		
			% create gif with transparent background
			[A,map] = rgb2ind(img,256);
			BGColor = double(A(1)); % background color to be set to be transparent
			if ii == 1
				%imwrite(A, map, videoOutput, 'gif', 'LoopCount', Inf, ...
				%		'DelayTime', DelayTime, 'TransparentColor', BGColor);
				imwrite(A, map, videoOutput, 'gif', 'LoopCount', Inf, 'DelayTime', DelayTime);
			else
				%imwrite(A, map, videoOutput, 'gif', 'WriteMode', 'append', ...
				%		'DelayTime', DelayTime, 'TransparentColor', BGColor);
				imwrite(A, map, videoOutput, 'gif', 'WriteMode', 'append', ...
						'DelayTime', DelayTime);
			end
	
		end
end
