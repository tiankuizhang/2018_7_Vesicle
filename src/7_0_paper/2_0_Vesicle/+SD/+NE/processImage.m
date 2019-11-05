function processImage(videoTime, simulationName)

	% create video from imges
	
		JPG = '../jpg';
		GIF = '../gif';
		PNG = '../png';

		imageNames = dir(fullfile(JPG,'*.jpg'));
		imageNames = {imageNames.name}';
	
		videoOutput = fullfile(GIF,[simulationName,'.gif']);
	
		ImageNum = length(imageNames);
		DelayTime = videoTime / ImageNum;
	
		for ii = 1:ImageNum
			img = imread(fullfile(JPG,imageNames{ii}));
	
			% save image with transparent background
			alphaChannel = all(img>150,3);
			imwrite(img, fullfile(PNG, [sprintf('%05d',ii), '.png']), ...
					'Alpha', double(~alphaChannel));
		
			% create gif with transparent background
			[A,map] = rgb2ind(img,256);
			BGColor = double(A(1)); % background color to be set to be transparent
			if ii == 1
				imwrite(A, map, videoOutput, 'gif', 'LoopCount', Inf, ...
						'DelayTime', DelayTime, 'TransparentColor', BGColor);
			else
				imwrite(A, map, videoOutput, 'gif', 'WriteMode', 'append', ...
						'DelayTime', DelayTime, 'TransparentColor', BGColor);
			end
	
		end
end
