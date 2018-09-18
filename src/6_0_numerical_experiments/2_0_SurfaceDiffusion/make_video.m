imageNames  = dir(fullfile('imges','*.jpg'));
imageNames = {imageNames.name}';

for ii=1:length(imageNames)
	img = imread(fullfile('imges',imageNames{ii}));
	alphaChannel = all(img>150,3);
	imwrite(img, fullfile('png', [sprintf('%05d',ii), '.png']), 'Alpha', double(~alphaChannel));
end

filename = fullfile('gif','surfaceDiffusion.gif');

for ii=1:length(imageNames)
	img = imread(fullfile('imges',imageNames{ii}));
	[A,map] = rgb2ind(img,256);
	BGColor = double(A(1)); % background color
	if ii == 1
		imwrite(A, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', 1, 'TransparentColor', BGColor);
	else
		imwrite(A, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 1, 'TransparentColor', BGColor);
	end
end

