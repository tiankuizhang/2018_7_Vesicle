% this file will use .jpg files in the ../images folder to create a video

imageNames = dir(fullfile('../images','*.jpg'));
imageNames = {imageNames.name}';

outputVideo = VideoWriter(fullfile('../videos','SurfaceDiffusion.avi'));
outputVideo.FrameRate = 1;
open(outpuVideo)

for ii = 1:length(imageNames)
	img = imread(fullfile(IMG,imageNames{ii}));
	writeVideo(outputVideo,img)
end

close(outputVideo)
