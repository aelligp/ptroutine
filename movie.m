%% Movie animation for thesis

writerObj = VideoWriter('Test.avi');
open(writerObj);
for n = 1 : 20000
  filename = sprintf('%d.png', n);
  thisimage = imread(filename);
  writeVideo(writerObj, thisimage);
end
close(writerObj);