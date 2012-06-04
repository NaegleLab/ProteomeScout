function bar_sites(y, x)

X = ['r', 'g', 'b', 'm'];
[a,b,c,d,e,f] =csvread2('for_matlab.csv', '\t');
j
for(i = 1:68)
  subplot(x,y,i);
  Y(1,:) = cell2mat(f(i,2:5));
  Y(2,:) = zeros(1,4);
  bar(Y);
  set(gca, 'FontName', 'Arial');
  set(gca, 'XTick', []);
  title([f{i,1}]);
  ax = axis;
  ax(1) = 0.6;
  ax(2) = 1.4;
  axis(ax);
end



