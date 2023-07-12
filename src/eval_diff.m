function y=eval_diff(x,pts,div)
c=1;
y=zeros(numel(x,1),1);
for i = 1:numel(pts)
    y=y+div(i)*c;
    c=c.*(x-pts(i));
end    

end