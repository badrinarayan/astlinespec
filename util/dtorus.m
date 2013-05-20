function d = dtorus(a,b)
d = 0.5 - abs(0.5 - mod((a-b),1.0));
end