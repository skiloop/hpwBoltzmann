
clear;
clc;
dlist = dir('*.dat');

len = size(dlist,1);
for i=1:len
    if(~ dlist(i).isdir)
        load(dlist(i).name);
    end
end