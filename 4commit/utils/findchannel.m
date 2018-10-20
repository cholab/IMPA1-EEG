function out=findchannel(chanlocs, in)
% out=findchannel(chanlocs, in)

if iscell(in)
out=size(in);
for n=1:max(size(in))
out(n)=find(ismember({chanlocs.labels},in(n)));
end
elseif ischar(in)
out=find(ismember({chanlocs.labels},{in}));
else
out={chanlocs(in).labels};
end
out=out(:);
