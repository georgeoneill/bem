function vnorms = normstri2vert(tri,trinorms,triarea,npoints)

verts = reshape(tri,numel(tri),1);
s = size(tri);
vnorms = zeros(3,npoints);
for ii = 1:npoints
    [idx, ~] = ind2sub(s,find(verts==ii));
    vnorms(:,ii) = sum(trinorms(idx,:).*(triarea(idx)*ones(1,3)));
end

vnorms = (vnorms./(ones(3,1)*vnorm(vnorms,1)))';
