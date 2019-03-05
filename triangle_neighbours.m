function neighbour_tri = triangle_neighbours(tris,npts)
sz = size(tris);
neighbour_tri = cell(1,npts);
for ii = 1:npts
    [idx, ~] = ind2sub(sz,find(tris==ii));
    neighbour_tri{ii} = idx;
end