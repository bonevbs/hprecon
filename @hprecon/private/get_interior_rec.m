function interior = get_interior_rec(P, interior)
  if (~isempty(P.Son1))
    interior = get_interior_rec(P.Son1, interior);
  end
  if (~isempty(P.Son2))
    interior = get_interior_rec(P.Son2, interior);
  end
  
  % throw an error if we try to append something that was already
  % eliminated
  if preconoption('checks')
    assert(isempty(intersect(interior, P.inter)), 'Found previously eliminated items in P.inter');
  end
  interior = [interior, P.inter];
end