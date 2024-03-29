function i = check_complements_rec(P, A, i)
if ~isempty(P.Son1)
  i = check_complements_rec(P.Son1, A, i);
end
if ~isempty(P.Son2)
  i = check_complements_rec(P.Son2, A, i);
end

if ~isempty(P.S)
  bound = P.get_boundary();
  inter = P.get_interior();
  fbound = P.fbound;
  finter = P.finter;
  S = A(bound,bound) - A(bound,inter)*(A(inter,inter)\A(inter,bound));
  S = S([finter,fbound], [finter,fbound]);
  fprintf('Relative error of Schur complement #%d: %e \n', i, norm(S - P.S, 'fro')/norm(S,'fro'));
  i = i + 1;
end
end