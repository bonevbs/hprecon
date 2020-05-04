function B = full(A)
  B = [full(A.A11), full(A.A12); full(A.A21), full(A.A22)];
end