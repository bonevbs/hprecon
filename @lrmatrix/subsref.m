function varargout = subsref(obj,s)

switch s(1).type
  case '()'
    if length(s) == 1
      if length(s.subs) == 2
        varargout{1} = lrmatrix(obj.U(s.subs{1},:), obj.V(s.subs{2},:));
      else
        varargout{1} = builtin('subsref',obj.U*obj.V',s);
      end
    else
      % Use built-in for any other expression
      [varargout{1:nargout}] = builtin('subsref',obj,s);
    end
  case '.'
    % Use built-in for any other expression
    [varargout{1:nargout}] = builtin('subsref',obj,s);
  otherwise
    error('Not a valid indexing expression')
end
end