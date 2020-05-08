function varargout = subsref(obj,s)

switch s(1).type
  case '()'
    if length(s) == 1
      if length(s.subs) == 2
        II = s.subs{1};
        JJ = s.subs{2};
        [Is, Ip] = sort(II); 
        [Js, Jp] = sort(JJ);
        iIp(Ip) = [1:length(II)];
        iJp(Jp) = [1:length(JJ)];
        NI = size(obj.A11,1);
        NJ = size(obj.A11,2);

        % maybe in the future this functionality will be moved to get(.)
        if isa(obj.A11, 'hss')
          A11 = get(obj.A11,Is(Is<=NI),Js(Js<=NJ));
        else
          A11 = full(obj.A11(Is(Is<=NI),Js(Js<=NJ)));
        end
        if isa(obj.A12, 'hss')
          A12 = get(obj.A12,Is(Is<=NI),Js(Js>NJ)-NJ);
        else
          A12 = full(obj.A12(Is(Is<=NI),Js(Js>NJ)-NJ));
        end
        if isa(obj.A21, 'hss')
          A21 = get(obj.A21,Is(Is>NI)-NI,Js(Js<=NJ));
        else
          A21 = full(obj.A21(Is(Is>NI)-NI,Js(Js<=NJ)));
        end
        if isa(obj.A22, 'hss')
          A22 = get(obj.A22,Is(Is>NI)-NI,Js(Js>NJ)-NJ);
        else
          A22 = full(obj.A22(Is(Is>NI)-NI,Js(Js>NJ)-NJ));
        end
        M = [A11, A12; A21, A22];
        varargout{1} = M(iIp, iJp);
      else
        varargout{1} = builtin('subsref',full(obj),s);
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