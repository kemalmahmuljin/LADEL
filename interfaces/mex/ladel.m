classdef ladel < handle
    %Matlab interface for the LADEL C library
    %
    %ladel methods
    %
    %
    %
    %
    
    properties (SetAccess = private, Hidden = true)
        ncol;
    end
    
    methods
        
        %% Constructor - Create a new solver instance
        function this = ladel(ncol)
            ladel_mex('init', ncol);
        end
        
        function delete(~)
            ladel_mex('delete');
        end
        
        function varargout = factorize(~, M, varargin)
            M = triu(M);
            if nargin == 3
                ordering = varargin{1};
                ladel_mex('factorize', M, ordering);
            else
                ladel_mex('factorize', M);
            end
            
            if nargout > 0
                if nargout == 2
                    [L, D] = ladel_mex('return');
                    varargout{1} = L;
                    varargout{2} = D;
                elseif nargout == 3
                    [L, D, p] = ladel_mex('return');
                    varargout{1} = L;
                    varargout{2} = D;
                    varargout{3} = p;
                else
                    error('Wrong number of output arguments for factorize');
                end
            end
        end
        
        function y = dense_solve(~, x)
            y = ladel_mex('solve', x);
        end
        
%         function y = symmetric_matvec(~, M, x)
%             y = ladel_mex('symmetric_matvec', M, x);
%         end
        
        function varargout = factorize_advanced(~, M, Mbasis, varargin)
            M = triu(M);
            Mbasis = triu(Mbasis);
            if nargin == 4
                ordering = varargin{1};
                ladel_mex('factorize_advanced', M, Mbasis, ordering);
            else
                ladel_mex('factorize_advanced', M, Mbasis);
            end
            
            if nargout > 0
                if nargout == 2
                    [L, D] = ladel_mex('return');
                    varargout{1} = L;
                    varargout{2} = D;
                elseif nargout == 3
                    [L, D, p] = ladel_mex('return');
                    varargout{1} = L;
                    varargout{2} = D;
                    varargout{3} = p;
                else
                    error('Wrong number of output arguments for factorize_advanced');
                end
            end
                    
        end

        function varargout = factorize_advanced_with_fixed_part(~, M, Mbasis, varargin)
            M = triu(M);
            Mbasis = triu(Mbasis);
            if nargin == 5
                ordering = varargin{1};
                fixed_num = varargin{2};
                if ordering ==2
                    error('factorize_advanced_with_fixed_part should be called with 7 input arguments if ordering==2');
                end
                ladel_mex('factorize_advanced_with_fixed_part', M, Mbasis, ordering, fixed_num);
            elseif nargin ==7
                ordering = varargin{1};
                fixed_num = varargin{2};
                beta = varargin{3};
                n = varargin{4};
                ladel_mex('factorize_advanced_with_fixed_part', M, Mbasis, ordering, fixed_num, beta, n);        
            elseif nargin == 8
                ordering = varargin{1};
                fixed_num = varargin{2};
                beta = varargin{3};
                n = varargin{4};
                REG = varargin{5};
                ladel_mex('factorize_advanced_with_fixed_part', M, Mbasis, ordering, fixed_num, beta, n, REG);        
            else
                error('Wrong number of input arguments for factorize_advanced_with_fixed_part');
            end
            
            if nargout > 0
                if nargout == 1
                    E = ladel_mex('return');
                    varargout{1} = E;
                elseif nargout == 2
                    [L, D] = ladel_mex('return');
                    varargout{1} = L;
                    varargout{2} = D;
                elseif nargout == 3
                    [L, D, p] = ladel_mex('return');
                    varargout{1} = L;
                    varargout{2} = D;
                    varargout{3} = p;
                elseif nargout == 4
                    [L, D, p, E] = ladel_mex('return');
                    varargout{1} = L;
                    varargout{2} = D;
                    varargout{3} = p;
                    varargout{4} = E;
                else
                    error('Wrong number of output arguments for factorize_advanced');
                end
            end
                    
        end
%
        function varargout = factorize_with_prior_basis(~, M, varargin)
            M = triu(M);         
            if nargin == 2
                ladel_mex('factorize_with_prior_basis', M);
            elseif nargin == 4
                beta = varargin{1};
                n = varargin{2};
                ladel_mex('factorize_with_prior_basis', M, beta, n);
            elseif nargin == 5
                beta = varargin{1};
                n = varargin{2};
                tol = varargin{3};
                ladel_mex('factorize_with_prior_basis', M, beta, n, tol);
            else 
                error('Wrong number of output arguments for factorize_with_prior_basis');
            end
            if nargout > 0
                if nargout == 1 
                    E = ladel_mex('return');
                    varargout{1} = E;
                elseif nargout == 2
                    [L, D] = ladel_mex('return');
                    varargout{1} = L;
                    varargout{2} = D;
                elseif nargout == 3
                    [L, D, p] = ladel_mex('return');
                    varargout{1} = L;
                    varargout{2} = D;
                    varargout{3} = p;
                elseif nargout == 4
                    [L, D, p, E] = ladel_mex('return');
                    varargout{1} = L;
                    varargout{2} = D;
                    varargout{3} = p;
                    varargout{4} = E;
                else
                    error('Wrong number of output arguments for factorize_with_prior_basis');
                end
            end 
        end
        
        function varargout = row_mod(~, row, varargin)
            
            if ~(nargin ==2 || nargin ==4 || nargin == 5)
                error('Wrong number of input arguments for row_mod.\n Use .row_mod(row) to delete a row (with index row) or .row_mod(row, w, diag_elem) to add a row w and with on the diagonal diag_elem.\n');
            end
            if nargin == 2
                ladel_mex('rowmod', int64(row), length(row));
            elseif nargin == 4
                w = varargin{1};
                diag_elem = varargin{2};
                ladel_mex('rowmod', int64(row), length(row), w, diag_elem);  
            else
                w = varargin{1};
                diag_elem = varargin{2};
                beta = varargin{3};
                ladel_mex('rowmod', int64(row), length(row), w, diag_elem, beta);
            end
            
            if nargout > 0
                if nargout == 1
                    E = ladel_mex('return');
                    varargout{1} = E;
                elseif nargout == 2
                    [L, D] = ladel_mex('return');
                    varargout{1} = L;
                    varargout{2} = D;
                elseif nargout == 3
                    [L, D, p] = ladel_mex('return');
                    varargout{1} = L;
                    varargout{2} = D;
                    varargout{3} = p;
                elseif nargout == 4
                    [L, D, p, E] = ladel_mex('return');
                    varargout{1} = L;
                    varargout{2} = D;
                    varargout{3} = p;
                    varargout{4} = E;
                else
                    error('Wrong number of output arguments for factorized_advanced');
                end
            end
        end
    end
    
end

