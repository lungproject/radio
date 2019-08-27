function[res] = bstrap(b,f,fun,x,varargin)  %#ok
% BSTRAP    'Univariate' or paired bootstrap 
% INPUTS   : b    - number of bootstrap samples 
%            f    - bootstrap sample size, fraction of total
%            fun  - function computing desired statistics, with output
%                   placed in p*q array or structure
%            x    - n*k data matrix, with observations in rows
%            y,.. - (optional) additional vectors/matrices, referenced 
%                    by fun, in the order they appear in fun
% OUTPUTS  : res  - statistic  values in actual and bootstrap samples, 
%                   where bootstrap samples are constructed by random-
%                   ly  selecting,  with  replacement,  m = floor(f*n) 
%                   rows of x, b times. (Note  that rows of y, z, etc. 
%                   are not resampled; all variables  that need  to be 
%                   jointly  resampled should be packed in the columns
%                   of x). res  is a p*q*(b+1) array  if fun outputs a
%                   p*q array, and a 1*(b+1)  structure  array if  fun 
%                   outputs a structure, with element 1 storing  value
%                   of fun in actual sample
% EXAMPLE  : See BSTRAP_DEMO 
% SEE ALSO : BSTATS, JKNIFE; JACKKNIFE, BOOTSTRP (Statistics Toolbox)
% AUTHOR   : Dimitri Shvorob, dimitri.shvorob@vanderbilt.edu, 4/15/07

if nargin < 1
   error('Input argument "b" is undefined')
end
if nargin < 2
   error('Input argument "f" is undefined') 
end  
if nargin < 3
   error('Input argument "fun" is undefined') 
end    
if nargin < 4    
   error('Input argument "x" is undefined')
end

if ~isnumeric(b) || ~isscalar(b) || b ~= floor(b) || b <= 0
   error('Input argument "b" is invalid')
end
if ~isnumeric(f) || ~isscalar(f) || f < 0 || f > 1
   error('Input argument "f" is invalid')
end
if ~ischar(fun) 
   error('Input argument "fun" must be a string')
end
if ~isnumeric(x)
   error('Input argument "x" must be numeric')
end

n = size(x,1);
np = length(find(x(:,end)==1));
nn = length(find(x(:,end)==0));
m = floor(f*n);
mp = floor(f*np);
mn = floor(f*nn);
xp = x(find(x(:,end)==1),:);
xn = x(find(x(:,end)==0),:);
evalString = 'resi = feval(fun,xboot';
for j = 1:(nargin - 4)
    evalString = [evalString ',varargin{' num2str(j) '}'];
end
for i = 0:b
    if ~i
%        xboot = x(1:n,:);          %#ok
       xpboot = xp(1:np,:);
       xnboot = xn(1:nn,:);
       xboot = [xpboot;xnboot];
       
    else
%        iboot = ceil(n*rand(n,1));
%        xboot = x(iboot(1:m),:);   %#ok 
       ipboot = ceil(np*rand(np,1));
       xpboot = xp(ipboot(1:mp),:);   %#ok 
       inboot = ceil(nn*rand(nn,1));
       xnboot = xn(inboot(1:mn),:);   %#ok 
       xboot = [xpboot;xnboot];
    end
    try
       eval([evalString ');'])
    catch
       error('Function "fun" could not be evaluated')  
    end
    if ~i
       if isstruct(resi) 
          outputStructure = true;
          f = fieldnames(resi);
          r = length(f);
       else
          outputStructure = false; 
          res = repmat(resi,[1 1 b+1]);
       end    
    end
    if outputStructure
       for j = 1:r
           eval(['res(i+1).' f{j} ' = resi.' f{j} ';']);
       end
    else
       res(:,:,i+1) = resi;
    end
end