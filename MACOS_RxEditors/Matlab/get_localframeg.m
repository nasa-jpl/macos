function e = get_localframeg(v,iaxes)
% function c = get_localframeg(v,iaxes)
% given 
%  v      a set of (column) vectors (ordered)
%  iaxes  axis indices corresponding to v 
% determine 
%  c a transformation matrix 
% where v(:,1) -> e(iaxes(1)) etc
% priority is given to the first vector - the rest are orthogonalized to
% the set previous.
% 
% only for 3d vectors.
%
% 9/3/2005 sws
% 10/20/2008 sws improve calculations for close vectors

if nargin < 2
  error('[get_localframeg] need 2 arguments')
end
nv = size(v);
ni = length(iaxes);
if ni ~= nv(2)
  error('[get_localframeg] length error')
end

e=zeros(3,3);
for ii=1:ni
  j = iaxes(ii);
  % add v - e'* e*v 
  x = v(:,ii) - e'*(e*v(:,ii));
  nx = norm(x);
  %if 1 == 1+nx
  if 1 == 1 + nx/2
    disp(['warning: ignoring case ' num2str(ii) ' iaxis ' num2str(j)])
    iaxes(ii) = 0;
  else
    if nx < 1
      % numbers are small, try again
      % - repeating always seems to help
      x2 = x/nx;
      x = x2 - e'*(e*x2);
      nx = norm(x);
    end
    e(j,:) = x'/nx;
  end
end
%e

% now get any remaining vectors
ne = sum( (e.*e)' );
irest = find(ne < 0.99);
nrest = length(irest);

e0 = eye(3);
if nrest > 2
  e = e0;
  return;
end
if nrest < 1
  kk = iaxes(3);
end
if nrest == 1
  kk = irest;
end
if  nrest == 2
  for ii = irest(1)
    x = e0 - e'*e*e0;
    nx = sum( x.*x );
    [xx,i] = max(nx);
    e(ii,:) = x(:,i)'/sqrt(nx(i));
  end
  kk = irest(2);
end
%kk
% now get kk from cross product
in = [1 2 3 1 2 3];
j = max(find( (in== kk) ));
%e(kk,:) = sws_cross( e(in(j-2),:), e(in(j-1),:));
e(kk,:) = cross( e(in(j-2),:), e(in(j-1),:));
  
% end
