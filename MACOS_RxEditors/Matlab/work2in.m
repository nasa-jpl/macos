% work2in	script file which translates imos optical
%		prescription in the workspace into a 
%		MACOS .in file
%	
%       The data is appended to the file fname where fname exists in the
% workspace. Note: fname must include the suffix (.in) if it is to be
% used by MACOS. If fname is not specified, then the optical prescription
% is written to (comp.in). If fname_overwrite exists and is 0, then fname
% will not be overwritten, else it will be deleted if it exists. 
%	If Tout is specified, then output coordinates
% are included in the .in file.  Otherwise, output coordinates
% are not included.
%	If nECoord is specified for each element, then for any
% element with nECoord>0, the proper columns of TElt (which must
% be specified) are written into the .in file.  TElt should contain
% the element coordinates for all of the elements.  TElt should be
% 6 rows by n columns, where n is sum(nECoord>0).
%
% Extra variables:
%  author - if it exists, this string is added to the file
%  notes  - if it exists, this string or set of stings is added
%
% see Osys_init.m for list of EltTypes and mappings to (Eltment,Surface) 
%
% Assumptions:
%  for the source:
%    zSource = 1.0D22
%    IndRef = 1
%    Extinc = 0
%    Flux = 1.0
%  for other elements:
%    Extinc = 1.0D22
% 
% Copyright 1992.  National Aeronautics and Space Administration,
% all rights reserved.
% Copyright 1993-2000.  California Institute of Technology, National
% Aeronautics and Space Administration.  All rights reserved.  US Government
% Sponsorship acknowledged.
% 
% THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
% INCLUDING BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTIBILITY
% AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL
% CALTECH BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR
% ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
% WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTUOUS
% ACTION, ARISING OUT OF OR IN CONNECTION WITH THE ACCESS, USE OF
% PERFORMANCE OF THIS SOFTWARE.
% 

% History
%  17Dec97  lneedels changed max(size( to length for matlab5
%  02Jun98  lneedels fixed bug in warning about filename already existing
% 9/04/2003 sws fixing up... folding in many changes by others
%                            fopen/fclose for new matlab
%                            baseunits, waveunits...
% 10/1/2003  sws added note variable
% 12/16/2003 sws don't define defaults for vectors; scalars are ok
% 4/13/2004  sws minor fix to GridType, ApType
% 4/16/2004  sws fixed bug in EltType(5); added some more
% 1/6/2005   sws create fElt,eElt if they don't exist
% 12/26/2005 sws added assumtions (doc) 
%  2/17/2006 sws fix for note, reorder 
%  9/16/2006 sws fix ObsType to allow all elements
%  9/17/2006 sws use Extinc, zSource, IndRef_Source  if available
%                    zernike added
%                    bug in zElt 
%                    ObsVec -- matrix??
% 5/15/2007          ObsVec now in work as a cell vector of vectors
% 6/22/2007 sws handle multiple zernikes        
% 6/29/2007     fix xObs handling
% 7/2/2007      ApVec, general element types
% 9/17/2007 sws update for new macos -- can handle long input lines!


% for some things, uses imos_v2s for formatting

work2workA;

if exist('fname')~=1 % fname does not exist in the workspace
  fname='comp.in';
end

% sws don't do this or it could cause problems when another system is
% built, as then nECoord will exist but might be the wrong size
% if exist('nECoord')~=1	%nECoord does not exist in the workspace
%  nECoord=-6*ones(nElt,1);
% end

if exist('Wavelen') ~= 1
  Wavelen=0.55e-6;
  WaveUnits = 'm';
end

if exist('Tout') ~= 1
  Tout=eye(7);
end

[tempn,tempm]=size(Tout);
if tempm ~= 7
  error('Tout must be n by 7!!')
end

% global that will allow us to overwrite the file or not
if exist('fname_overwrite') ~= 1
  fname_overwrite = 1;
end

if exist(fname) == 2;	%fname already exists as file
  if fname_overwrite
    disp(['deleting ' fname])
    delete(fname);
  else
    warning([fname ' already exists!!!'])
    return
  end
end

if exist('nECoord')==1	
  tempm=find(nECoord>0);
  if length(tempm)>0
    temptot=sum(nECoord(tempm));
    if exist('TElt')~=1	%TElt does not exist in the workspace
      error('No element coordinates are specified!!!');
    end
    [tempnn,tempmm]=size(TElt);
    if tempmm ~= temptot
      disp(['TElt should have first dimension ' num2str(temptot) ])
      disp(['but has dimension ' num2str(tempmm)])
      error('TElt is the wrong size!!');
    end
    tempm=1;	%to keep track of where in TElt you are
  end
end % of nECoord check  
  
disp(['[work2in] writing ' fname])
fname_id = fopen(fname,'w');

%credits
fprintf(fname_id,'%% %s\n',fname);
if exist('author') == 1
  fprintf(fname_id,'%% Author: %s\n', author);
else
  fprintf(fname_id,'%% Author:  (fill in)\n');
end  
if exist('note') == 1
  % allow multiline notes
  [m,n]= size(note);
  for i=1:m
    fprintf(fname_id,'%% %s\n', note(i,:));    
  end
end
fprintf(fname_id,'%% %s\n\n',date);

% input plane
fprintf(fname_id,'BaseUnits=  %s\n',BaseUnits);
fprintf(fname_id,'WaveUnits=  %s\n',WaveUnits);
fprintf(fname_id,'ChfRayDir=%s\n',imos_v2s(ChfRayDir));
fprintf(fname_id,'ChfRayPos=%s\n',imos_v2s(ChfRayPos));
if exist('zSource')
  fprintf(fname_id,' zSource=   %s\n', imos_v2s(zSource));
else
  fprintf(fname_id,' zSource=   1.000000000D+22\n');
end
if exist('IndRef_Source')
  fprintf(fname_id,'  IndRef=   %s\n', imos_v2s(IndRef_Source));
else
  fprintf(fname_id,'  IndRef=   1.000000000D+00\n');
end
if exist('Extinc_Source')
  % use if it exists
  fprintf(fname_id,'  Extinc=   %1.14e\n',Extinc_Source);
else
  fprintf(fname_id,'  Extinc=   0.000000000D+00\n');
end
fprintf(fname_id,' Wavelen=  %s\n',imos_v2s(Wavelen));
fprintf(fname_id,'    Flux=   1.000000000D+00\n');
if GridType==1
  fprintf(fname_id,'GridType=   Circular\n');
else
  %fprintf(fname_id,'GridType=   Rectangular\n');% sws bug?
  fprintf(fname_id,'GridType=   %d\n', GridType);
end
fprintf(fname_id,'Aperture= %s\n',imos_v2s(Aperture));
fprintf(fname_id,'Obscratn= %s\n',imos_v2s(Obscratn));
fprintf(fname_id,'nGridpts=   %1.0f\n',nGridpts);

% check orthogonality
dotp = ChfRayDir*xGrid';
if abs(dotp) > 1e-6
  disp(['warning; ChfRayDir and xGrid not orgthogonal: ' ...
	num2str(acos(dotp)*180/pi) ' deg (fixing)'])
end
xGrid = xGrid - dotp*ChfRayDir;
if norm(xGrid) < eps
  error('xgrid paralell to ChfRayDir')
end
xGrid = xGrid/norm(xGrid);
fprintf(fname_id,'   xGrid= %s\n',imos_v2s(xGrid));
if ~exist('yGrid')
    yGrid=-ChfRayDir*crmat(xGrid');	%from trace.m
    yGrid = yGrid/norm(yGrid);
end
% check y orthogonality
dotp = yGrid*[ChfRayDir' xGrid'];
if any(abs(dotp) > 1e-6)
  disp(['warning; [ChfRayDir xGrid] not orgthogonal to yGrid: ' ...
	num2str(acos(dotp)*180/pi) ' deg (fixing)'])
  yGridn = (crmat(ChfRayDir)*xGrid')';
  % ensure we have the same general direction
  if yGridn*yGrid' < 0
    yGrid = -yGridn;
  else
    yGrid = yGridn;
  end
  yGrid = yGrid/norm(yGrid);
  clear yGridn;
end
clear dotp;

fprintf(fname_id,'   yGrid= %s\n',imos_v2s(yGrid));
fprintf(fname_id,'    nElt=  %3.0f\n',nElt);

% segments
% nSeg, width, gap, SegXgrid, SegCoord
flag_segments = 0;
if exist('nSeg')
  if nSeg > 0
    flag_segments = 1;
    fprintf(fname_id,'    nSeg=  %3.0f\n',nSeg);    
    fprintf(fname_id,'   width= %s\n',imos_v2s(width));
    fprintf(fname_id,'     gap= %s\n',imos_v2s(gap));
    fprintf(fname_id,'SegXgrid= %s\n',imos_v2s(SegXgrid));
    
    fprintf(fname_id,'SegCoord= %3.0f %3.0f %3.0f\n',SegCoord(1,:));
    for ii=2:nSeg
      fprintf(fname_id,'          %3.0f %3.0f %3.0f\n',SegCoord(ii,:));
    end
  end
end

if exist('ApType')
  % list of apertures
  if any(ApType>1)
    error('ApType > 1 not handled yet')
  end
  flag_ApType = ApType > 0;
  % an element may have an ApVec, length may vary
  map_ApVec = flag_ApType.*cumsum(flag_ApType); 
end

% get list of zern elements, and mapping to reduced set
% flag_zern = (EltType == 14) | (EltType == 15);
jj = imos_matchstr('Zernike',Surface);
flag_zern = zeros(size(EltType));
flag_zern(jj) = 1;
% iElt_zern = find( flag_zern );
map_zern = flag_zern.*cumsum(flag_zern);

if exist('nObs')
  % get some mappings
  xx = cumsum(nObs);
  map_obs = [xx - nObs xx];% (start-1, end) for each iElt
end
if exist('xObs')
  % present if either nObs or ApType
  x = (nObs>0) + (ApType > 0);
  map_xobs = x.*cumsum(x);% location in xObs
  nxObs = size(xObs);
  if max(map_xobs) ~= nxObs(:,1)
    disp(['need xObs for ' num2str(max(map_xobs)) ' elements'])
    disp(find(x))
    disp(['size of xObs is ' num2str(nxObs)])
    error('[work2in] xObs wrong size!');
  end
end

% elements
for iElt=1:nElt		% loop over elements
  fprintf(fname_id,'\n');
  fprintf(fname_id,'    iElt=    %1.0f\n',iElt);

  if exist('EltName')
    fprintf(fname_id,' EltName= %s\n',char(EltName(iElt)));
  else
    fprintf(fname_id,' EltName= \n');
  end

  switch(EltType(iElt))
   case(1)
    fprintf(fname_id,' Element= Reflector\n');
    fprintf(fname_id,' Surface= Conic\n');
   case(2)
    fprintf(fname_id,' Element= Reflector\n');
    fprintf(fname_id,' Surface= Flat\n');
    % bob - 2/07/00
   case(3)
    fprintf(fname_id,' Element= FocalPlane\n');
    fprintf(fname_id,' Surface= Flat\n');
    
   case(4)
    fprintf(fname_id,' Element= Reference\n');
    fprintf(fname_id,' Surface= Conic\n');
    
   case(5)% sws
    fprintf(fname_id,' Element= Segment\n');
    fprintf(fname_id,' Surface= Conic\n');
    
   case(8)% sws
    fprintf(fname_id,' Element= Refractor\n');
    fprintf(fname_id,' Surface= Conic\n');
   case(9)
    fprintf(fname_id,' Element= Obscuring\n');
    fprintf(fname_id,' Surface= Conic\n');
   case(10)% sws
    fprintf(fname_id,' Element= Return\n');
    fprintf(fname_id,' Surface= Conic\n');
   case(11)% sws
    fprintf(fname_id,' Element= Reflector\n');
    fprintf(fname_id,' Surface= Conic\n');
    
   case(13)
    fprintf(fname_id,' Element= NSReflector\n');
    fprintf(fname_id,' Surface= Conic\n');

   case(14)% sws
    fprintf(fname_id,' Element= Reflector\n');
    fprintf(fname_id,' Surface= Zernike\n');
   case(15)% sws
    fprintf(fname_id,' Element= Segment\n');
    fprintf(fname_id,' Surface= Zernike\n');
    
   otherwise
    Osys_init;
    if EltType(iElt) <= EltType_max % max old EltType
                          % leave for macos...
      fprintf(fname_id,' EltType= %1.0f\n',EltType(iElt));
    else
      % assume element*100 + surf
      jje = floor(0.5+EltType(iElt)/100);
      jjs = EltType(iElt) - jje*100;
      fprintf(fname_id,[' Element= ' EltName_table{jje}   '\n']);
      fprintf(fname_id,[' Surface= ' SrfName_table{jjs} '\n']);
    end
  end

  if (exist('KrElt') & exist('KcElt'))
    fprintf(fname_id,'   KrElt= %s\n',imos_v2s(KrElt(iElt)));
    fprintf(fname_id,'   KcElt= %s\n',imos_v2s(KcElt(iElt)));
    % create fElt and eElt
    fElt(iElt)          = -KrElt(iElt)./(1+sqrt(-KcElt(iElt)));
    eElt(iElt)          = sqrt(-KcElt(iElt));
  else
    fprintf(fname_id,'    fElt= %s\n',imos_v2s(fElt(iElt)));
    fprintf(fname_id,'    eElt= %s\n',imos_v2s(eElt(iElt)));
  end

  fprintf(fname_id,'  psiElt= %s\n',imos_v2s(psiElt(iElt,:)));
  fprintf(fname_id,'  VptElt= %s\n',imos_v2s(VptElt(iElt,:)));
  fprintf(fname_id,'  RptElt= %s\n',imos_v2s(RptElt(iElt,:)));

  if EltType(iElt) == 6	%HOE
    if exist('h1HOE') ~= 1
      error('h1HOE does not exist for the Holographic Optical Element');
    end
    if exist('h2HOE') ~= 1
      error('h2HOE does not exist for the Holographic Optical Element');
    end
    if exist('OrderHOE') ~= 1
      error('OrderHOE does not exist for the Holographic Optical Element');
    end
    if exist('WaveHOE') ~= 1
      error('WaveHOE does not exist for the Holographic Optical Element');
    end

    % fixme: format 
  fprintf(fname_id,'   h1HOE= %1.14e %1.14e %1.14e\n',...
		h1HOE(iElt,1),h1HOE(iElt,2),h1HOE(iElt,3));
  fprintf(fname_id,'   h2HOE= %1.14e %1.14e %1.14e\n',...
		h2HOE(iElt,1),h2HOE(iElt,2),h2HOE(iElt,3));
    fprintf(fname_id,'  OrderHOE= %1.14e\n',OrderHOE(iElt));
    fprintf(fname_id,'  WaveHOE= %1.14e\n',WaveHOE(iElt));
  end

  fprintf(fname_id,'  IndRef= %s\n',imos_v2s(IndRef(iElt)));
  if exist('Extinc')
    % use if it exists
    fprintf(fname_id,'  Extinc=  %s\n',imos_v2s(Extinc(iElt)));
  else
    fprintf(fname_id,'  Extinc=   1.000000000D+22\n');
  end

  % possible zernike
  % if (EltType(iElt) == 14) | (EltType(iElt) == 15)
  if flag_zern(iElt)
    %  14: (Reflector Zernike)
    % ZernType, ZernCoef, (p,x,y,z,l)Mon
    k = map_zern(iElt);
    %[iElt k    size(ZernType)]
    fprintf(fname_id,'ZernType=  %s\n',ZernType{k});
    ii = 0;
    jj = length(ZernCoef{k});
    while jj>0
      kk = min(6,jj);
      if ii==0
	fprintf(fname_id,'ZernCoef= %s\n',imos_v2s(ZernCoef{k}(ii+ (1:kk))));
      else
	fprintf(fname_id,'          %s\n',imos_v2s(ZernCoef{k}(ii+ (1:kk)))); 
      end
      ii = ii + kk;
      jj = jj-kk;
    end

    % sometimes print these -- probably a special case for sab 
    % GridFile= none                    
    % nGridMat=   99
    % GridSrfdx=   1.416291010000000D+01
    if exist('GridFile')
      fprintf(fname_id,'GridFile=  %s\n',GridFile{k});
      fprintf(fname_id,'nGridMat= %1.0f\n',nGridMat(k));
      fprintf(fname_id,'GridSrfdx= %s\n',imos_v2s(GridSrfdx(k)));
    end
    
    fprintf(fname_id,'    pMon= %s\n',imos_v2s(pMon(k,:)));
    fprintf(fname_id,'    xMon= %s\n',imos_v2s(xMon(k,:)));
    fprintf(fname_id,'    yMon= %s\n',imos_v2s(yMon(k,:)));
    fprintf(fname_id,'    zMon= %s\n',imos_v2s(zMon(k,:)));
    fprintf(fname_id,'    lMon= %s\n',imos_v2s(lMon(k,:)));
    
  end
  
  if exist('nObs')
    fprintf(fname_id,'    nObs=    %d\n',nObs(iElt));
    if nObs(iElt) > 0
      % each element can have several pairs of (type/vec) lines
      % xx=cumsum(nObs);
      % ot = xx(iElt);
      ot0 = map_obs(iElt);
      for iObs=1:nObs(iElt)
	ii = find(ObsTypeVal_table == ObsType(ot0+iObs));
	if length(ii) < 1
	  % disp(['ObsType ' num2str(ObsType(iElt,iObs))])
	  disp(['ObsType ' num2str(ObsType(ot0+iObs))])
	  error('ObsType not found')
	end
        fprintf(fname_id,[' ObsType=   ' ObsTypeName_table{ii} '\n']);
        % fprintf(fname_id,'  ObsVec= %s\n',imos_v2s(ObsVec(iElt,iObs,:)));
        % fprintf(fname_id,'  ObsVec= %s\n',imos_v2s(ObsVec(ot0+iObs,:)));
        fprintf(fname_id,'  ObsVec= %s\n',imos_v2s(ObsVec{ot0+iObs}));
      end
    end
  end
  if exist('xObs')
    % expect to have 1/element
    % either nObs>0 or ApType > 0
    %fprintf(fname_id,'    xObs= %s\n',imos_v2s(xObs(iElt,:)));
    if map_xobs(iElt) > 0
      fprintf(fname_id,'    xObs= %s\n',imos_v2s(xObs(map_xobs(iElt),:)));
    end
  end

  if exist('ApType')
    if ApType(iElt) == 0
      fprintf(fname_id,'  ApType=  None\n');
    else
      % fprintf(fname_id,'  ApType=   %d\n',ApType(iElt));
      fprintf(fname_id,['  ApType=   ' ApType_string{iElt} '\n']);
      if exist('ApVec')
        fprintf(fname_id,'   ApVec= %s\n',imos_v2s( ApVec{map_ApVec(iElt)}));
      end
    end
  else
    fprintf(fname_id,'  ApType=  None\n');
  end

  fprintf(fname_id,'    zElt= %s\n',imos_v2s(zElt(iElt)));
  
  if exist('PropType')
    if exist('PropTypeName_table')
      fprintf(fname_id,'PropType=  %s\n', ...
	  PropTypeName_table{PropType(iElt)});
    else
      switch(PropType(iElt))
      case(1)
        fprintf(fname_id,'PropType=   Geometric\n');
      otherwise
        fprintf(fname_id,'PropType=    %d\n',PropType(iElt));
      end
    end
  else
    fprintf(fname_id,'PropType=   Geometric\n');
  end

  if exist('nECoord')~=1	
    fprintf(fname_id,' nECoord=   -6.0\n');
  else
    fprintf(fname_id,' nECoord=   %1.0f\n',nECoord(iElt));
    if nECoord(iElt)>0	%print TElt
      if 1
        % new macos -- full precision
        fprintf(fname_id,'    TElt= %s',imos_v2s(TElt(1,tempm)));
        if nECoord(iElt)>1
          for tempnn=1:nECoord(iElt)-1
            fprintf(fname_id,', %s',imos_v2s(TElt(1,tempm+tempnn)));
          end
        end % end print first row of TElt
        fprintf(fname_id,'\n');
        for tempmm=2:6
          fprintf(fname_id,'          %s',imos_v2s(TElt(tempmm,tempm)));
          if nECoord(iElt)>1
            for tempnn=1:nECoord(iElt)-1
              fprintf(fname_id,', %s',imos_v2s(TElt(tempmm,tempm+tempnn)));
            end
          end
          fprintf(fname_id,'\n');
        end % end print remaining rows of TElt
      else
        % %1.14 is too much precision, 
        % need 6x < 120 chars for old macos
        fprintf(fname_id,'    TElt= %1.11e',TElt(1,tempm));
        if nECoord(iElt)>1
          for tempnn=1:nECoord(iElt)-1
            fprintf(fname_id,', %1.11e',TElt(1,tempm+tempnn));
          end
        end % end print first row of TElt

        fprintf(fname_id,'\n');
        for tempmm=2:6
          fprintf(fname_id,'          %1.11e',TElt(tempmm,tempm));
          if nECoord(iElt)>1
            for tempnn=1:nECoord(iElt)-1
              fprintf(fname_id,', %1.11e',TElt(tempmm,tempm+tempnn));
            end
          end
          fprintf(fname_id,'\n');
        end % end print remaining rows of TElt
      end % precision
      tempm=tempm+nECoord(iElt);
    end % end print TElt
  end

end % end loop over elements

% output coordinates
fprintf(fname_id,'\n');
fprintf(fname_id,'nOutCord=    %1.0f\n',tempn);

fprintf(fname_id,'    Tout= %s\n',imos_v2s(Tout(1,:)));
if tempn > 1
  for j=2:tempn
    fprintf(fname_id,'          %s\n',imos_v2s(Tout(j,:)));
  end
end

fclose(fname_id);
disp(['[work2in] done'])
% end
