% Purpose: This m-files is used to calculate rigidbody sensitivities, this
%          includes entire instrument rigidbody which is equivalent to
%          source rigidbody sensitivities. All sensitivities are in global
%          coordinates. The sensitivities may or may not include the stop,
%          and or the FSM; this is determined by user in the initialization
%          file. For more information contact Luis Marchen (818)393-7235 or
%          lmarchen@jpl.nasa.gov
%
% Record:
%       Date            Author               Description
%  ===============  ==================     ===============
%     12/13/2008      Luis F. Marchen       Orignial code
%
% Define Variables and Constants:
% evalopd      ---Element at which OPD is evaluated

function [dxdu d0du dFdU]=makeRBSens(param)
%==========================================================================
% Define size of full sensitivity (matrix dwdu + d0du)
%========================================================================== 
[row1 col1]=size(param.chfrayPos); %computes size of chfrayPos for chiefray intersection
Sx=zeros(row1,row1*6); % row-1 => row on February 18th 2010
Sy=zeros(row1,row1*6);
Sz=zeros(row1,row1*6);
%==========================================================================
% Calculate gain matrix to control FSM, to be used only if param.ifFSM = 'On'
%==========================================================================
if strcmpi(param.ifFSM,'On');
[Gm dif pertpos]=makeGM(param); %generate FSM gain matrix
end

%==========================================================================
% Beamwalk Sensitivity Matrix Calculation
%==========================================================================
[row2 col2]=size(param.rbSrf); %computes size of chfrayPos for chiefray intersection
for jj=1:row2;
    for ii=1:6;

        CARG1=' ';
        CARG2=' ';
        DARG=zeros(9,1);
        IARG=zeros(9,1);
        LARG=zeros(9,1);
        RARG=zeros(9,1);
        zip=0;
        zip1=zeros(6,1);
        zip2=zeros(3,1);
        chfRayRef = [0.000000000e+00  0.000000000e+00  1.000000000e+00 ...
            0.000000000e+00  0.0  0.0];

%--------------------------------------------------------------------------
% Loading MACOS prescription
%--------------------------------------------------------------------------
        command='old';
        CARG1=param.Rx;
        callsmacos;
%--------------------------------------------------------------------------
% Enforcing system stop
%--------------------------------------------------------------------------
        if strcmpi(param.ifSTOP,'On');
            command='stop';
            CARG1='elt';
            IARG(1)=param.iSTOP;
            DARG(1,1)=0.0;
            DARG(2,1)=0.0;
            callsmacos;
        end
%--------------------------------------------------------------------------
% Defining iterator to save data in matrix 1:ielt by dxji
%--------------------------------------------------------------------------
        J=(jj-1)*6+ii;
%--------------------------------------------------------------------------
% Computing nominal chiefray intersection
%--------------------------------------------------------------------------
        for k=1:row1;%param.iEvalSrf; changed this to row on Tuesday February 16th 2010 (when param.iEvalSrf is smaller than param.rbSrf)
            command='spot';
            IARG=zeros(9,1);
            IARG(1)=param.chfrayPos(k,:); %raypos on elements of interest;
            CARG1='tout';
            callsmacos;
            raypos_nom(:,k)=rayPos(1:3,1); %nominal chiefray position on all elements
        end
        
        %chiefray intersection at image plane for FSM pointing
            clear rayPos
            command='spot';
            IARG=zeros(9,1);
            IARG(1)=param.iEvalSrf; %raypos on elements of interest;
            CARG1='tout';
            callsmacos;
            raypos_FSMnom=rayPos(1:3,1); %nominal chiefray position on all elements

%--------------------------------------------------------------------------
% Plot nominal OPD
%--------------------------------------------------------------------------
        if strcmpi(param.ifPLOT,'On')
            command='opd';
            IARG=zeros(9,1);
            IARG(1)=param.iOPD; %raypos on elements of interest;
            CARG1='tout';
            callsmacos;
            OPDnom=OPD;%OPD(1:24,1:24);
            
            if strcmpi(param.units,'mm')
               figure(1), graph(OPDnom./1e-6), axis image%, niceaxes,
               title('Nominal WF'),
               xlabel(sprintf('RMS = %0.4g,    PV = %0.4g nm', stat2d(OPDnom./1e-6))), drawnow
            elseif strcmpi(param.units,'m');
                   figure(1), graph(OPDnom./1e-9), axis image%, niceaxes,
                   title('Nominal WF'),
                   xlabel(sprintf('RMS = %0.4g,    PV = %0.4g nm', stat2d(OPDnom./1e-9))), drawnow
            end
            clear OPD

        end
%--------------------------------------------------------------------------
% Perturb elements in param.rbSrf, and make sure that the units of
% translations are consistent with the units used in MACOS prescription
%--------------------------------------------------------------------------
        if isempty(find(param.rbSrf(jj,:))); %to perturb source as well
            kk=1;
        elseif ~isempty(find(param.rbSrf(jj,:)));
            kk=find(param.rbSrf(jj,:))';
        end

        for u=kk
            command='perturb';
            IARG=zeros(9,1) ;
            IARG(u)=param.rbSrf(jj,u); %element/s to be perturbed
            CARG1='global';
            DARG=zeros(9,1);
            DARG(ii)=param.pvec(ii);
            callsmacos;
        end
%--------------------------------------------------------------------------
% Enforcing system stop at element 14 after perturbation
%--------------------------------------------------------------------------
        if strcmpi(param.ifSTOP,'On');;
            command='stop';
            CARG1='elt';
            IARG(1)=param.iSTOP;
            DARG(1,1)=0.0;
            DARG(2,1)=0.0;
            callsmacos;
        end
%--------------------------------------------------------------------------
% Plot perturbed OPD to see what it looks like
%--------------------------------------------------------------------------
        if strcmpi(param.ifPLOT,'On')
            command='opd';
            IARG=zeros(9,1);
            IARG(1)=param.iOPD; %raypos on elements of interest;
            CARG1='tout';
            callsmacos;
            OPD=OPD;%OPD(1:24,1:24);
            if strcmpi(param.units,'mm')
               figure(2), graph((OPD-OPDnom)./1e-6), axis image%, niceaxes,
               title(['Residual RMS OPD (Perturbed iElt=' num2str(param.rbSrf(jj,1)) ',DOF=' num2str(ii) ')']),
               xlabel(sprintf('RMS = %0.4g,    PV = %0.4g nm', stat2d((OPD-OPDnom)./1e-6))), drawnow
            elseif strcmpi(param.units,'m');
               figure(2), graph((OPD-OPDnom)./1e-9), axis image%, niceaxes,
               title(['Residual RMS OPD (Perturbed iElt=' num2str(param.rbSrf(jj,1)) ',DOF=' num2str(ii) ')']),
               xlabel(sprintf('RMS = %0.4g,    PV = %0.4g nm', stat2d((OPD-OPDnom)./1e-9))), drawnow
            end
        end
%==========================================================================
% Command to calibrate after perturbation, and to calculate actuations
%==========================================================================
        if strcmpi(param.ifFSM,'On');
            Dpf=0;
            Dp=0;
            clear raypos_fsm
            
            for uu=1:4; %iterates fsm controls specified times
                clear raypos_fsm
                command='spot';
                IARG=zeros(9,1);
                IARG(1)=param.iEvalSrf;
                CARG1='tout';
                callsmacos;
                raypos_fsm=rayPos(1:3,1);

                dp=raypos_fsm-raypos_FSMnom % made change March 15th 2010
                Dp=(Gm*dp);
                Dpf=Dpf+(Gm*dp);
                for yy=1:3;
                    command='perturb';
                    IARG=zeros(9,1) ;
                    IARG(1)=param.iFSM; %element to be perturbed
                    CARG1='global';
                    DARG=zeros(9,1);
                    DARG(yy)=-Dp(yy);
                    callsmacos;
                end
            end
%--------------------------------------------------------------------------
% Determine how accurate the FSM steers the beam, if the difference is
% greater than 1e-9 it displays a warning
%--------------------------------------------------------------------------
            if sqrt(sum(dp.^2))>1e-9
                warning(['FSM Failed to Repoint Chiefray at Element' param.iEvalSrf])
            end
        end;
%--------------------------------------------------------------------------
% Enforcing system stop at element 14 after perturbation
%--------------------------------------------------------------------------
        if strcmpi(param.ifSTOP,'On');;
            command='stop';
            CARG1='elt';
            IARG(1)=param.iSTOP;
            DARG(1,1)=0.0;
            DARG(2,1)=0.0;
            callsmacos;
        end        
%--------------------------------------------------------------------------
% Spot command to trace chiefray
%--------------------------------------------------------------------------
        clear rayPos
        for k=1:row1;%param.iEvalSrf; changed this to row on Tuesday February 16th 2010 (when param.iEvalSrf is smaller than param.rbSrf)
            command='spot';
            IARG=zeros(9,1);
            IARG(1)=param.chfrayPos(k,:); %raypos on elements of interest;
            CARG1='tout';
            callsmacos;
            raypos_pert(:,k)=rayPos(1:3,1);
        end
%--------------------------------------------------------------------------
% If units in MACOS file are mm we convert deltas from mm to meters, if the
% units are meters we do not perform any convertions. The point is to get
% the final sensitivity in units of Rad/m and m/m.
%--------------------------------------------------------------------------
        if strcmpi(param.units,'mm'); %if units are mm a unit convertion is made
            raypos_nom_m=raypos_nom./1000;   %convert from mm to meters by dividing by 1000
            raypos_pert_m=raypos_pert./1000; %convert from mm to meters by dividing by 1000
            pvec=param.pvec;
            pvec(4:6)=pvec(4:6)./1000; %converts translations from mm to m
        elseif strcmpi(param.units,'m');
            raypos_nom_m=raypos_nom; 
            raypos_pert_m=raypos_pert;
            pvec=param.pvec;
        end

        d=(raypos_pert_m-raypos_nom_m)./pvec(ii);
        if jj>1 %changed March 11th 2010, need to subtract 1 from right element 0 element does not count
        if ii>3;
            d(ii-3,jj-1)=d(ii-3,jj-1)-1;
        end
        end

        Sx(:,J)=d(1,:)';
        Sy(:,J)=d(2,:)';
        Sz(:,J)=d(3,:)';

    end
end

dFdU=[Sx; Sy;Sz];
[rr cc]=size(dFdU);

d0du=dFdU(:,1:2); %This is the source sensitivity matrix
dxdu=dFdU(:,7:cc);%This is the rigidbody sensitivity matrix

if strcmpi(param.ifSAVE,'On') 
   if strcmpi(param.cType,'FSM')
      save('JitterBWFSMC_dxdu.xls','dxdu','-ASCII')
      save('JitterRBFSMC_d0du.xls','d0du','-ASCII')
   elseif strcmpi(param.cType,'SM')
       save('ThermalBWSMC_dxdu.xls','dxdu','-ASCII')
       save('ThermalRBSMC_d0du.xls','d0du','-ASCII')
   else
       save('JitterBWNC_dxdu.xls','dxdu','-ASCII')
       save('JitterRBNC_d0du.xls','d0du','-ASCII')
   end 
end

%Save element names in a text file for beam walk
write_fname='bwEltList.txt'
macos_fname=[param.Rx '.in']
[OpticalData,numElt,List_EltNames]= ReadRxFile(macos_fname) 
writeRxEltNames(write_fname,OpticalData,param.bwElt);%this writes element names in a text file

writeKcElt_fname='KcElt_file.txt';
writeKcEltValues(writeKcElt_fname,OpticalData,param.bwElt) %this writes the conic constant for each element
return
