% This codes is just an example of how the source is written to a text file
% it is not required by the OpticalNamedNode or the OptNodeDlList code.


field_names = fieldnames(Elements{1}); %Field names of structure, or public fields of object
param_names = strvcat(field_names);      %parameter labels (names)
for pnum=1:length(param_names);       %loops over all parameters of current element
    param_label = param_names(pnum,:)
    param_label=strtrim(param_label);
    switch param_label
        case 'BaseUnits'
            if exist('BaseUnits') ;
                temp_param=eval(param_names(pnum,:)); %param value
            elseif ~exist('BaseUnits');
                temp_param=eval(['Elements{1}.' param_names(pnum,:)]);
            end
            fprintf(fname_id, 'BaseUnits=  %s\n',temp_param);
            clear temp_param;

        case 'WaveUnits'
            if exist('WaveUnits');
                temp_param=eval(param_names(pnum,:)); %param value
            elseif ~exist('WaveUnits');
                temp_param=eval(['Elements{1}.' param_names(pnum,:)]);
            end
            fprintf(fname_id, 'WaveUnits=  %s\n',temp_param);
            clear temp_param;

        case 'ChfRayDir'
            if exist('ChfRayDir');
                temp_param=eval(param_names(pnum,:));
                [row col]=size(temp_param);
            elseif ~exist('ChfRayDir');
                temp_param=eval(['Elements{1}.' param_names(pnum,:)]);
                [row col]=size(temp_param);
            end;
            fprintf(fname_id, 'ChfRayDir=');
            for ii=1:col;
                if abs(temp_param(1,ii))<1 && abs(temp_param(1,ii))~=0;
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.11e'),'e-0','D-'));
                else
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.11e'),'e+0','D+'));
                end;
            end
            clear temp_param;

        case 'ChfRayPos'
            if exist('ChfRayPos');
                temp_param=eval(param_names(pnum,:));
                [row col]=size(temp_param);
            elseif ~exist('ChfRayPos');
                temp_param=eval(['Elements{1}.' param_names(pnum,:)]);
                [row col]=size(temp_param);
            end
            fprintf(fname_id, '\nChfRayPos=');
            for ii=1:col;
                if abs(temp_param(1,ii))<1 && abs(temp_param(1,ii))~=0;
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.11e'),'e-0','D-'));
                else
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.11e'),'e+0','D+'));
                end;
            end
            clear temp_param;

        case 'zSource'
            if exist('zSource');
                temp_param=eval(param_names(pnum,:));
            elseif ~exist('zSource');
                temp_param=eval(['Elements{1}.' param_names(pnum,:)]);
            end
            fprintf(fname_id, '\n  zSource=')
            if abs(temp_param)<1 && abs(temp_param~=0);
                fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.11e'),'e-0','D-'));
            else
                fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.11e'),'e+0','D+'));
            end;
            clear temp_param;

        case 'IndRef'
            if exist('IndRef');
                temp_param=eval(param_names(pnum,:));
            elseif ~exist('IndRef');
                temp_param=eval(['Elements{1}.' param_names(pnum,:)]);
            end
            fprintf(fname_id, '\n   IndRef=')
            if abs(temp_param)<1 && abs(temp_param~=0);
                fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.11e'),'e-0','D-'));
            else
                fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.11e'),'e+0','D+'));
            end;
            clear temp_param;
        case 'Extinc'
            if exist('Extinc');
                temp_param=eval(param_names(pnum,:));
            elseif ~exist('Extinc');
                temp_param=eval(['Elements{1}.' param_names(pnum,:)]);
            end
            fprintf(fname_id, '\n   Extinc=')
            if abs(temp_param)<1 && abs(temp_param~=0);
                fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.11e'),'e-0','D-'));
            else
                fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.11e'),'e+0','D+'));
            end;
            clear temp_param;

        case 'Wavelen'
            if exist('Wavelen');
                temp_param=eval(param_names(pnum,:));
            elseif ~exist('Wavelen');
                temp_param=eval(['Elements{1}.' param_names(pnum,:)]);
            end
            fprintf(fname_id, '\n  Wavelen=')
            if abs(temp_param)<1 && abs(temp_param~=0);
                fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.11e'),'e-0','D-'));
            else
                fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.11e'),'e+0','D+'));
            end;
            clear temp_param;
        case 'Flux'
            if exist('Flux');
                temp_param=eval(param_names(pnum,:));
            elseif ~exist('Flux');
                temp_param=eval(['Elements{1}.' param_names(pnum,:)]);
            end
            fprintf(fname_id, '\n     Flux=')
            if abs(temp_param)<1 && abs(temp_param~=0);
                fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.11e'),'e-0','D-'));
            else
                fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.11e'),'e+0','D+'));
            end;
            clear temp_param;

        case 'GridType'
            if exist('GridType');
                temp_param=eval(param_names(pnum,:)); %param value
            elseif ~exist('GridType');
                temp_param=eval(['Elements{1}.' param_names(pnum,:)]);
            end
            fprintf(fname_id, '\n GridType=  %s',temp_param);
            clear temp_param;

        case 'Aperture'
            if exist('Aperture');
                temp_param=eval(param_names(pnum,:));
            elseif ~exist('Aperture');
                temp_param=eval(['Elements{1}.' param_names(pnum,:)]);
            end
            fprintf(fname_id, '\n Aperture=')
            if abs(temp_param)<1 && abs(temp_param~=0);
                fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.11e'),'e-0','D-'));
            else
                fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.11e'),'e+0','D+'));
            end;
            clear temp_param;

        case 'Obscratn'
            if exist('Obscratn');
                temp_param=eval(param_names(pnum,:));
            elseif ~exist('Obscratn');
                temp_param=eval(['Elements{1}.' param_names(pnum,:)]);
            end
            fprintf(fname_id, '\n Obscratn=')
            if abs(temp_param)<1 && abs(temp_param~=0);
                fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.11e'),'e-0','D-'));
            else
                fprintf(fname_id, '  %s',strrep(num2str(temp_param,'%1.11e'),'e+0','D+'));
            end;
            clear temp_param;

        case 'nGridpts'
            if exist('nGridpts');
                temp_param=eval(param_names(pnum,:));
            elseif ~exist('nGridpts');
                temp_param=eval(['Elements{1}.' param_names(pnum,:)]);
            end
            fprintf(fname_id, '\n nGridpts=')
            fprintf(fname_id, '  %1.0f',temp_param);
            clear temp_param;

        case 'xGrid'
            if exist('xGrid');
                temp_param=eval(param_names(pnum,:));
                [row col]=size(temp_param);
            elseif ~exist('xGrid');
                temp_param=eval(['Elements{1}.' param_names(pnum,:)]);
                [row col]=size(temp_param);
            end;
            fprintf(fname_id, '\n    xGrid=');
            for ii=1:col;
                if abs(temp_param(1,ii))<1 && abs(temp_param(1,ii))~=0;
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.11e'),'e-0','D-'));
                else
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.11e'),'e+0','D+'));
                end;
            end
            clear temp_param;
        case 'yGrid'
            if exist('yGrid');
                temp_param=eval(param_names(pnum,:));
                [row col]=size(temp_param);
            elseif ~exist('yGrid');
                temp_param=eval(['Elements{1}.' param_names(pnum,:)]);
                [row col]=size(temp_param);
            end;
            fprintf(fname_id, '\n    yGrid=');
            for ii=1:col;
                if abs(temp_param(1,ii))<1 && abs(temp_param(1,ii))~=0;
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.11e'),'e-0','D-'));
                else
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.11e'),'e+0','D+'));
                end;
            end
            clear temp_param;

        case 'SegXGrid'
            if exist('SegXGrid');
                temp_param=eval(param_names(pnum,:));
                [row col]=size(temp_param);
            elseif ~exist('SegXGrid');
                temp_param=eval(['Elements{1}.' param_names(pnum,:)]);
                [row col]=size(temp_param);
            end;
            fprintf(fname_id, '\n SegXGrid=');
            for ii=1:col;
                if abs(temp_param(1,ii))<1 && abs(temp_param(1,ii))~=0;
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.11e'),'e-0','D-'));
                else
                    fprintf(fname_id, '  %s',strrep(num2str(temp_param(1,ii),'%1.11e'),'e+0','D+'));
                end;
            end
            clear temp_param;
        case 'nElt'
            if exist('nElt');
                temp_param=eval(param_names(pnum,:));
            elseif ~exist('nElt');
                temp_param=eval(['Elements{1}.' param_names(pnum,:)]);
            end
            fprintf(fname_id, '\n     nElt=')
            fprintf(fname_id, '  %1.0f\n',temp_param);
            clear temp_param;
    end
end
    