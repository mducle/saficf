function [xdat,ydat,edat,axeshandle] = saficf_getdata(handle)

% Reference: fitgetdata.m from the custom fitting routines used on ID20, by SPC, SBW.

if ~exist('handle')
  disp('Click on figure with spectra');
  waitforbuttonpress;
  handle = gcf;
end

% Checks that the current figure is an mslice cut
if strcmp(get(get(0,'CurrentFigure'),'Tag'),'plot_cut')
  grphdat = findobj(get(0,'CurrentFigure'),'Type','line');
  if isempty(grphdat); disp(['No line or surface objects in: ',handle]); return; end		
  if(length(grphdat)>2)
    ask = questdlg(sprintf(['Warning: More than one dataset in Mslice window.\n Is this correct,' ...
                            'and do you want to fit all datasets?']),'SAFiCF Warning','No');
    if strcmp(ask,'No')
      xdat = get(grphdat(1),'XData'); ydat = get(grphdat(1),'YData'); 
      uldat = get(grphdat(2),'YData'); udat = []; ldat = [];
      for iEdat = 1:9:length(uldat)
        udat = [udat; uldat(iEdat)];
        ldat = [ldat; uldat(iEdat+1)];
      end
      edat = (udat-ldat)./2;
      axeshandle = gca;
    elseif strcmp(ask,'Cancel')
      xdat=[]; ydat=[]; edat=[]; axeshandle=gca;
    else
      for idataset = 1:length(grphdat/2)
        xdat{idataset} = get(grphdat(2*idataset-1),'XData'); ydat{idataset} = get(grphdat(2*idataset-1),'YData');
        uldat = get(grphdat(2*idataset),'YData'); udat = []; ldat = [];
        for iEdat = 1:9:length(uldat) udat = [udat; uldat(iEdat)]; ldat = [ldat; uldat(iEdat+1)]; end
        edat{idataset} = (udat-ldat)./2;
        axeshandle = gca;
      end
    end
  else
    xdat = get(grphdat(1),'XData'); ydat = get(grphdat(1),'YData'); 
    uldat = get(grphdat(2),'YData'); udat = []; ldat = [];
    for iEdat = 1:9:length(uldat)
      udat = [udat; uldat(iEdat)];
      ldat = [ldat; uldat(iEdat+1)];
    end
    edat = (udat-ldat)./2;
    axeshandle = gca;
  end
else
  % Checks that the current object is a line, errorbar or something else.
  if strcmp(get(handle,'Type'),'line')|strcmp(get(handle,'Type'),'hggroup');  % hggroup == errorbar
    linehandle = handle;
  else
    grphdat = [findobj(handle,'Type','hggroup');findobj(handle,'Type','line')];
    if isempty(grphdat)
      disp(['No line or surface objects in: ',handle]);
      return
    end		
    % Assigns handle to first line/errorbar child object with a tag.
    linehandle = grphdat(1);
    for ind_gd = 1:length(grphdat);
      if isempty(get(grphdat(ind_gd),'Tag'))
        linehdl = grphdat(ind_gd);
        break
      end
    end;     
  end;

  xdat = get(linehandle,'XData');
  ydat = get(linehandle,'YData');

  if findobj(linehandle,'Type','hggroup')
    udat = get(linehandle,'UData');
    ldat = get(linehandle,'LData');
    edat = (udat+ldat)./2;
  else
    edat = zeros(size(xdat));
  end

  axeshandle = gca;
end
