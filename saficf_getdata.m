function [xdat,ydat,edat,axeshandle] = saficf_getdata()

% Reference: fitgetdata.m from the custom fitting routines used on ID20, by SPC, SBW.

disp('Click on figure with spectra');
waitforbuttonpress;
handle = gco;

% Checks that the current object is a line, errorbar or something else.
if strcmp(get(handle,'Type'),'line')|strcmp(get(handle,'Type'),'hggroup');  % hggroup == errorbar
  linehandle = handle;
else
  grphdat = [findobj(handle,'Type','line');findobj(handle,'Type','hggroup')];
  if isempty(lines)
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
