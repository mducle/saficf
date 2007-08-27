function [right,datpos] = saficf_buttonpress();
% Returns logical 1 if a right button is pressed, otherwise returns position of nearest
% data point to mouse position.
%
% Syntax: [right,data] = saficf_buttonpress()
%
% Based on mf_btprs.m from the mfit4 sources.
% By Duc Le 2007 - duc.le@ucl.ac.uk

%waitforbuttonpress;
[xdat,ydat,edat,handle] = saficf_getdata(gca);
mouse_pos = get(gca,'CurrentPoint');
[diff ind_diff] = min((xdat-mouse_pos(1)).^2);  % ind_diff is the index of closest point

if strcmp(get(gcf,'SelectionType'),'alt')       % alt == either shift or right click.
  right = 1;
  datpos = [0 0 0];
else
  right = 0;
  datpos = [xdat(ind_diff) ydat(ind_diff) edat(ind_diff)];
end
