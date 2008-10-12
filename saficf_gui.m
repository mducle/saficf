function varargout = saficf_gui(varargin)
% saficf_gui - Creates and runs a GUI for the crystal field fitting program, SAFiCF
% 
% At present this function takes no arguments and returns no output

% Function declarations!
%
% Callbacks - All have arguments (hObject,eventdata) 
% --------------------------------------------------------------------------------
% cb_symm()        - Gets allowed CF Pars from ptgp() and sets visible Levels buttons
% cb_lvls(_mult_)  - Plots guessed levels from Allowed Levels buttons of multiplicity _mult_
% cb_erange(ind)   - Resets the y-range of the Levels axes
% cb_fitengy()     - Estimates CF parameters from energy levels
% cb_trans()       - Estimates transition intensities from energy levels
% cb_lvlbt()       - For manual input of energy levels
% cb_varentry(var) - Assigns the value in an edit box to the specified variable
% cb_close()       - Makes sure to close all figure windows created
%
% Utilities
% --------------------------------------------------------------------------------
% DisplayCFPars(hObject,B)      - Displays CF parameters in hObject's String
% SetAllowedLevels(hObject,J,B) - Sets visible allowed levels buttons
%

%------------------------------------------------------------------------------------------------------------------------------------%
% Global variables
%------------------------------------------------------------------------------------------------------------------------------------%
small = 1e-6;
symm_string = {'Triclinic','Monoclinic','Orthorhombic','Hexagonal','Cubic','------','dhcp','P6mmc','D46h','------',...
               'Ci','Cs','C1','C2','C3','C4','C6','C2h','C3h','C4h','C6h','C2v','C3v','C4v','C6v','D2','D3','D4','D6',...
               'D2h','D3h','D4h','D6h','D2d','D3d','S4','S6','T','Th','Td','O','Oh'};
numlvls = [0 0 0 0]; efit = []; J = 0; parsfit = {}; allowed = {}; T = 2; H = 0; L=0; S=0; Jdir = [];
UN = 'Units'; NM = 'normalized';

%------------------------------------------------------------------------------------------------------------------------------------%
% Creates the Figure
%------------------------------------------------------------------------------------------------------------------------------------%
FigSizePixels = [750 550 750 550]; bgcolor = get(0,'defaultUicontrolBackgroundColor'); fgcolor = get(0,'defaultUicontrolForegroundColor');
hFig        = figure('Visible','off','Name','SAFiCF','Position',[100 100 FigSizePixels(1:2)],'Color',bgcolor,'CloseRequestFcn',@cb_close);
hSpectra    = axes('Parent',hFig,'Position',[0.08,0.62,0.36,0.36],'XColor',fgcolor,'YColor',fgcolor,'Box','on');
hPhysProp   = axes('Parent',hFig,'Position',[0.08,0.20,0.36,0.36],'XColor',fgcolor,'YColor',fgcolor,'Box','on');
hLevels     = axes('Parent',hFig,'Position',[0.50,0.58,0.30,0.35],'ButtonDownFcn',{@cb_lvls,0},'XColor',fgcolor,'YColor',fgcolor,'Box','on');

hErangeTxt  = uicontrol(hFig,'Style','text','String','Erange',UN,NM,'Position',[365 520 50 15]./FigSizePixels);
hEminEdt    = uicontrol(hFig,'Style','edit','String','0',     UN,NM,'Position',[420 520 55 20]./FigSizePixels,'Callback',{@cb_erange,1});
hErngSepTxt = uicontrol(hFig,'Style','text','String','-',     UN,NM,'Position',[480 520 20 15]./FigSizePixels);
hEmaxEdt    = uicontrol(hFig,'Style','edit','String','0',     UN,NM,'Position',[505 520 55 20]./FigSizePixels,'Callback',{@cb_erange,2});
hErngMeVTxt = uicontrol(hFig,'Style','text','String','meV',   UN,NM,'Position',[560 520 30 15]./FigSizePixels);

hAllowedPan = uipanel('Parent',hFig,'Title','Allowed States','Position',[0.82 0.72 0.16 0.24]);
hSingletBut = uicontrol(hAllowedPan,'Style','pushbutton','String','Singlet',UN,NM,'Position',[0.1 0.75 0.8 0.2],'Visible','off','Callback',{@cb_lvls,1});
hDoubletBut = uicontrol(hAllowedPan,'Style','pushbutton','String','Doublet',UN,NM,'Position',[0.1 0.51 0.8 0.2],'Visible','off','Callback',{@cb_lvls,2});
hTripletBut = uicontrol(hAllowedPan,'Style','pushbutton','String','Triplet',UN,NM,'Position',[0.1 0.28 0.8 0.2],'Visible','off','Callback',{@cb_lvls,3});
hQuartetBut = uicontrol(hAllowedPan,'Style','pushbutton','String','Quartet',UN,NM,'Position',[0.1 0.05 0.8 0.2],'Visible','off','Callback',{@cb_lvls,4});
hLvlsButs   = {hSingletBut hDoubletBut hTripletBut hQuartetBut}; lvlsstr = {'Singlet','Doublet','Triplet','Quartet'};

hEstTranBut = uicontrol(hFig,'Style','pushbutton','String','Estimate Transitions',UN,NM,'Position',[620 370 120 25]./FigSizePixels,'FontSize',9,'Callback',{@cb_trans});
hEstCFParBt = uicontrol(hFig,'Style','pushbutton','String','Estimate CF Pars',    UN,NM,'Position',[620 345 120 25]./FigSizePixels,'FontSize',9,'Callback',{@cb_fitengy});

hJTxt       = uicontrol(hFig,'Style','text','String','J=',UN,NM,'Position',[360 270 50 15]./FigSizePixels);
hJEdt       = uicontrol(hFig,'Style','edit','String','0', UN,NM,'Position',[410 270 50 20]./FigSizePixels,'Callback',{@cb_symm});
hSymmTxt    = uicontrol(hFig,'Style','text','String','Symmetry',UN,NM,'Position',[470 270 100 15]./FigSizePixels);
hSymmPop    = uicontrol(hFig,'Style','popupmenu','String',symm_string,'Value',6,UN,NM,'Position',[580 270 140 25]./FigSizePixels,'Callback',{@cb_symm}); 
hFitCFBut   = uicontrol(hFig,'Style','pushbutton','String','Fit CF Pars',UN,NM,'Position',[370 200 155 50]./FigSizePixels);
hGetSpecBut = uicontrol(hFig,'Style','pushbutton','String','Get Spectra',UN,NM,'Position',[550 200 155 50]./FigSizePixels,'Callback',{@cb_getspec});

hSuscBut    = uicontrol(hFig,'Style','pushbutton','String','M(T)',  UN,NM,'Position',[30 40 70 25]./FigSizePixels,'Callback',{@cb_physprop,1});
hInvSuscBut = uicontrol(hFig,'Style','pushbutton','String','1/M(T)',UN,NM,'Position',[100 40 70 25]./FigSizePixels,'Callback',{@cb_physprop,2});
hHeatCapBut = uicontrol(hFig,'Style','pushbutton','String','Cp(T)', UN,NM,'Position',[170 40 70 25]./FigSizePixels,'Callback',{@cb_physprop,3});
hMagBut     = uicontrol(hFig,'Style','pushbutton','String','M(H)',  UN,NM,'Position',[240 40 70 25]./FigSizePixels,'Callback',{@cb_physprop,4});
hTcalcTxt   = uicontrol(hFig,'Style','text',      'String','Tcalc', UN,NM,'Position',[20 20 50 20]./FigSizePixels); 
hTcalcEdt   = uicontrol(hFig,'Style','edit',      'String','2',     UN,NM,'Position',[75 15 55 25]./FigSizePixels,'Callback',{@cb_varentry,'T'}); 
hTcalcKTxt  = uicontrol(hFig,'Style','text',      'String','K',     UN,NM,'Position',[135 20 20 20]./FigSizePixels); 
hBcalcTxt   = uicontrol(hFig,'Style','text',      'String','Bcalc', UN,NM,'Position',[175 20 50 20]./FigSizePixels); 
hBcalcEdt   = uicontrol(hFig,'Style','edit',      'String','0',     UN,NM,'Position',[230 15 55 25]./FigSizePixels,'Callback',{@cb_varentry,'H'}); 
hBcalcTTxt  = uicontrol(hFig,'Style','text',      'String','Tesla', UN,NM,'Position',[290 20 50 20]./FigSizePixels); 

hCFParsPan  = uipanel('Parent',hFig,'Title','Crystal Field Parameters','Position',[0.48 0.025 0.50 0.30]);
hCFParsTxt  = uicontrol(hCFParsPan,'Style','text',      'String','',           UN,NM,'Position',[0.02 0.02 0.7 0.95]);
hSaveParBut = uicontrol(hCFParsPan,'Style','pushbutton','String','Save Pars',  UN,NM,'Position',[0.75 0.80 0.2 0.16]);
hLoadParBut = uicontrol(hCFParsPan,'Style','pushbutton','String','Load Pars',  UN,NM,'Position',[0.75 0.64 0.2 0.16]);
hSaveSpeBut = uicontrol(hCFParsPan,'Style','pushbutton','String','Save Spec',  UN,NM,'Position',[0.75 0.42 0.2 0.16]);
hSaveLvlBut = uicontrol(hCFParsPan,'Style','pushbutton','String','Save Levels',UN,NM,'Position',[0.75 0.26 0.2 0.16]);
hSavePhyBut = uicontrol(hCFParsPan,'Style','pushbutton','String','Save Phys',  UN,NM,'Position',[0.75 0.10 0.2 0.16]);

hdlLvls     = figure('Visible','off');
hdlLvlsEdt  = uicontrol(hdlLvls,'Style','edit','String','0',UN,NM,'Position',[0.1 0.1 0.8 0.5],'Callback',{@cb_lvlbt});
hdlLvlsTxt  = uicontrol(hdlLvls,'Style','text','String','',UN,NM,'Position',[0.1 0.65 0.8 0.3]);

%------------------------------------------------------------------------------------------------------------------------------------%
%  Initialization tasks
%------------------------------------------------------------------------------------------------------------------------------------%
set(hFig,'Visible','on');

%------------------------------------------------------------------------------------------------------------------------------------%
%  Callbacks
%------------------------------------------------------------------------------------------------------------------------------------%
function cb_symm(hObject,eventdata)                         % Determines allowed energy levels from symmetry and J
  try 
    allowed = ptgp(symm_string{get(hSymmPop,'val')});
    J = str2num(get(hJEdt,'String'));
    DisplayCFPars(hCFParsTxt,allowed);
    if(J~=0)
      SetAllowedLevels(hLvlsButs,allowed);
    end
  catch
  end
end % cb_symm
%------------------------------------------------------------------------------------------------------------------------------------%
function cb_lvls(hObject,eventdata,multiplicity)
  if(multiplicity>0)
    hFp = get(hFig,'Position'); 
    set(hdlLvls,'Name',[lvlsstr{multiplicity} ' value']); 
    set(hdlLvls,'Position',[hFp(1)+hFp(3) hFp(2)+hFp(4)-150 300 150]); 
    set(hdlLvls,'Visible','on');
    figstr = sprintf('Enter energy of %s in meV.\nOr click on graph to left.',lvlsstr{multiplicity});
    set(hdlLvlsTxt,'String',figstr);
    set(hdlLvlsTxt,'Tag',num2str(multiplicity));
  elseif(strcmp(get(hdlLvls,'Visible'),'on'))
    multiplicity = str2num(get(hdlLvlsTxt,'Tag'));
    if(isempty(findobj('Tag','SAFiCF_LevelsLine')))
      for i=1:multiplicity; line('XData',[0.1 0.4],'YData',[1 1].*(i-1)*(diff(get(hObject,'YLim'))/100),'Tag','SAFiCF_LevelsLine'); end;
      efit = [zeros(1,multiplicity)];
      numlvls(multiplicity) = numlvls(multiplicity) - 1;
      set(hLvlsButs{multiplicity},'String',[lvlsstr{multiplicity},'(',num2str(numlvls(multiplicity)),')']);
    else
      if(numlvls(multiplicity)>0)
        pos = get(hObject,'CurrentPoint');
        numlvls(multiplicity) = numlvls(multiplicity) - 1;
        set(hLvlsButs{multiplicity},'String',[lvlsstr{multiplicity},'(',num2str(numlvls(multiplicity)),')']);
        for i=1:multiplicity
          line('XData',[0.1 0.4],'YData',[1 1].*(pos(1,2)+(i-1)*(diff(get(hObject,'YLim'))/100)),'Tag','SAFiCF_LevelsLine');
        end
        efit = [efit ones(1,multiplicity).*pos(1,2)];
      else
        set(hLvlsButs{multiplicity},'Visible','off');
      end
    end
    set(hdlLvls,'Visible','off');
  end
end % cb_lvls
%------------------------------------------------------------------------------------------------------------------------------------%
function cb_erange(hObject,eventdata,index)
  limits = get(hLevels,'YLim');
  limits(index) = str2num(get(hObject,'String'));
  set(hLevels,'Ylim',limits);
end % cb_erange
%------------------------------------------------------------------------------------------------------------------------------------%
function cb_fitengy(hObject,eventdata)
  if(J~=0 && length(efit)==(2*J+1))
    if(~isempty(parsfit))
      symmval = get(hSymmPop,'val');
      if(symmval==5 || symmval>37);   % Cubic!
        constraints = zeros(15); constraints(8,4) = 5; constraints(13,9)=-21; 
        parsfit = fitengy(J,parsfit{1},parsfit{2},parsfit{3},efit,[4 9],constraints);
      else
        parsfit = fitengy(J,parsfit{1},parsfit{2},parsfit{3},efit);
      end
    elseif(~isempty(allowed))
      hwb = waitbar(0,'Estimating CF parameters','Name','Estimating CF parameters');
      parsfit = norm2stev(J,{2.*(rand(1,5)-0.5).*allowed{1} 2.*(rand(1,9)-0.5).*allowed{2} 2.*(rand(1,13)-0.5).*allowed{3}});
      constraints = zeros(15); constraints(8,4) = 5; constraints(13,9)=-21; 
      for iter = 1:20
        symmval = get(hSymmPop,'val');
        if(symmval==5 || symmval>37);   % Cubic!
          [parsfit,minsq] = fitengy(J,parsfit{1},parsfit{2},parsfit{3},efit,[4 9],constraints);
        else
          [parsfit,minsq] = fitengy(J,parsfit{1},parsfit{2},parsfit{3},efit);
        end
        if(sqrt(minsq)>10)
          parsfit = norm2stev(J,{2.*(rand(1,5)-0.5).*allowed{1} 2.*(rand(1,9)-0.5).*allowed{2} 2.*(rand(1,13)-0.5).*allowed{3}});
        elseif(sqrt(minsq)<small)
          break;
        end
        waitbar(iter/20,hwb);
      end % iter
    close(hwb);
    end
    enew = eig(cf_hmltn(J,parsfit)); tE = enew; enew = enew-min(enew); degen = [0 0 0 0];
    while(tE)
      % Works out degeneracy of the levels
      idegen = find(abs(tE-tE(1))<(small*10));
      tE(idegen) = [];
      degen(length(idegen)) = degen(length(idegen))+1;
    end
    % Plots fitted levels
    iE = 1; delete(findobj('Tag','SAFiCF_FitsLine'));
    for i=1:length(degen)
      for id = 1:degen(i);
        for idegen=1:i
          line('XData',[0.5 0.8],'YData',[1 1].*(enew(iE)+(idegen-1)*(diff(get(hLevels,'YLim'))/100)),'Tag','SAFiCF_FitsLine','Color','r');
          iE = iE+1;
        end
      end
    end
    set(hLevels,'Xlim',[0 1]);
  end
  DisplayCFPars(hCFParsTxt,parsfit);
end % cb_fitengy
%------------------------------------------------------------------------------------------------------------------------------------%

%------------------------------------------------------------------------------------------------------------------------------------%
function cb_trans(hObject,eventdata)
  if(isempty(parsfit)); cb_fitengy(hEstCFParBt,[]); end
  trans = cflvls(cf_hmltn(J,parsfit),str2num(get(hTcalcEdt,'String')),[0 1])
  % TODO Draw transition matrix elements on!
end % cb_trans
%------------------------------------------------------------------------------------------------------------------------------------%
function cb_lvlbt(hObject,eventdata)
  if(strcmp(get(hdlLvls,'Visible'),'on'))
    multiplicity = str2num(get(hdlLvlsTxt,'Tag'));
    PlotMultiplet(multiplicity);
    set(hdlLvls,'Visible','off');
  end
end % cb_figbt
%------------------------------------------------------------------------------------------------------------------------------------%
function cb_varentry(hObject,eventdata,variable)            % Assigns the value in an edit box to the specified variable
  if(strcmp(variable,'H')) 
    ask = questdlg('Do you want to change the field direction or L and S?','Change Field?','Direction','L and S','Neither','Neither'); 
    if(strcmp(ask,'Neither'))
    elseif(strcmp(ask,'Direction'))
      Jdir = []; 
      getLSJdir;
    elseif(strcmp(ask,'L and S'))
      L=0; S=0; 
      getLSJdir;
    end
  end
  eval([variable '=' get(hObject,'String') ';']);
end % cb_figbt
%------------------------------------------------------------------------------------------------------------------------------------%
function cb_physprop(hObject,eventdata,calcflag)            % Calculates the physical properties and plots it on lower left figure
  if(J~=0 && ~isempty(parsfit))
    delete(findobj('Tag','SAFiCF_physproplines')); legend(hPhysProp,'off'); hold(hPhysProp,'off'); legendstr = {};
    if(calcflag==3)           % Calculates heat capacity
      if(H>0)
        getLSJdir;
        if(length(H)>1)
          for iH = 1:length(H)
            hplt = plot(hPhysProp,T,heatcap(eig(cf_hmltn(J,parsfit) + zeeman_hmltn_dir([L S J],H(iH),Jdir)),T)); 
	    legendstr{iH} = ['Cp(T,H=' num2str(H(iH)) 'T)']; hold(hPhysProp,'on'); set(hplt,'Tag','SAFiCF_physproplines');
          end
	  hold(hPhysProp,'off'); legend(hPhysProp,legendstr);
        else
	  hplt = plot(hPhysProp,T,heatcap(eig(cf_hmltn(J,parsfit) + zeeman_hmltn_dir([L S J],H,Jdir)),T)); 
	  set(hplt,'Tag','SAFiCF_physproplines');
        end
      else
        hplt = plot(hPhysProp,T,heatcap(eig(cf_hmltn(J,parsfit)),T)); set(hplt,'Tag','SAFiCF_physproplines');
      end  % H>0
      xlabel(hPhysProp,'T(K)'); ylabel(hPhysProp,'Cp (J/mol/K)');
    elseif(calcflag<5)
      getLSJdir;
      moment = mag([L S J],H,Jdir,T,cf_hmltn(J,parsfit));
      if(calcflag==4)         % Calculates the magnetisation
        if(length(T)>1)
          for iT = 1:length(T)
            hplt = plot(hPhysProp,H,moment(:,iT)); hold(hPhysProp,'on');
            legendstr{iT} = ['M(H,T=' num2str(T(iT)) 'K)']; set(hplt,'Tag','SAFiCF_physproplines'); 
          end
	  hold(hPhysProp,'off'); legend(hPhysProp,legendstr);
        else
          hplt = plot(hPhysProp,H,moment); set(hplt,'Tag','SAFiCF_physproplines');
        end
        xlabel(hPhysProp,'H(T)'); ylabel(hPhysProp,'M(\mu_B/atom)');
      else
        if(calcflag==1)       % Calculates the susceptibility
          if(length(H)>1); for iH=1:length(H); susc{iH} = moment(iH,:)./H(iH); end; else susc = moment./H; end; suscstr = 'M/H';
        elseif(calcflag==2)   % Calculates the inverse susceptibility
          if(length(H)>1); for iH=1:length(H); susc{iH} = H(iH)./moment(iH,:); end; else susc = H./moment; end; suscstr = 'H/M';
        end
        if(length(H)>1)
          for iH = 1:length(H)
            hplt = plot(hPhysProp,T,susc{iH}); hold(hPhysProp,'on');
	    legendstr{iH} = [suscstr '(T,H=' num2str(H(iH)) 'T)']; set(hplt,'Tag','SAFiCF_physproplines');
          end
	  hold(hPhysProp,'off'); legend(hPhysProp,legendstr);
        else
          hplt = plot(hPhysProp,T,susc); set(hplt,'Tag','SAFiCF_physproplines');
        end
        xlabel(hPhysProp,'T(K)'); ylabel(hPhysProp,suscstr);
      end % calcflag==4
    end
  end
end % cb_physprop
%------------------------------------------------------------------------------------------------------------------------------------%
function cb_close(hObject,eventdata)
  if(ishandle(hdlLvls)); close(hdlLvls); end;
  closereq;
end % cb_close
%------------------------------------------------------------------------------------------------------------------------------------%
function cb_getspec(hObject,eventdata)
  set(hdlLvls,'Name','Obtain spectra from figure.'); 
  hFp = get(hFig,'Position'); set(hdlLvls,'Position',[hFp(1)+hFp(3)/3 hFp(2)+hFp(4)/3 300 150]); 
  set(hdlLvls,'Visible','on'); set(hdlLvlsTxt,'String','Please click on a figure with data of the spectra to be fitted.');
  set(hdlLvlsEdt,'Visible','off'); waitforbuttonpress; hSpec = gcf; set(hdlLvlsEdt,'Visible','on'); set(hdlLvls,'Visible','off');
  [xdat,ydat,edat] = saficf_getdata(hSpec);
  errorbar(hSpectra,xdat,ydat,edat); xlabel(hSpectra,'Energy Transfer (meV)'); ylabel(hSpectra,'Intensity (arb. units)');
  mx = max(ydat(find(xdat>(max(xdat)/5)))); lg = log10(mx); lg = ceil(abs(lg))*sign(lg); set(hSpectra,'YLim',[0 ceil(mx/10^lg)*10^lg]);
end % cb_getspec

%------------------------------------------------------------------------------------------------------------------------------------%
%  Utility functions
%------------------------------------------------------------------------------------------------------------------------------------%
function DisplayCFPars(hObject,B)                           % Displays CF parameter values
  strout = '';
  if(size(B,2)==3 && size(B,1)~=3); B=B'; end;
  for k = 1:3
    for q = 1:(2*(2*k)+1)
      for i = 1:size(B,2)
        if(B{k,i}(q)~=0) 
          strout = [strout,'B',num2str(k*2),num2str(q-1-2*k)];
          if(size(B,2)>1) strout = [strout,'(',num2str(i),')=',num2str(B{k,i}(q)),' ']; 
          else strout = [strout,' = ',num2str(B{k,i}(q)),'    ']; end
        end
      end
      if(size(B,2)>1) strout = [strout,'    ']; end
    end
    strout = [strout,sprintf('\n')];
  end
  set(hObject,'String',strout);
end % DisplayCFPars
%------------------------------------------------------------------------------------------------------------------------------------%
function SetAllowedLevels(hObject,B)                        % Sets visible allowed levels buttons
  if(J~=0)
    tE = eig(cf_hmltn(J,B));
    numlvls = [0 0 0 0];
    efit = [];
    delete(findobj('Tag','SAFiCF_LevelsLine'));
    while(tE)
      idegen = find(abs(tE-tE(1))<small);
      tE(idegen) = [];
      numlvls(length(idegen)) = numlvls(length(idegen))+1;
    end
    for i=1:length(numlvls)
      if(numlvls(i)>0) set(hLvlsButs{i},'String',[lvlsstr{i},'(',num2str(numlvls(i)),')']); set(hLvlsButs{i},'Visible','on'); 
      else set(hLvlsButs{i},'Visible','off'); end
    end
  end
end % SetAllowedLevels
%------------------------------------------------------------------------------------------------------------------------------------%
function PlotMultiplet(multiplicity)                        % Plots a degenerate level in the energy levels figure
  set(0,'CurrentFigure',hFig);
  set(hFig,'CurrentAxes',hLevels);
  if(isempty(findobj('Tag','SAFiCF_LevelsLine')))
    for i=1:multiplicity; line('XData',[0.1 0.4],'YData',[1 1].*(i-1)*(diff(get(hLevels,'YLim'))/100),'Tag','SAFiCF_LevelsLine'); end;
    efit = [zeros(1,multiplicity)];
    numlvls(multiplicity) = numlvls(multiplicity) - 1;
    set(hLvlsButs{multiplicity},'String',[lvlsstr{multiplicity},'(',num2str(numlvls(multiplicity)),')']);
  else
    if(numlvls(multiplicity)>0)
      pos = str2num(get(hdlLvlsEdt,'String'));
      lylim = get(hLevels,'YLim');
      if(pos>lylim(2)); set(hLevels,'YLim',[lylim(1) 10*ceil(pos/10)]); end;
      numlvls(multiplicity) = numlvls(multiplicity) - 1;
      set(hLvlsButs{multiplicity},'String',[lvlsstr{multiplicity},'(',num2str(numlvls(multiplicity)),')']);
      for i=1:multiplicity
        line('XData',[0.1 0.4],'YData',[1 1].*(pos+(i-1)*(diff(get(hLevels,'YLim'))/100)),'Tag','SAFiCF_LevelsLine');
      end
      efit = [efit ones(1,multiplicity).*pos];
    else
      set(hLvlsButs{multiplicity},'Visible','off');
    end
  end
end % PlotMultiplet
%------------------------------------------------------------------------------------------------------------------------------------%
function getLSJdir()                                        % Gets the field direction and quantum numbers
  if(L==0 || S==0)
    LS = inputdlg({'Total orbital quantum number L','Total spin quantum number S'},'Input L and S',1);
    L = str2num(LS{1}); S = str2num(LS{2}); 
  end
  if(isempty(Jdir))
    Jd = inputdlg({'Hx','Hy','Hz'},'Input magnetic field direction',1); Jdir = [str2num(Jd{1}) str2num(Jd{2}) str2num(Jd{3})]; 
  end
end
%------------------------------------------------------------------------------------------------------------------------------------%
end % Figure main function
