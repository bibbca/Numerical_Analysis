%% Dialog Box
clear all

prompt = {'Anonymous Function = ',...
    'Left Endpoint',...
    'Right Endpoint'...
    'Tolerance',...
    'Would you like a plot? [Y/N]'};            % Allows us to read words.
dlg_title = ('Caleb Bibbs Input');              % Sets Dialog Box Title
size_wind = [1 50; 1 50; 1 50; 1 50; 1 50];     % Window size
defaultans = {'@(x) exp(3.*x).*sin(2.*x)',...
    '1','3','10^-5','Y'};                       % Preprogramed Defaults for Convenience.
options.Resize = 'on';

promptInputs = inputdlg(prompt,dlg_title,size_wind, defaultans,options); % Creates our Dialog Box.
Function = str2func(promptInputs{1});           % Collects the Anonymous Function
LeftEnd = str2double(char(promptInputs(2)));    % Collects the Left Endpoint
RightEnd = str2double(char(promptInputs(3)));   % Collects the Right Endpoint
Tolerance = str2num(promptInputs{4});           %#ok<ST2NM> Collects the Tolerance
plotChoice = char(promptInputs(5));             % Collects the plot choice

%% Endpoint Check
if LeftEnd > RightEnd
    disp('Check your bounds next time, but I guess I can fix it for you this once.');
    Switch = LeftEnd;
    LeftEnd = RightEnd;
    RightEnd = Switch;
end

%% Modular Commands
[h_last, Approx] = romberg (LeftEnd, RightEnd,  Tolerance, Function);
[I,h_final,x]= adaptive_quad (LeftEnd, RightEnd,  Tolerance, Function, plotChoice);

%% Outputs:
disp(['For Romberg our Approximation is: ',num2str(Approx),' with ',num2str(h_last),' as our smallest stepsize.']);
disp(['For Adaptive Curvature our Approximation is: ',num2str(I),' with ',num2str(h_final),' as our smallest stepsize.']);