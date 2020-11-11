
close all
clear all
clc

ui = SavedUserInputs(mfilename);                                            % Instantiate the Saved User Inputs Class
ui = ui.NewData(false);                                                     % Indicate whether new data should be selescted by the user

saveFile = ui.SaveFile('DemoMovie','mp4');                                  % GUI to get file path to save variable / object
if isempty(saveFile), return, end                                           % End script if there isn't an input for it to use


Z = peaks;                                                                  % Create object to be recorded
surf(Z)
axis tight manual
ax = gca;
ax.NextPlot = 'replaceChildren';


loops = 40;

movie = Fig2Movie('Water Feature Sim', loops);                              % Create Movie recorder

for ii = 1: loops
    
    X = sin(ii*pi/10)*Z;                                                     % Update Figure
    surf(X,Z)
    drawnow
    
    
    movie = movie.Add_Frame(gcf);                                           % Record Figure
    
end


FrameRate = 5;                                                              % Frames per second
movie.MakeVideo(saveFile, FrameRate)                                        % Make Movie


disp("Done!!")

