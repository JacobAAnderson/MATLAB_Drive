% Demo NOAA Tides
% Jacob Anderson
% 8/15/2019

% Get Tide Data from noaa at a specified staion and time duration

% station     --> NOAA station ID [ doulbe ], find at  https://tidesandcurrents.noaa.gov/products.html
% start       --> Beging of time duration [ datetime ]
% stop        --> End of time duration [ datetime ]
% plotResults --> Optianal input to indicat that the interpolated results should be plotted [ boolean ]

close all
clear all
clc

stationID = 9410079;                        % NOAA station ID
start = datetime('14-July-2016 12:45:07');  % Beging of time duration
stop  = datetime('15-july-2016 15:30:56');  % End of time duration --> 'hour, minuts, seconds' will be discarded 
                                            %                           but are acceptible as part of the input

tide = NOAA_tides( stationID, start, stop, true);
