classdef LogFile
    %LOGFILE Contents of ALS log file
    %   There are methods for extracting specific information for each run
    
    properties
        C % cell array with data from log file
        dir % directory from initialization
    end
    
    methods
        function obj = LogFile(directory, prefix)
            % LogFile(directory, prefix) constuctor
            %   Default directory ../Feb2018/
            %   Default prefix Feb2018
            %   Full file name will be [directory prefix '.log']
            
            %% Pre Initialization %%
            % Any code not using output argument (obj)
            if nargin == 0
                directory = '../Feb2018/';
                prefix = 'Feb2018';
            elseif nargin == 1
                prefix = 'Feb2018';
            end
            
            %% Post Initialization %%
            filename = [directory prefix '.log'];
            fid = fopen(filename);
            format = ['%*s%*s%s%s%s%s%*s%*s%n%*s%*s%*s%*s%*s%*s%*s%s' ...
                '%*s%d%*s%*s%s%s%*s%*s%*s%n%*s%n%*s%*s%*s%*s%*s%n%*s' ...
                '%*n%*s%n%*s%n%*s%n%*s'];
            obj.C = textscan(fid, format, 'Delimiter', char(009), ...
                'Headerlines', 2);
            obj.dir = directory;
            fclose(fid);
        end
        
        function name = dataFile(obj, run)
            % DATAFILE(run) filename for the given run
            name = char(obj.C{1}(run));
        end
        
        function name = fullDataFile(obj, run)
            % FULLDATAFILE(run) file name for the given run with path
            name = [obj.dir obj.dataFile(run)];
        end
        
        function c = comment(obj, run)
            %COMMENT(run) comment entered for run
            c = char(obj.C{2}(run));
        end
        
        function d = date(obj, run)
            %DATE(run) date on which data was taken (string)
            d = char(obj.C{3}(run));
        end
        
        function t = time(obj, run)
            %TIME(run) time for run : string
            t = char(obj.C{4}(run));
        end
        
        function g = gain(obj, run)
            %GAIN(run) log of gain for run;
            g = obj.C{5}(run);
        end
        
        function s = scanType(obj, run)
            %SCANTYPE(run) type of scan; t2d, etc.
            s = char(obj.C{6}(run));
        end
        
        function g = grating(obj, run)
            %GRATING(run) selected grating: lines per mm
            g = obj.C{7}(run);
        end
        
        function f = filter1(obj, run)
            %FILTER1(run) first filter
            f = char(obj.C{8}(run));
        end
        
        function f = filter2(obj, run)
            %FILTER2(run) second filter
            f = char(obj.C{9}(run));
        end
        
        function m = mono(obj, run)
            %MONO(run) wavelength in nm : float
            m = obj.C{10}(run);
        end
        
        function s = sample(obj, run)
            %SAMPLE(run) sample angle in deg : float
            s = obj.C{11}(run);
        end
        
        function d = detector(obj, run)
            %DETECTOR(run) detector angle in deg : float
            d = obj.C{12}(run);
        end
        
        function x = samplex(obj, run)
            %SAMPLEX(run) x position of sample in mm : float
            x = obj.C{13}(run);
        end
        
        function y = sampley(obj, run)
            %SAMPLEY(run) y position of sample in mm : float
            y = obj.C{14}(run);
        end
        
        function z = samplez(obj, run)
            %SAMPLEZ(run) z position of sample in mm : float
            z = obj.C{15}(run);
        end
        
    end
    
end

