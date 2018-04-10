classdef Run
    %RUN Data from a particular run
    %   Detailed explanation goes here
    
    properties
        ang
        diode
        m3
        current
        spectrum
    end
    
    methods
        function obj = Run(filename)
            % Run(filename) constuctor, read data file given full path
            
            if nargin == 0
                C{1}=[];
                C{2}=[];
                C{3}=[];
                C{4}=[];
            else
                fid = fopen(filename);
                C = textscan(fid,'%n%n%n%n','headerlines',1);
                fclose(fid);
            end
            obj.ang = C{1};
            obj.diode = C{2};
            obj.m3 = C{3};
            obj.current = C{4};
            obj.spectrum = [obj.ang obj.diode];
        end
    end
    
end

