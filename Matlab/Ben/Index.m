classdef Index
    %INDEX Index of refraction from a .nk file
    %   constructor with material as string
    
    properties
        wavelength % in nm
        index      % complex number by linear interpolation
    end
    
    methods
        function obj = Index(material, hlines)
            %argument is material name as a string
            filename = [material '.nk'];
            fid = fopen(filename);
            C=textscan(fid, '%n%n%n','headerlines',hlines);
            fclose(fid);
            obj.wavelength=C{1}/10;
            obj.index = C{2}+C{3}*1i;
        end
        
        function n = at(obj,lambda)
            %N returns complex index of refraction at
            %  wavelength lambda in nm.
            n = interp1(obj.wavelength, obj.index, lambda);
        end
    end
    
end

