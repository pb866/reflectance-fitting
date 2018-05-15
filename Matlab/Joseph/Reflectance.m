classdef Reflectance
    %REFLECTANCE Reflectance for a given wavelength
    %   Gain is taken from the log file
    
    properties
        spect % reflectance spectrum
        lambda % wavelength in nm
        fracs  % fraction s polarization at this wavelength
    end
    
    methods
        function obj = Reflectance(log, runs, i, i0, dark)
            % Reflectance constructor, create reflectance data
            %     Args: log -- Log object with log file
            %           runs -- array of Run data
            %           i -- array of i run numbers
            %           i0 -- run number of i0 data
            %           dark -- array of dark currents
            
            if nargin == 0
                i=[];
            end
            ang=[];
            refl=[];
            obj.lambda = log.mono(i(1));
            obj.fracs = polarfind(obj.lambda);
            ir=interp1(runs(i0).ang, runs(i0).diode,obj.lambda);
            r0=(ir-dark(log.gain(i0)))/10^log.gain(i0);
            for j = 1:length(i)
                ang=[ang; runs(i(j)).ang];
                gain = log.gain(i(j));
                rf = (runs(i(j)).diode-dark(gain))/10^log.gain(i(j));
                refl=[refl; rf/r0];
            end
            % Combine angle and reflectance and sort
            obj.spect=sortrows([ang refl]);
        end
        
        function obj = filter(obj, minang, maxang)
            %TRUNCATE truncate the data to be between minimum
            %         and maximum angles.
            fltr = obj.spect(:,1)>=minang & obj.spect(:,1)<=maxang;
            obj.spect = obj.spect(fltr,:);
        end
        
        function obj = plus(obj1, obj2)
            %PLUS overloads + operator to compine two Reflectance objects
            if(abs(obj1.lambda-obj2.lambda)>0.1)
                error('wavelengths don''t agree')
            end
            obj = obj1;
            obj.spect = sortrows([obj1.spect; obj2.spect]);
        end
    end
    
end

