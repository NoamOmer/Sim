classdef GradPulse
    % GradPulse: Master class of gradient pulses which implements a simple
    % rectangular pulse
    
    properties % everything is public
        StartTime = 0 % in mus
        TotalTime = 1 % in mus
        Amplitude = 0 % in T/m
        Axis = 1% 1, X-axis, 2, Y-axis, 3, Z-axis
        Moment = 0 % in T/m s
        s = 10; % GRT in mus
        GradColor = 'k';
    end
    
    methods
        % simple constuctor
        function obj = CalcAmpl(obj,StartTime,TotalTime,Amplitude,Axis)
            obj.StartTime=StartTime;
            obj.TotalTime=TotalTime;
            obj.Amplitude=Amplitude;
            obj.Axis=Axis;
            obj.Moment = Amplitude * TotalTime/1e6;
        end;
        function obj = CalcMoment(obj,TotalTime,StartTime,Moment,Axis)
            obj.StartTime=StartTime;
            obj.TotalTime=TotalTime;
            obj.Amplitude=Moment/ (TotalTime/1e6);
            obj.Axis=Axis;
            obj.Moment = Moment;
        end;
        % intstep integration step, e.g. 1 ms (this should normally be fine enough)
        % returns a vector containing the gradient shape on timesteps intstep
        function gshape = GradShape(obj,intstep) % in T/m
            gshape = obj.Amplitude * ones(1,ceil(obj.TotalTime/intstep));
        end;
        function display(obj)
            disp(['Moment ' num2str(obj.Moment) '  Ampl ' num2str(obj.Amplitude) '  Start  ' num2str(obj.StartTime) '   Totaltime ' num2str(obj.TotalTime)  ' Axis ' num2str(obj.Axis)]);
        end;
    end
    
end

