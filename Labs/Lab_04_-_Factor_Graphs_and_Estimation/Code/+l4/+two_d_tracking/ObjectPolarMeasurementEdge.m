classdef ObjectPolarMeasurementEdge < g2o.core.BaseUnaryEdge
    % ObjectPolarMeasurementEdge summary of ObjectPolarMeasurementEdge
    %
    % This class stores an edge which represents the factor for observing
    % the vertex position using a range bearing sensors.
    %
    % The measurement model is
    %
    %    z_(k+1)=h[x_(k+1)]+w_(k+1)
    %
    % The measurements are r_(k+1) and beta_(k+1) and are given as follows.
    % The sensor is at (sx, sy, stheta).
    %
    %    dx = x(k+1) - sx; dy = y(k+1) - sy
    %
    %    r(k+1) = sqrt(dx^2+dy^2)
    %    beta(k+1) = atan2(dy, dx) - stheta
    %
    % The error term
    %    e(x,z) = z(k+1) - h[x(k+1)]
    %
    % However, remember that angle wrapping is required, so you will need
    % to handle this appropriately in compute error.

    properties(Access = protected)
        % The x,y and theta of the sensor
        sensorPose; % (3x1 double vector)
    end
    
    methods(Access = public)
    
        function obj = ObjectPolarMeasurementEdge()
            % ObjectPolarMeasurementEdge for ObjectPolarMeasurementEdge
            %
            % Syntax:
            %   obj = ObjectPolarMeasurementEdge()
            %
            % Description:
            %   Creates an instance of the ObjectPolarMeasurementEdge object.
            %   This is an observation of the particle's position.
            %
            % Outputs:
            %   obj - (handle)
            %       An instance of a ObjectPolarMeasurementEdge
            
            % Call the base constructor, specifying that this has a single
            % vertex and the dimension of the measurement is 2
            obj = obj@g2o.core.BaseUnaryEdge(2);

            % Set to a default value
            obj.sensorPose = zeros(3, 1); 
        end
        
        function setSensorPose(obj, sensorPose)
            % SETSENSORPOSE Set the sensor pose.
            %
            % Syntax:
            %   obj.setSensorPose(sensorPose);
            %
            % Description:
            %   The range bearing sensor, sits in a given position and is
            %   oriented in a particular direction. This function sets the
            %   pose of the sensor.
            %
            % Inputs:
            %   sensorPose - (3x1 double vector)
            %       The field are (x,y,theta) [theta in radians]

            obj.sensorPose = sensorPose;
        end
        
        function computeError(obj)
            % computeError Compute the error for the edge.
            %
            % Syntax:
            %   obj.computeError();
            %
            % Description:
            %   Compute the value of the error, which is the difference
            %   between the predicted and actual range-bearing measurement.
            %

            % Compute the error
            % warning('ObjectPolarMeasurementEdge.computeError: complete implementation')
            dx = obj.edgeVertices{1}.x(1) - obj.sensorPose(1);
            dy = obj.edgeVertices{1}.x(3) - obj.sensorPose(2);
            dXY = [dx; dy];
            r = norm(dXY); % = sqrt(dx^2 + dy^2)
            
            beta = atan2(dy, dx) - obj.sensorPose(3);



            obj.errorZ(1) = obj.z(1) - r;
            obj.errorZ(2) = g2o.stuff.normalize_theta(obj.z(2) - beta);
        end
        
        function linearizeOplus(obj)
            % linearizeOplus Compute the Jacobian of the error in the edge.
            %
            % Syntax:
            %   obj.linearizeOplus();
            %
            % Description:
            %   Compute the Jacobian of the error function with respect to
            %   the vertex.
            %
            

            dx = obj.edgeVertices{1}.x(1) - obj.sensorPose(1);
            dy = obj.edgeVertices{1}.x(3) - obj.sensorPose(2);
            dXY = [dx; dy];

            r = norm(dXY); % = sqrt(dx^2 + dy^2)
            
            % warning('ObjectPolarMeasurementEdge.linearizeOplus: complete implementation')
            j11 = -dx / r; % dr/dx
            j12 = 0; %  dr/d(xdot)
            j13 = -dy / r; % dr/dy
            j14 = 0; % dr/d(ydot)

            j21 = dy / (r^2); % d(beta)/dx
            j22 = 0; % d(beta)/d(xdot)
            j23 = -dx / (r^2); % d(beta)/dy
            j24 = 0; % d(beta)/d(ydot)

            obj.J{1} = [j11 j12 j13 j14; j21 j22 j23 j24];
        end        
    end
end