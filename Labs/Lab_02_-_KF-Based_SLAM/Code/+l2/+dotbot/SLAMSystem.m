classdef SLAMSystem < ebe.slam.SLAMSystem

    properties(Access = protected)

        % Although in the lectures we presented three types of states
        % (estimated, predicted, partial) from an implementation point of
        % view it's easier just to store the one set
        x; % The state vector
        P; % The covariance matrix
        u; % The control input

        % Map stores landmark ID with the indices in the state vector
        landmarkIDStateVectorMap;

        % The linear system used to predict the state and the observation
        systemModel;

        % The map
        map

        % Store of the mean and covariance values
        timeStore;
        xStore;
        PStore;

        % For activity 4
        muckUp;
    end

    methods(Access = public)

        % Constructor for the SLAM system
        function obj = SLAMSystem(config)

            % Call base class
            obj@ebe.slam.SLAMSystem(config);

            % Set up the discrete time system for prediction
            obj.systemModel = l2.dotbot.SystemModel(config);

            obj.muckUp = false;

            % Set up the event handlers
            obj.registerEventHandler('init', @obj.handleInitializationEvent); % Initialization
            obj.registerEventHandler('null_obs', @obj.handleNoUpdate); % No observation
            obj.registerEventHandler('gps', @obj.handleGPSObservationEvent); % GPS observation
            obj.registerEventHandler('slam', @obj.handleSLAMObservationEvent); % SLAM observation
            obj.registerEventHandler('bearing', @obj.handleBearingObservationEvent); % Bearing observation
            obj.registerEventHandler('odom', @obj.handleUpdateOdometryEvent); % Odometry update

            % Set the name
            obj.setName('SLAMSystem');
        end

        % Start the SLAM system and set up the state vector
        % and the dictionary which maps landmark ID to the state vector index.
        % Returns true when the system was starts successfully, false otherwise.
        function success = start(obj)
            start@ebe.slam.SLAMSystem(obj); % Call the base class start method
            obj.timeStore = []; % Store the time
            obj.xStore = zeros(l2.dotbot.SystemModel.NP, 0); % Store the state
            obj.PStore = zeros(l2.dotbot.SystemModel.NP, 0); % Store the covariance

            % Set the dictionary which maps landmark ID to coefficient in the state estimate.
            obj.landmarkIDStateVectorMap = configureDictionary("uint32", "double");

            % Get the map data
            if (isfield(obj.config, 'map'))
                obj.map = obj.config.map;
            end

            success = true;
        end

        % Return the current platform estimate and the covariance
        function [x,P] = platformEstimate(obj)
            x = obj.x(1:l2.dotbot.SystemModel.NP); % From 1 to number of platform states (NP initially = 2)
            P = obj.P(1:l2.dotbot.SystemModel.NP, 1:l2.dotbot.SystemModel.NP);  % From 1 to number of platform states (NP initially = 2)
        end
        
        % Return the stored time, state and covariance values
        function [T, X, PX] = platformEstimateHistory(obj)
            T = obj.timeStore; % Time
            X = obj.xStore; % State
            PX = obj.PStore; % Covariance
        end
        
        % Return the current landmark estimates and the covariance values
        function [x, P, landmarkIds] = landmarkEstimates(obj)

            landmarkIds = keys(obj.landmarkIDStateVectorMap); % Get the landmark IDs
            numberOfLandmarks = numel(landmarkIds); % Get the number of landmarks based on the IDs
           
            % Initialize the arrays as NaN values to store the landmark estimates and the covariance values
            % x has elements equal to the number of landmark states (NL = 2) and the number of landmarks
            % P has elements equal to (NL x NL x number of landmarks)
            x = NaN(l2.dotbot.SystemModel.NL, numberOfLandmarks);  
            P = NaN(l2.dotbot.SystemModel.NL, l2.dotbot.SystemModel.NL, numberOfLandmarks); 
            
            % For each landmark
            for l = 1 : numberOfLandmarks
                % Get the landmark ID
                landmarkId = landmarkIds(l);
                % Get the offset and the index of the landmark in the state vector
                offset = lookup(obj.landmarkIDStateVectorMap, landmarkId);
                idx = offset + [1;2];
                % Store the landmark estimate and the covariance values
                x(:, l) = obj.x(idx);
                P(:, :, l) = obj.P(idx, idx);
            end
            
        end

        % Return the stored time, state and covariance values.
        function [T, X, PX] = estimateHistory(obj)
            T = obj.timeStore;
            X = obj.xStore;
            PX = obj.PStore;
        end

        function muckUpCovarianceMatrix(obj, muckUp)
            obj.muckUp = muckUp;
        end
    end

    methods(Access = protected)

        function success = handleNoPrediction(obj)
            % if there is no prediction, keep the state and the covariance matrix as they are
            obj.x = obj.x;
            obj.P = obj.P;
            success = true;
        end

        function success = handleNoUpdate(obj, ~)
            % if there is no update, keep the state and the covariance matrix as they are
            obj.x = obj.x;
            obj.P = obj.P;
            success = true;
        end

        function success = handlePredictForwards(obj, dT)

            NP = l2.dotbot.SystemModel.NP;

            [obj.x(1:NP), FXd, QXd] = obj.systemModel.predictState(obj.x(1:NP), obj.u, dT);

            % This is an inefficient way to do it:
            %
            % FS = eye(numel(obj.x));
            % FS(1:NP, 1:NP) = FXd;
            % QS = eye(numel(obj.x));
            % QS(1:NP,1:NP) = QXd;
            % obj.P = FS * obj.P * FS' + QS
            %
            % The more efficient approach is exploit the structure of FS
            % and QS and expand by hand

            % Multiply the top left block for the platform state
            obj.P(1:NP,1:NP) = FXd * obj.P(1:NP, 1:NP) * FXd' + QXd;

            % Do the platform landmark-prediction blocks
            obj.P(1:NP, NP+1:end) = FXd * obj.P(1:NP, NP+1:end);
            obj.P(NP+1:end, 1:NP) = obj.P(1:NP, NP+1:end)';

            success = true;
        end

        function success = handleInitializationEvent(obj, event)
            obj.x = event.data;
            obj.P = event.covariance;
            obj.initialized = true;
            success = true;
        end

        function success = handleGPSObservationEvent(obj, event)
            [zPred, Hx, R] = obj.systemModel.predictGPSObservation(obj.x(1:l2.dotbot.SystemModel.NP));

            % Expand to the full state
            HS = zeros(2, numel(obj.x));
            HS(:, 1:l2.dotbot.SystemModel.NP) = Hx;

            % Kalman Filter Update
            nu = event.data - zPred;
            C = obj.P * HS';
            S = HS * C + R;
            W = C / S;
            obj.x = obj.x + W * nu;
            obj.P = obj.P - W * S * W';
            success = true;
        end

        function success = handleBearingObservationEvent(obj, event)

            x = obj.x;
            P = obj.P;

            % Update each measurement separately
            for s = 1 : numel(event.info)
                sensor = obj.map.sensors.bearing.sensors(event.info(s));
                [zPred, Hx, R] = obj.systemModel.predictBearingObservation(x(1:l2.dotbot.SystemModel.NP), sensor.position, sensor.orientation);

                % Expand to full state
                HS = zeros(1, numel(obj.x));
                HS(:, 1:l2.dotbot.SystemModel.NP) = Hx;
                
                % Kalman Filter Update
                nu = event.data(s) - zPred;
                nu = atan2(sin(nu), cos(nu));
                C = P * HS';
                S = HS * C + R;
                W = C / S;
                x = x + W * nu;
                P = P - W * S * W';
            end
            obj.x = x;
            obj.P = P;
            success = true;
        end

        % Handle a set of measurements of landmarks
        function handleSLAMObservationEvent(obj, event)
            % assert: check if the condition is true. If not, display an error message.
            assert(obj.stepNumber == event.eventGeneratorStepNumber)

            % Store useful values
            NL = l2.dotbot.SystemModel.NL; % Number of landmark states
            NP = l2.dotbot.SystemModel.NP; % Number of platform states

            % Get the list of landmarks we know about
            knownLandmarkIDs = obj.landmarkIDStateVectorMap.keys();

            % Find the set of known landmarks in the observation set.
            %
            % Parameters:
            % event.info is the set of landmark IDs observed.
            % knownLandmarkIDs is the set of landmark IDs we know about.
            %
            % Their intersection is set of landmark IDs we know about which we have re-observed.
            [existingLandmarks, idx] = intersect(event.info, knownLandmarkIDs);

            % FIRST: Update the known landmarks
            % numel(existingLandmarks) is the number of known landmarks
            for o = 1 : numel(existingLandmarks)

                % The following two lines look up for each known landmark and
                % what the index of that landmark state is in the state vector
                % obj.landmarkIDStateVectorMap is a dictionary that maps landmark IDs to the index in the state vector
                % existingLandmarks(o) is the landmark ID
                offset = lookup(obj.landmarkIDStateVectorMap, existingLandmarks(o));
                % Add the number of landmark states (NL=2) to the offset to get the index of the landmark in the state vector
                landmarkIdx = offset + (1:NL); 

                % Call the observation model; this returns the predicted observation (zPred)
                % and the observation Jacobian (Hx) for the platform and the landmark (Hm)
                [zPred, Hx, Hm, ~] = ...
                    obj.systemModel.predictSLAMObservation(obj.x(1:NP), ...
                    obj.x(landmarkIdx));

                % ACTIVITY 3:
                %
                % Create the SLAM system observation matrix HS which
                % has to be of the right size and composition to update the
                % landmark and platform.
                HS = zeros(2, numel(obj.x));
                HS(:, 1:2) = -eye(2);
                HS(:, landmarkIdx) = eye(2);
                % Use HS to implement the Kalman filter update for the SLAM
                % system
                % parameters
                Pcov = obj.P;
                Rcov = event.covariance();

                v = zPred - HS * obj.x;
                C = Pcov * HS';
                S = HS * C + Rcov;
                W = C / S;
                obj.x = obj.x + W * v;
                obj.P = Pcov - W * S * W';

            end

            % SECOND: Add the remaining (new) landmarks
            % The difference between the set of observed landmarks and
            % the set of known landmarks is the set of new landmarks
            [newLandmarks, idx] = setdiff(event.info, existingLandmarks);

            % numel(newLandmarks) is the number of new landmarks
            for o = 1 : numel(newLandmarks)

                %% Book keeping in the SLAM system: tie landmark IDs with their position in the state vector
                % Get the length of the state vector. Initially this is the number of platform states (NP=2)
                % but it will grow as we add landmarks
                offset = length(obj.x); 
                % Add the number of landmark states (NL=2) to the offset to get the index of the landmark in the state vector
                landmarkIdx = offset + (1:l2.dotbot.SystemModel.NL);
                % Insert the landmark ID (newLandmarks(o)) and the index of that landmark (offset) 
                % into the dictionary that maps landmark IDs to the index in the state vector
                obj.landmarkIDStateVectorMap = insert(obj.landmarkIDStateVectorMap, newLandmarks(o), offset);
                % obj.x = [obj.x; event.data(:, idx(o))]; % Add the new landmark to the state vector

                % ACTIVITY 2: Insert here the correct code to compute the J
                % and K matrices needed to augment the state with the
                % landmark. 

                l = length(obj.x); % is the length of x before adding this landmark

                % J should look like this: 
                j1 = eye(l);
                j2 = [1, zeros(1, l-1)];
                j3 = [0, 1, zeros(1, l-2)];
                J = [j1; j2; j3];
                % J = [1 0 0 0...0
                %      0 1 0 0...0
                %      0 0 1 0...0
                %      0 .. .. ..0
                %      1 0 0 0...0
                %      0 1 0 0...0]
                
                % K should look like this:
                k1 = [zeros(l,1); 1; 0];
                k2 = [zeros(l,1); 0; 1];
                K = [k1, k2];
                % K = [0 0]
                %     [0 0]
                %     [0 0]
                %     [...]
                %     [1 0]
                %     [0 1]
                obj.x = J * obj.x + K * event.data(:, idx(o));
                obj.P = J * obj.P * J' + K * event.covariance() * K';

                %% Ideal Solution
                % J = zeros(offset + NL, offset);
                % J(1:offset, 1:offset) = eye(offset);
                % J(landmarkIdx, 1:NL) = eye(NL);                
                % K = zeros(offset + NL, NL);
                % K(landmarkIdx, 1:NL) = eye(NL);
                % obj.x = J * obj.x + K * event.data(:, idx(o));
                % obj.P = J * obj.P * J' + K * event.covariance() * K';
            end

            % ACTIVITY 4
            if (obj.muckUp == true)
                for k = 1 : 2 : length(obj.P)
                    obj.P(k:k+1, 1:k-1) = 0;
                    obj.P(k:k+1, k+2:end) = 0;
                end
            end
        end

        function success = handleUpdateOdometryEvent(obj, event)
                assert(obj.stepNumber == event.eventGeneratorStepNumber);
                obj.u = event.data;
                success = true;
            end

        function storeStepResults(obj)
            % Store the estimate for the future
            obj.timeStore(:, obj.stepNumber + 1) = obj.currentTime;
            obj.xStore(:, obj.stepNumber + 1) = obj.x(1:l2.dotbot.SystemModel.NP);
            obj.PStore(:, obj.stepNumber + 1) = diag(obj.P(1:l2.dotbot.SystemModel.NP, ...
                1:l2.dotbot.SystemModel.NP));

        end
    end
end