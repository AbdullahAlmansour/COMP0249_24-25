% Dynamic object example. This shows how to build a graph to estimate an
% object which has both position and velocity

import ebe.graphics.FigureManager;
import g2o.core.*;
import l4.two_d_tracking.*;

% Some parameters
numberOfTimeSteps = 1000; % Number of time steps
dT = 1; % Time step (s)
sigmaR = 100; % Measurement noise
sigmaQ = 0.01; % Process noise

% Work out the state transition equations
F0=[1 dT; 0 1]; % Top left block of state transition matrix
Q0=[dT^3/3 dT^2/2;dT^2/2 dT] * sigmaQ; % top left and bottom right blocks of process noise matrix

F = [F0 zeros(2); zeros(2) F0]; % Full state transition matrix
Q = [Q0 zeros(2); zeros(2) Q0]; % Full process noise matrix
R = eye(2) * sigmaR; % Measurement noise matrix

% Measurement matrix to extract x and y
H = [1 0 0 0; % to extract x
    0 0 1 0]; % to extract y

% Work out the information matrices
omegaR = inv(R); % Information matrix for the measurement
omegaQ = inv(Q); % Information matrix for the process noise

% Ground truth array
trueX = zeros(4, numberOfTimeSteps); % Ground truth state (x, xdot, y, ydot) times numberOfTimeSteps
z = zeros(2, numberOfTimeSteps); % Measurement array (x, y) times numberOfTimeSteps

% First timestep
trueX(2, 1) = 0.1; % Initial x velocity
trueX(4, 1) = -0.1; % Initial y velocity
% and the velocities are zeros (initially)
z(:, 1) = H * trueX(:, 1) + sqrtm(sigmaR) * randn(2, 1);

% Now predict the subsequent steps
% For each time step, we predict the next state and take a measurement
for k = 2 : numberOfTimeSteps
    % x_k = F * x_{k-1} + v
    v = sqrtm(Q) * randn(4, 1);
    trueX(:, k) = F * trueX(:, k - 1) + v;
    w = sqrtm(sigmaR) * randn(2, 1);
    z(:, k) = H * trueX(:, k) + w; %z(1) = x + w(1) and z(2) = y + w(2)
end

% Create the graph
graph = SparseOptimizer(); % SparseOptimiser: provides methods to for underlying graph structure and optimisation
algorithm = GaussNewtonOptimizationAlgorithm(); % GaussNewtonOptimisationAlgorithm: provides the optimisation algorithm to be used
graph.setAlgorithm(algorithm); % Set the optimisation algorithm to be used

% This array contains the set of vertices for the target state over time
v = cell(numberOfTimeSteps, 1); % Cell array of size numberOfTimeSteps x 1

% Now create the vertices and edges

for n = 1 : numberOfTimeSteps
    
    % Create the first object state vertex
    v{n} = ObjectStateVertex(); % ObjectStateVertex: a vertex that contains the state we want to estimate

    % Set the initial estimate.
    v{n}.setEstimate(zeros(4, 1)); % setEstimate(initialEstimate) where initialEstimate is a 4x1 vector of zeros

    % Another way to do it would be to initialize from pairs of
    % measurements. For example:    
    %if (n == 1)
    %    dz = zeros(2, 1);
    %else
    %    dz = z(:, n) - z(:, n-1);
    %end
    %dzdt = dz / dT;
    %v{n}.setEstimate([z(1, n); dzdt(1); z(2, n); dzdt(2)]);
    
    % Added the vertex to the graph.
    graph.addVertex(v{n}); % addVertex(vertex) adds the vertex to the graph
    
    % If this isn't the first vertex, add the process model edge
    if (n > 1)
        % Create the process model edge
        processModelEdge = ObjectProcessModelEdge();
        % Connect the edge to the vertices v{n-1} and v{n}
        processModelEdge.setVertex(1, v{n-1});
        processModelEdge.setVertex(2, v{n});
        % Set the measurement value, the state transition matrix and the
        processModelEdge.setMeasurement([0;0;0;0]);
        % Set the information matrix
        processModelEdge.setF(F);
        % Set the information matrix
        processModelEdge.setInformation(omegaQ);
        % Add the edge to the graph
        graph.addEdge(processModelEdge);
    end
    

    % To simulate the GPS dropout, the code down here should be disabled
    % during the specified timesteps.

    % Create the measurement edge
    e = ObjectGPSMeasurementEdge();
    
    % Link it so that it connects to the vertex we want to estimate
    e.setVertex(1, v{n}); 
    
    % Set the measurement value and the measurement covariance
    e.setMeasurement(z(:,n));
    e.setInformation(omegaR);
    
    % Add the edge to the graph; the graph now knows we have these edges
    % which need to be added
    graph.addEdge(e);
end

% Graph construction complete

% Initialise the optimization. This is done here because it's a bit
% expensive and if we cache it, we can call optimize multiple times later.
graph.initializeOptimization();

% Create some output as we go
x = zeros(4, numberOfTimeSteps);

% First copy the state values - these are the prior we set the graph to
for n = 1 : numberOfTimeSteps
    x(:, n) = v{n}.estimate();
end

% Plot the prior and the truth
ebe.graphics.FigureManager.getFigure('Lab 04 Activity 02 Scenario');
clf

% Plot. Note that we capture the line handle. This is overkill in this
% case, but it's a useful habit to get into for labelling graphs.
gH(1)=plot(x(1, :), x(3,:));
hold on
gH(2)=plot(trueX(1, :), trueX(3, :));

% Optimize the graph
tic
graph.optimize(5000);
toc

% Extract the optimized state estimate and plot
for n = 1 : numberOfTimeSteps
    x(:, n) = v{n}.estimate();
end
gH(3)=plot(x(1, :), x(3, :), 'LineWidth', 2);

gH(4)=plot(z(1,:),z(2,:));

% Generate the legend
legend(gH, {'Prior', 'Truth', 'Optimized', 'Observation'});
title([num2str(graph.chi2())])
drawnow

