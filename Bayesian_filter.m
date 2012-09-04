%Example adapted from Michael A. Goodrich and StudentDave
% URL: http://students.cs.byu.edu/~cs470ta/goodrich/fall2005/MATLAB/BayesEstimator.html

s=[3;5]; % Actual state of our system.
N = 100; % No. of iteration.

figure(1);clf;
%figure(2);clf;

figure(1);
h = plot(s(1), s(2), 'ro');  % Plot the state of the system as a red circle.
set(h,'markersize',6,'linewidth',3); % make the circle big.
axis([0,10,0,10]); % Set the scale of the plot so that we can see the origin.
hold off;
hold on
n = 2*randn(2, 100); % Create a 100-sample gaussian noise sequence with a standard deviation of 2.
x = zeros(2, 100);
for i = 1:N
    x(:,i) = s + n(:, i);  % Add the noise to the true state to create 100 observations of the true state.
    plot(x(1,i),x(2,i),'k.');
    pause
end;


% Positions where our system can be in.
Sa = [2:0.05:4];
Sb = [4:0.05:6];
L = length(Sa);

% Pr is prior belief
% Po is posterior belief
Pr = ones(L,L); % Initialize the table to all ones
Po = ones(L,L);
Pr = Pr/sum(sum(Pr)); % Turn the table into a pmf by dividing by the sum.
Po = Po/sum(sum(Po)); % Each entry is now 1/60.
%Pr=0*Pr;Pr(2,2)=1;
figure(1);clf;mesh(Po), axis([0 40 0 40 0 0.015])
pause

[a,b]=find(Po==max(max(Po)));  % Pull out the indices at which Po achieves its max to start.
sest=[Sa(a);Sb(b)];  % The best estimate of the true state to start.
figure(1);
clf
figure(2);
clf
subplot(211); plot(1,sest(1)); hold on;
line([1,N],[s(1),s(1)]); % Draw a line at the location of the x dimension.
subplot(212); plot(1,sest(2)); hold on;
line([1,N],[s(2),s(2)]); % Draw a line at the location of the y dimension.
pause
K=[4,0;0,4]; % covariance matrix for making a 2-D gaussian
for (n=2:length(x));
    Pr=Po; %store the posterior to the prior.
    m=0*Pr;   
    %likelihood
    % look at each location, assume that the given location is
    % is there quail is, and get the likelihood of the data x(:,n) assuming
    % 2-d gaussian noise
    for i=1:length(Pr)
       for j=1:length(Pr)
           me=[Sa(i);Sb(j)];
           m(i,j) = 1/sqrt((2*pi)^2*det(K)) * exp(-(x(:,n)-me)'*inv(K)*(x(:,n)-me)/2); %Compute likelihood           
           m(i,j) = m(i,j) * Pr(i,j); % Combine this likelihood with the prior   
       end;
    end;
    Po=m/sum(sum(m)); %normalize this distribution to make it a proper probability distribution.
    figure(1);mesh(Po), axis([0 40 0 40 0 0.015]) %plot it
    figure(2);
    [a,b]=find(Po==max(max(Po)));  % Get the peak value; it's most likely location of the Quail.
    sest=[Sa(a);Sb(b)];  %A store the coordinates of this location in the bushes
    subplot(211);plot(n,sest(1),'k.');axis([0 N 2 4 ])
    subplot(212); plot(n,sest(2),'k.');axis([0 N 4 6 ])
    pause
end;  
subplot(211); hold off;
subplot(212); hold off;