%% Tutorial 4: Basic plotting with Matlab
% Here is an example of some basic plotting tools. You can checkout more in the documentation,
% e.g. type 'help plot' or online
%%

%% Basic plot
x=0:0.1:10; % vector of x values between 0 and 10 with spacing 0.1
figure;
set(gca,'FontSize',15); % make the default font a bit bigger
plot(x,sin(x),'k-'); % plot sin(x)
% put labels on axis:
xlabel('x');  
ylabel('y');

%% add another curve to above plot
hold on; % use this otherwise previous curve gets overwritten
plot(x,cos(x),'ro--');
% to specify x and y limits use:
xlim([0 5]);  % x-limits
ylim([-1 1]); % y-limits
% add title and legend
title('Exampleplot');
legend('sin(x)','cos(x)');


%% Using subplots
figure;
subplot(1,2,1);  % first subplot
set(gca,'FontSize',15); 
semilogy(x,exp(x),'m^-.'); % semilogarithmic plot
xlabel('x');  
ylabel('y');

subplot(1,2,2);  % first subplot
set(gca,'FontSize',15); 
loglog(x,x.^3,'gs:'); % loglog plot
xlabel('x');  
ylabel('y');

