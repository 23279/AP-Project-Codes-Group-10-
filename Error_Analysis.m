%% 
clear
clc
close all

%% -----------------------------
% STEP 1. Fixed Parameters
%% -----------------------------
alpha = 1;
k1 = 10;
k2 = 50;

%% -----------------------------
% STEP 2. Initial Conditions
%% -----------------------------
A0 = 0;
V0 = 0;
initialConditions = [A0 V0];

%% -----------------------------
% STEP 3. Simulation Settings
%% -----------------------------
startTime = 0;
endTime = 400;
perturbation = 0.05;

%% -----------------------------
% STEP 4. Subject Data
%% -----------------------------
% Columns:
% ID r1 r2 r3 r4 Dose1 Dose2

subjects = [
%ID   r1       r2      r3       r4       vDay1  vDay2
    1     0.019    0.26    0.0140   0.038    0      28;   % Moderna
    3     0.0062   0.40    0.0072   0.027    0      21;   % Pfizer
    4     0.011    0.35    0.0079   0.026    0      21;   % Pfizer
    5     0.045    0.17    0.0140   0.025    0      28;   % Moderna
    19    0.110    0.04    0.0064   0.028    0      21;   % Pfizer
    20    0.026    0.23    0.0140   0.030    0      28;   % Moderna
    21    0.015    0.28    0.0150   0.032    0      28;   % Moderna
    23    0.063    0.12    0.0140   0.023    0      28;   % Moderna
    25    0.030    0.21    0.0140   0.031    0      28;   % Moderna
    28    0.022    0.25    0.0150   0.033    0      28;   % Moderna
    30    0.066    0.12    0.0140   0.022    0      28;   % Moderna
    40    0.026    0.22    0.0150   0.029    0      28;   % Moderna
    41    0.013    0.26    0.0063   0.028    0      21;   % Pfizer
    42    0.0089   0.41    0.0160   0.043    0      28;   % Moderna
    44    0.021    0.18    0.0061   0.028    0      21;   % Pfizer
    64    0.012    0.36    0.0140   0.047    0      28;   % Moderna
    65    0.0048   0.46    0.0130   0.073    0      28;   % Moderna
    76    0.025    0.21    0.0061   0.029    0      21;   % Pfizer
    79    0.035    0.17    0.0140   0.026    0      28;   % Moderna
    94    0.028    0.21    0.0130   0.034    0      28;   % Moderna
    100   0.013    0.29    0.0160   0.033    0      28;   % Moderna
    101   0.0079   0.32    0.0120   0.074    0      28;   % Moderna
    139   0.024    0.22    0.0150   0.028    0      28;   % Moderna
    146   0.021    0.26    0.0069   0.027    0      21;   % Pfizer
    147   0.047    0.15    0.0140   0.026    0      28;   % Moderna
    173   0.039    0.17    0.0140   0.029    0      28;   % Moderna
    175   0.0057   0.43    0.0120   0.073    0      28;   % Moderna
    179   0.011    0.32    0.0110   0.067    0      28;   % Moderna
    200   0.011    0.34    0.0140   0.053    0      28;   % Moderna
    215   0.010    0.21    0.0053   0.029    0      21;   % Pfizer
    216   0.016    0.19    0.0053   0.030    0      28;   % Moderna
    226   0.0061   0.49    0.0069   0.027    0      21;   % Pfizer
% ...
% <<< Paste the remaining rows from your original code here >>>
];

numberOfSubjects = size(subjects,1);

parameterNames = {'r1','r2','r3','r4','alpha','k1','k2'};

increase14 = zeros(numberOfSubjects,7);
decrease14 = zeros(numberOfSubjects,7);
increase400 = zeros(numberOfSubjects,7);
decrease400 = zeros(numberOfSubjects,7);

options = odeset('RelTol',1e-6,'AbsTol',1e-8);

%% -----------------------------
% STEP 5. Loop Through Subjects
%% -----------------------------
for subject = 1:numberOfSubjects

    r1 = subjects(subject,2);
    r2 = subjects(subject,3);
    r3 = subjects(subject,4);
    r4 = subjects(subject,5);

    vaccineDays = [subjects(subject,6) subjects(subject,7)];

    parameters = [r1 r2 r3 r4 alpha k1 k2];

    % ---------- Baseline ----------
    rhs = @(t,y)[ ...
        parameters(1)*y(2) + parameters(2)*y(1)*y(2) + y(1)*(parameters(3)-parameters(4)*y(1));
        parameters(5)*(double(t>=vaccineDays(1)&t<vaccineDays(1)+1)+double(t>=vaccineDays(2)&t<vaccineDays(2)+1))...
        - parameters(6)*y(2)/(parameters(7)+y(2)+eps)];

    [t,y] = ode45(rhs,[startTime endTime],initialConditions,options);

    day14 = vaccineDays(2)+14;

    baseline14 = interp1(t,y(:,1),day14);
    baseline400 = interp1(t,y(:,1),400);

    % ---------- Change one parameter ----------
    for parameter = 1:7

        increased = parameters;
        increased(parameter)=increased(parameter)*1.05;

        rhs = @(t,y)[ ...
            increased(1)*y(2)+increased(2)*y(1)*y(2)+y(1)*(increased(3)-increased(4)*y(1));
            increased(5)*(double(t>=vaccineDays(1)&t<vaccineDays(1)+1)+double(t>=vaccineDays(2)&t<vaccineDays(2)+1))...
            - increased(6)*y(2)/(increased(7)+y(2)+eps)];

        [t,y]=ode45(rhs,[startTime endTime],initialConditions,options);

        A14=interp1(t,y(:,1),day14);
        A400=interp1(t,y(:,1),400);

        increase14(subject,parameter)=100*(A14-baseline14)/max(baseline14,eps);
        increase400(subject,parameter)=100*(A400-baseline400)/max(baseline400,eps);

        decreased = parameters;
        decreased(parameter)=decreased(parameter)*0.95;

        rhs = @(t,y)[ ...
            decreased(1)*y(2)+decreased(2)*y(1)*y(2)+y(1)*(decreased(3)-decreased(4)*y(1));
            decreased(5)*(double(t>=vaccineDays(1)&t<vaccineDays(1)+1)+double(t>=vaccineDays(2)&t<vaccineDays(2)+1))...
            - decreased(6)*y(2)/(decreased(7)+y(2)+eps)];

        [t,y]=ode45(rhs,[startTime endTime],initialConditions,options);

        A14=interp1(t,y(:,1),day14);
        A400=interp1(t,y(:,1),400);

        decrease14(subject,parameter)=100*(A14-baseline14)/max(baseline14,eps);
        decrease400(subject,parameter)=100*(A400-baseline400)/max(baseline400,eps);

    end

end

%% STEP 6. Calculate Mean and SEM
meanIncrease14 = mean(increase14);
meanDecrease14 = mean(decrease14);
meanIncrease400 = mean(increase400);
meanDecrease400 = mean(decrease400);

semIncrease14 = std(increase14)/sqrt(numberOfSubjects);
semDecrease14 = std(decrease14)/sqrt(numberOfSubjects);
semIncrease400 = std(increase400)/sqrt(numberOfSubjects);
semDecrease400 = std(decrease400)/sqrt(numberOfSubjects);

%% STEP 7. Plot Results
figure

subplot(1,2,1)
bar([meanIncrease14' meanDecrease14'])
hold on
errorbar((1:7)-0.15,meanIncrease14,semIncrease14,'.k')
errorbar((1:7)+0.15,meanDecrease14,semDecrease14,'.k')
xticklabels(parameterNames)
title('14 Days After Final Dose')
ylabel('% Change')

subplot(1,2,2)
bar([meanIncrease400' meanDecrease400'])
hold on
errorbar((1:7)-0.15,meanIncrease400,semIncrease400,'.k')
errorbar((1:7)+0.15,meanDecrease400,semDecrease400,'.k')
xticklabels(parameterNames)
title('Day 400')
ylabel('% Change')