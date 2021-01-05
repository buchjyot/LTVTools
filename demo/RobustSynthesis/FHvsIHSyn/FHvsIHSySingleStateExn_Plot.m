%% Create Data
data_case = 1;
switch data_case
    case 1
        wcgainFH = [0.0344, 0.0356, 0.0366, 0.0402, 0.0439, 0.0428, 0.0366, 0.0369, 0.056, 0.056];
        Tall = [1, 3, 10, 50, 70, 100, 200, 250, 300, 500];
    case 2
        wcgainFH = [0.0344, 0.0356, 0.0366, 0.0402, 0.0428, 0.056, 0.056];
        Tall = [1, 3, 10, 50, 100, 300, 500];
end
NT = length(Tall);
wcgainIH = 0.0563*ones(NT,1);

%% Plot
figure;clf;
plot(Tall,wcgainFH,'-bo',Tall,wcgainIH,'r','LineWidth',2);
grid on;box on;xlim([Tall(1) Tall(end)]);
xlabel('Horizon (T) (sec.)');
ylabel('Worst-Case Induced L_2 Gain');
legend('Finite Horizon Synthesis','Infinite Horizon Synthesis','Location','southeast')