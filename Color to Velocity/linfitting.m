% This script is for internal use, helps with calculating and evaluating
% linear models for color_mapping_function_fast. Note that these models
% assume max velocity of 1 for easy conversion later


%% Run this section first

% Load color bars from chosen velocity map file and plot 3D representation

load('F:\Data\220124 single subject\20210930\Velocity Maps\Velocity Map - A0287_old.mat')

color_bar = scanner_parameters.cmap;
%%
v = VideoReader("D:\Data\230530 Rat Vasc Reac\US\2023-05-30 13-48-08.mkv"); % SMI
v = VideoReader("D:\Data\231103 Rat SCAR\US\2023-11-03 14-30-29.mkv"); % CDI
vidFrame = read(v,100);
color_bar = vidFrame(221:379, 1533:1550,:);

%%
% velocity_spectrum = linspace(1,-1,size(color_bar,1));
velocity_spectrum = 1:size(color_bar,1);
% color_bar = rgb2hsv(color_bar);

figure
% x:R y:G z:B
% plot3(color_bar(:,:,1), color_bar(:,:,2), color_bar(:,:,3), "LineStyle", 'none', 'Marker', '.');
% x:G y:B z:v
plot3(color_bar(:,:,1), color_bar(:,:,3), (velocity_spectrum), "LineStyle", 'none', 'Marker', '.');


%% Run after first section

% Linear model 1 (only blue)

indlist1 = 66:80;

c_vect = reshape(color_bar(indlist1,:,:), prod(size(color_bar(indlist1,:,:), [1,2])), []);
v_vect = repmat(velocity_spectrum(indlist1)', size(c_vect, 1)/length(velocity_spectrum(indlist1)),1);

mdl1 = fitlm(double(c_vect), v_vect, 'y~x1+x2+x3-1')
figure, plot(mdl1)

% defines next intercept and threshold value
thresh = squeeze(mean(color_bar(66,:,:)))';
next_intercept = predict(mdl1, thresh);

% max(c_vect(:,3))

%% Run after previous section

% Linear model 2 (red and green)
indlist2 = 1:66;

c_vect = reshape(color_bar(indlist2,:,:), prod(size(color_bar(indlist2,:,:), [1,2])), []);
v_vect = repmat(velocity_spectrum(indlist2)', size(c_vect, 1)/length(velocity_spectrum(indlist2)),1);

% uses intercept and threshold from previous
mdl2 = fitlm(double(c_vect)-thresh, v_vect-next_intercept, 'y~x1+x2+x3-1')
figure, plot(mdl2)

% % defines next intercept and threshold value
% thresh2 = squeeze(mean(color_bar(66,:,:)))';
% next_intercept2 = predict(mdl1, thresh2);
% 
% %% Run after previous section
% 
% % Linear model 3 (red and green)
% indlist = 1:53;
% 
% c_vect = reshape(color_bar(indlist,:,:), prod(size(color_bar(indlist,:,:), [1,2])), []);
% v_vect = repmat(velocity_spectrum(indlist)', size(c_vect, 1)/length(velocity_spectrum(indlist)),1);
% 
% % uses intercept and threshold from previous
% mdl3 = fitlm(double(c_vect)-thresh2, v_vect-next_intercept2, 'y~x1+x2+x3-1')
% figure, plot(mdl3)


%% Model for CHI

color_bar = load("CHI Color Map.mat").chi_map;
velocity_spectrum = linspace(1,-1,size(color_bar,1));

c_vect = reshape(color_bar, prod(size(color_bar, [1,2])), []);
v_vect = repmat(velocity_spectrum', size(c_vect, 1)/length(velocity_spectrum),1);

mdlCHI = fitlm(double(c_vect), v_vect, 'y~x1+x2+x3')
figure
subplot(2,1,1), plot(mdlCHI)

c_span = squeeze(mean(color_bar, 2));
v_span = predict(mdlCHI, c_span);
subplot(2,1,2), plot3(c_span(:,1), c_span(:,2), v_span, "LineStyle", 'none', 'Marker', '.', 'MarkerSize',20)
hold on
plot3(color_bar(:,:,1), color_bar(:,:,2), (velocity_spectrum), "LineStyle", 'none', 'Marker', '.', 'MarkerSize',0.2);


%% Run after establishing your models

indlist = 1:80;

c_vect = reshape(color_bar(indlist,:,:), prod(size(color_bar(indlist,:,:), [1,2])), []);
v_vect = repmat(velocity_spectrum(indlist)', size(c_vect, 1)/length(velocity_spectrum(indlist)),1);


% Plot results of combined linear models

% choose color  to plot here
color_choice = [1,3];

c_span = squeeze(mean(color_bar, 2));
% c_span = [c_span; c_span; c_span]';
v_span_1 = predict(mdl1, c_span); % make sure c_span index matches model inputs
v_span_2 = predict(mdl2, c_span-thresh);
% v_span_3 = predict(mdl3, c_span);

coeffs = mdl2.Coefficients{:,1};

% for ii = color_choice
%     figure;
%     plot(v_span_2+next_intercept, c_span(:,ii), 'LineStyle', 'none', 'Marker','.');
%     hold on
%     plot(v_span_1, c_span(:,ii), 'LineStyle', 'none', 'Marker','.')
%     plot(v_vect, c_vect(:,ii), 'LineStyle', 'none', 'Marker','.')
%     legend
% end

total_v = [v_span_2(indlist2)+next_intercept; v_span_1(indlist1)];
actual_v = [velocity_spectrum(indlist2)'; velocity_spectrum(indlist1)'];

error = sqrt(mean((actual_v-total_v).^2));


f1 = figure;
set(gcf, 'Position', [400 200 200 200]);
% scatter3(d_mdl, v_mdl, predict(mdl1, [d_mdl, v_mdl]), 'red','.');
scatter3(c_span(indlist2,color_choice(1)), c_span(indlist2,color_choice(2)), v_span_2(indlist2)+next_intercept,  'Marker','o');
hold on

scatter3(c_span(indlist1,color_choice(1)), c_span(indlist1,color_choice(2)), v_span_1(indlist1), 'Marker','o')
h = scatter3(c_vect(:,color_choice(1)), c_vect(:,color_choice(2)), v_vect, 'Marker','o', 'MarkerEdgeAlpha', 0.2)
h.SizeData = h.SizeData/9;


text(150, 50, 0.75, sprintf('RMSE = \n%.3f', error), 'HorizontalAlignment','center')
zlabel({'Velocity Index'},'Interpreter','tex');
ylabel({'Blue Value'},'Interpreter','tex');
xlabel({'Red Value'},'Interpreter','tex');

ylh = get(gca,'ylabel');
set(ylh, 'Rotation',-27, 'Position',[-50,0,-0.15])

xlh = get(gca,'xlabel');
xlp = get(xlh, 'Position');
set(xlh, 'Rotation',25, 'Position',[0,-50,-0.15])

pos = get(gca, 'Position');

pos = pos+[-0.01 0.07 0, -0.07];
set(gca, 'CameraPositionMode', 'manual')
set(gca, 'CameraPosition',[-1.031636574981015e+03,-1.299676100945465e+03,4.336001514328491]);
set(gca, 'CameraTargetMode','manual');
set(gca, 'CameraTarget',[105.5000, 127.5000, 0.5013]);
set(gca, 'Position', pos)

% IMPORTANT FOR 3D VECTOR GRAPHICS!!!!
set(f1,'renderer','Painters')
print(f1, '-dsvg', 'C:\Users\Denis\OneDrive for Business\data\New stuff\colormap.svg');
saveas(f1, 'C:\Users\Denis\OneDrive for Business\data\New stuff\colormap.eps', 'epsc');

% figure;
% % plot3(c_span(indlist2,color_choice(1))+thresh(color_choice(1)), c_span(indlist2,color_choice(2))+thresh(color_choice(2)), v_span_2(indlist2)+next_intercept,  'LineStyle', 'none', 'Marker','.');
% % plot3(c_span(:,color_choice(1))+coeffs(color_choice(1))*thresh(color_choice(1)), c_span(:,color_choice(2))+coeffs(color_choice(2))*thresh(color_choice(2)), v_span_2+next_intercept,  'LineStyle', 'none', 'Marker','.');
% scatter3(c_span(indlist2,color_choice(1)), c_span(indlist2,color_choice(2)), v_span_2(indlist2)+next_intercept,  'Marker','.');
% hold on
% 
% scatter3(c_span(indlist1,color_choice(1)), c_span(indlist1,color_choice(2)), v_span_1(indlist1), 'Marker','.')
% h = scatter3(c_vect(:,color_choice(1)), c_vect(:,color_choice(2)), v_vect, 'Marker','o', 'MarkerEdgeAlpha', 0.2)
% h.SizeData = h.SizeData/9;
% legend
% % mean(sum()
%% Run after main, you need an old comparator as well

% calculate error between new run and old run

v_new = velocity(:,:,38);
load('C:\Users\Denis\Documents\Data\220124 single subject\20210930\Velocity Maps\Velocity Map - A0287_old.mat');
v_old = velocity(:,:,38);
sum(sum(abs(abs(v_old)-v_new)))/sum(v_new(:)>0)

%% Plot the color pixel values from chosen png image and compare to colorbar from the velocity map

load('C:\Users\Denis\Documents\Data\220124 single subject\20210930\Velocity Maps\Velocity Map - A0287_old.mat');

image0038 = imread('C:\Users\Denis\Documents\Data\220124 single subject\20210930\Images\A0287\image0038.png');
cfull = double(image0038(min(ROI_idx(:,2)):max(ROI_idx(:,2)), min(ROI_idx(:,1)):max(ROI_idx(:,1)),:));
cfull = rgb2hsv(cfull);
c1 = zeros(size(cfull, [1 2]));
c2 = zeros(size(cfull, [1 2]));
c3 = zeros(size(cfull, [1 2]));

v = velocity(:,:,38); % make sure this matches the image you load
v = double(v(min(ROI_idx(:,2)):max(ROI_idx(:,2)), min(ROI_idx(:,1)):max(ROI_idx(:,1))));

c1(v~=0)=cfull(v~=0);
c2(v~=0)=cfull(find(v~=0)+numel(c1));
c3(v~=0)=cfull(find(v~=0)+2*numel(c1));

figure
plot3(c1(:), c2(:), c3(:), "LineStyle", 'none', 'Marker', '.');

%% Run after section above, plots color against velocity

figure
plot3(c2(:), c3(:), abs(v(:)), "LineStyle", 'none', 'Marker', '.');




%% Run after first section

% Middle linear model, no utility that I've found

indlist = 57:63;

c_vect = reshape(color_bar(indlist,:,:), prod(size(color_bar(indlist,:,:), [1,2])), []);
v_vect = repmat(velocity_spectrum(indlist)', size(c_vect, 1)/length(velocity_spectrum(indlist)),1);

mdl = fitlm(double(c_vect(:,1:3)), v_vect)
figure, plot(mdl)


%% Run after first section

% my failed attempt at a piecewise linear function, the x-threshold does
% not optimize (Jacobian issues)

indlist = 1:80;

color_bar = scanner_parameters.cmap;
velocity_spectrum = linspace(1,-1,size(color_bar,1));


c_vect = double(reshape(color_bar(indlist,:,:), prod(size(color_bar(indlist,:,:), [1,2])), []));
v_vect = repmat(velocity_spectrum(indlist)', size(c_vect, 1)/length(velocity_spectrum(indlist)),1);

% model = @(P,x) P(1) + P(2)*x(1) + P(3)*x(2) + P(4)*x(5) + ...
%             plusfun(P(5)-x(3))*(P(6)*x(1) + P(7)*x(2) + P(8)*x(5))

beta = nlinfit(c_vect, v_vect, @pwise, [20 1 1 1 1 1 1 1]);

% t_span = 1:
%% Run after establishing your piecewise model

% Plot result of failed piecewise

c_span = mean(color_bar, 2);
% c_span = [c_span; c_span; c_span]';
v_span = pwise(beta, c_span);
figure;
plot(v_span, c_span(:,3), 'LineStyle', 'none', 'Marker','.');
hold on
plot(v_vect, c_vect(:,3), 'LineStyle', 'none', 'Marker','.')



