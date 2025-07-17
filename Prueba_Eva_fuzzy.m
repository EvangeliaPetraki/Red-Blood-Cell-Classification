%% Cargamos la imagen (variable I)

clear all
close all

% ----- 
[I,map] = imread('capture_1548244948896.jpg');
[I,map] = imread('capture_1548244759635.jpg');
[I,map] = imread('capture_1548244406801.jpg');
[I,map] = imread('capture_1548244339396.jpg');
[I,map] = imread('capture_1548244150176.jpg');

%validation
[I,map] = imread('capture_1548243990843.jpg');
[I,map] = imread('capture_1548244554243.jpg');


% figure; 
% imshow(I);
% title(['Original Image']);

I_green = im2double(I(:,:,2));
[counts, bins] = imhist(I_green);

% Smooth the histogram to reduce noise
smooth_counts = smooth(counts, 10);

% Find all peaks
[peaks, locs] = findpeaks(smooth_counts);

% Use the first peak (left-most one)
[~, idx_main_peak] = min(locs);  % since it's the left-most
main_peak_height = peaks(idx_main_peak);
main_peak_loc = locs(idx_main_peak);

% Compute Half-Max
half_max = main_peak_height / 2;

% Search to the right of the main peak for where we fall below half max
right_part = smooth_counts(main_peak_loc:end);
right_bins = bins(main_peak_loc:end);
idx_half_cross = find(right_part <= half_max, 1, 'first');

% Get the threshold intensity
fwhm_threshold = right_bins(idx_half_cross)

T = fwhm_threshold * 255/100

BW = im2bw(I(:,:,2),T);
% figure;imshow(BW);
% title(['Binarized Image']);

% Invertir la imagen
BW = imcomplement(BW);
% figure;imshow(BW);
% title(['Binarized Image']);

%% Excentricidad
% Obtener las propiedades de las regiones conectadas
stats = regionprops(BW, 'Eccentricity', 'Area', 'Solidity', 'BoundingBox');

cuantas = size(stats);
excentricidad = [];
solidity =[];
area =[];
bounding_box =[];

for i=1:cuantas(1,1),
    excentricidad = [excentricidad stats(i).Eccentricity;];
    solidity = [solidity stats(i).Solidity;];
    area = [area stats(i).Area;];
    bounding_box= [bounding_box stats(i).BoundingBox;];
end

% figure; stem(excentricidad);
% title(['Excentricidad']);
% 
% 
% figure; stem(solidity);
% title(['Solidity']);
% 
% 
% figure; stem(area);
% title(['Area']);


%% Etiquetar y contar
L = bwlabel(BW);
numero_de_celulas = max(max(L));

Tamano = [];
stats = regionprops(L, 'Area');
Tamano = [stats.Area];


%% Erosionar
SE = strel('disk',3); %structural element: Disco de radio 3
BWE = imerode(BW,SE);

%% Dilatar
SE = strel('disk',3); %structural element: Disco de radio 3
BWED = imdilate(BWE,SE);


%% Erosionar
BWEDE = imerode(BWED,SE);

%% Etiquetar y contar
L = bwlabel(1-BWEDE);
numero_de_celulas = max(max(L));

%%Eliminar elementos pequeños
Tamano_Corte = 30000;

% figure;imshow(BW);

for i=1:numero_de_celulas,
    if Tamano_Corte>Tamano(1,i);
    % if Tamano_Corte>stats(i).Area;
        [aux1,aux2] = find(L==i);
        for j=1:size(aux1),
            BW(aux1(j),aux2(j)) = 0;
        end
    end
end

% figure;imshow(BW);

BW = imerode(BW,SE);

% figure;imshow(BW);
BW = imdilate(BW, SE);
% figure;imshow(BW);

BW = imerode(BW, SE);
% figure;imshow(BW);

L = bwlabel(BW);
numero_de_celulas = max(max(L))

Tamano = [];
stats = regionprops(L, 'Area', 'Eccentricity', 'Solidity', 'BoundingBox');
Tamano = [stats.Area];

% figure;imshow(BW);

L = bwlabel(BW);
% numero_de_celulas = max(max(L))

Tamano = [];
stats = regionprops(L, 'Area', 'Eccentricity', 'Solidity', 'BoundingBox');
Tamano = [stats.Area];

numero_de_celulas = length(stats);

cuantas = size(stats);
excentricidad = [];
solidity =[];
area =[];
bounding_box =[];

for i=1:cuantas(1,1),
    excentricidad = [excentricidad stats(i).Eccentricity;];
    solidity = [solidity stats(i).Solidity;];
    area = [area stats(i).Area;];
    bounding_box= [bounding_box stats(i).BoundingBox;];
end

figure; stem(area);
title(['Area']);

figure; stem(solidity);
title(['Solidity']);

figure; stem(excentricidad);
title(['Excentricity']);

% avg_area = mean(area)
% avg_eccentricity = mean(excentricidad)
% avg_solidity = mean(solidity)

% Add inputs
% Compute dynamic fuzzy membership parameters
range_area = [min(area), max(area)]
range_solidity = [min(solidity), max(solidity)]
range_eccentricity = [min(excentricidad), max(excentricidad)]

fis = mamfis('Name','PathologyDetection');

% --- Area ---
fis = addInput(fis, range_area, 'Name', 'Area');

fis = addMF(fis, 'Area', 'trapmf', ...
    [0 0 50 100], 'Name', 'Small');
fis = addMF(fis, 'Area', 'trapmf', ...
    [90 200 500 1000], 'Name', 'Medium');
fis = addMF(fis, 'Area', 'trapmf', ...
    [900 2500 3500 4000], 'Name', 'Large');
fis = addMF(fis, 'Area', 'trapmf', ...
    [3800 5000 7000 15000], 'Name', 'Very Large');

% --- Solidity ---
fis = addInput(fis, range_solidity, 'Name', 'Solidity');

fis = addMF(fis, 'Solidity', 'trapmf', ...
    [0 0 0.3 0.65], 'Name', 'Low');
fis = addMF(fis, 'Solidity', 'trapmf', ...
    [0.6 0.7 1 1], 'Name', 'High');

% --- Eccentricity ---
fis = addInput(fis, range_eccentricity, 'Name', 'Eccentricity');

fis = addMF(fis, 'Eccentricity', 'trapmf', ...
    [0 0 0.5 0.7], 'Name', 'Low');
fis = addMF(fis, 'Eccentricity', 'trapmf', ...
    [0.6 0.9 1 1], 'Name', 'High');

% --- Output ---
fis = addOutput(fis, [0 2], 'Name', 'Cell Type');
fis = addMF(fis, 'Cell Type', 'trimf', [0 0.2 0.3], 'Name', 'None');
fis = addMF(fis, 'Cell Type', 'trimf', [0.25 1 1.6], 'Name', 'Cells');
fis = addMF(fis, 'Cell Type', 'trimf', [1.6 2 2], 'Name', 'Pathological Cells');

ruleList = [ 
    %Area, Solidity, Excentricity, Cell Type ? ? 
    1 1 1 1 1 1;
    1 2 2 1 1 1;
    1 1 2 1 1 1;
    1 2 1 1 1 1;
    2 1 1 1 1 1;
    2 1 2 1 1 1;
    2 2 1 1 1 1;
    2 2 2 2 1 1;
    3 1 1 1 1 1;
    3 1 2 3 1 1;
    3 2 1 1 1 1;
    3 2 2 3 1 1;
    4 1 1 1 1 1;
    4 1 2 1 1 1;
    4 2 1 1 1 1;
    4 2 2 1 1 1;
    ];

fis = addRule(fis, ruleList);

cell_type_scores = zeros(numero_de_celulas, 1);

for i = 1:numero_de_celulas
    input = [area(i), solidity(i), excentricidad(i)];
    cell_type_scores(i) = evalfis(fis, input);
    % cell_type_scores(i)
end

cell_type_scores;

figure, imshow(I), title('Detected Cells by Type');
hold on;


cells_total =0;
cells_path=0;
for i = 1:numero_de_celulas
    val = round(cell_type_scores(i));
    if val == 1
        rectangle('Position', stats(i).BoundingBox, 'EdgeColor', 'r', 'LineWidth', 2); % Type A
        cells_total = cells_total+1;
    elseif val == 2
        rectangle('Position', stats(i).BoundingBox, 'EdgeColor', 'g', 'LineWidth', 2); % Type B
        cells_path = cells_path+1;
        cells_total = cells_total+1;
    % elseif val == 0
    %     rectangle('Position', stats(i).BoundingBox, 'EdgeColor', 'b', 'LineWidth', 2); % Type B
    end
end
hold off;


fprintf('Total cells: %i \n', cells_total );
fprintf('Pathological Cells: %i', cells_path);
