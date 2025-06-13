%% Cargamos la imagen (variable I)

clear all
close all

% [I,map] = imread('capture_1548243990843.jpg');
% [I,map] = imread('capture_1548244948896.jpg');
% [I,map] = imread('capture_1548244759635.jpg');
% [I,map] = imread('capture_1548244554243.jpg');
% [I,map] = imread('capture_1548244406801.jpg');
% [I,map] = imread('capture_1548244339396.jpg');
[I,map] = imread('capture_1548244150176.jpg');


figure; 
imshow(I);
title(['Original Image']);

%% Vemos el nivel de intensidad de cada componente RGB
figure; imshow(I(:,:,1)); % Rojo
title(['Red Channel']);

figure; imshow(I(:,:,2)); % Verde
title(['Green Channel']);

figure; imshow(I(:,:,3)); % Azul
title(['Blue Channel']);
%% Representamos los histogramas de cada canal
figure; imhist(I(:,:,1)); % Rojo
title(['Red histogram']);

figure; imhist(I(:,:,2)); % Verde
title(['Green histogram']);

figure; imhist(I(:,:,3)); % Azul
title(['Blue histogram']);

%% Binariamos la imagen usando la informacion del canal verde
% [T,EM] = graythresh(I(:,:,2))

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
figure;imshow(BW);
title(['Binarized Image']);

% Invertir la imagen
BW = imcomplement(BW);
figure;imshow(BW);
title(['Binarized Image']);

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

figure; stem(excentricidad);
title(['Excentricidad']);


figure; stem(solidity);
title(['Solidity']);


figure; stem(area);
title(['Area']);


% Acceder a la excentricidad del primer objeto (ajusta el índice si tienes múltiples objetos)

% Mostrar la excentricidad

%% Etiquetar y contar
L = bwlabel(BW);
numero_de_celulas = max(max(L));

Tamano = [];
stats = regionprops(L, 'Area');
Tamano = [stats.Area];

figure; scatter(excentricidad,Tamano);
title(['Relación de características']);
xlabel('Excentricidad');
ylabel('Area');

figure; scatter(solidity,Tamano);
title(['Relación de características']);
xlabel('Solidity');
ylabel('Area');

%% Erosionar
SE = strel('disk',3); %structural element: Disco de radio 3
BWE = imerode(BW,SE);
figure;imshow(BWE);
title(['After removing small noise for the first time']); %expands white regions, merging nearby objects.
%L = bwlabel(BWEDE);
%max(max(L))

I2=I;
for i=1:numero_de_celulas,
    [aux3,aux4] = find(L==i);
    for j=1:size(aux3),
        I2(aux3(j),aux4(j),1) = 255;
        I2(aux3(j),aux4(j),2) = 0;
        I2(aux3(j),aux4(j),3) = 0;
    end
end
figure; imshow(I2)
title(['Cells located after removing small noise for the first time']);


%% Dilatar
SE = strel('disk',3); %structural element: Disco de radio 3
BWED = imdilate(BWE,SE);
figure;imshow(BWED);
title(['Dilated white regions']); %expands white regions, merging nearby objects.

I2=I;
for i=1:numero_de_celulas,
    [aux3,aux4] = find(L==i);
    for j=1:size(aux3),
        I2(aux3(j),aux4(j),1) = 255;
        I2(aux3(j),aux4(j),2) = 0;
        I2(aux3(j),aux4(j),3) = 0;
    end
end
figure; imshow(I2)
title(['Cells located after dilating white regions']);


%% Erosionar
BWEDE = imerode(BWED,SE); % removes small noise.
figure;imshow(BWEDE);
title(['Erosion of small noise second time']);
%L = bwlabel(BWEDE);
%max(max(L))

I2=I;
for i=1:numero_de_celulas,
    [aux3,aux4] = find(L==i);
    for j=1:size(aux3),
        I2(aux3(j),aux4(j),1) = 255;
        I2(aux3(j),aux4(j),2) = 0;
        I2(aux3(j),aux4(j),3) = 0;
    end
end
figure; imshow(I2)
title(['Cells located after removing small noise']);


%% Etiquetar y contar
L = bwlabel(1-BWEDE);
numero_de_celulas = max(max(L));
figure; imshow(1-BWEDE); 
title(['After counting (again)']);
figure; imshow(BWEDE); 
title(['After counting (again)']);

%%Eliminar elementos pequeños
Tamano_Corte = 700;

for i=1:numero_de_celulas,
    if Tamano_Corte>Tamano(1,i);
    % if Tamano_Corte>stats(i).Area;
        [aux1,aux2] = find(L==i);
        for j=1:size(aux1),
            BW(aux1(j),aux2(j)) = 0;
        end
    end
end

BW2=BW;
%% Eliminar los que tocan el borde y son mejores a 1000
Tamano_borde = 800;


%%Eliminar cuando la excentricidad es menor a 0,8
Tamano_Excentricity = 0.8;
Tamano_Superior_Cell= 800;
Solidity_max = 0.8;
% Etiquetar y contar
L = bwlabel(BW);
numero_de_celulas = max(max(L));

Tamano = [];
stats = regionprops(L, 'Area', 'Eccentricity', 'Solidity', 'BoundingBox');
Tamano = [stats.Area];

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

avg_area = mean(area)
avg_eccentricity = mean(excentricidad)
avg_solidity = mean(solidity)

%Poner a cero el valor si lo cumple
for i=1:numero_de_celulas,
    % if ((Tamano_Superior_Cell>Tamano(1,i))|(Tamano_Excentricity>excentricidad(1,i))|(Solidity_max>solidity(1,i)));
    %if (Tamano_Superior_Cell>Tamano(1,i))
    if (area(i) < 0.8 * avg_area)
    % if (area(i) < 0.8 * avg_area || solidity(i) < 0.9 * avg_solidity)

        [aux3, aux4] = find(L == i);
        [aux3,aux4] = find(L==i);
        for j=1:size(aux3),
            BW(aux3(j),aux4(j)) = 0;
        end
    end
end

figure; imshow(BW);
title(['After removing small things']);

for i=1:numero_de_celulas,
    if (excentricidad(i) < 0.8*avg_eccentricity)
        [aux3, aux4] = find(L == i);
        [aux3,aux4] = find(L==i);
        for j=1:size(aux3),
            BW(aux3(j),aux4(j)) = 0;
        end
    end
end 
figure; imshow(BW);
title(['After removing "not oval enough"']);

BW_failed = false(size(BW));  % Black canvas

for i=1:numero_de_celulas,
    if (solidity(i) < 0.8 * avg_solidity)
        [y, x] = find(L == i);
        for j = 1:length(y)
            BW_failed(y(j), x(j)) = 1;  % Paste into black image
        end
        [aux3, aux4] = find(L == i);
        [aux3,aux4] = find(L==i);
        for j=1:size(aux3),
            BW(aux3(j),aux4(j)) = 0;
        end
    end

end
figure; imshow(BW);
title(['After removing not solid enough']);


figure; imshow(BW_failed);
title('Failed cells before reprocessing');

BW_rescued = false(size(BW));

Tamano_Cell_Patologico_min=2000;
Tamano_Cell_Patologico_max=3000;
numero_de_celulas_patologicas_f=0;
posi_patologica_f=[];

L_failed = bwlabel(BW_failed);

numero_de_celulas_f = max(max(L));
Tamano = [];

for i=1:numero_de_celulas_f,
    aux = find(L_failed==i);
    aux2 = size(aux);
    Tamano = [Tamano aux2(1,1)];
end

figure; stem(Tamano);
%Poner a cero el valor si lo cumple
for i=1:numero_de_celulas_f,
    if ((Tamano_Cell_Patologico_min<Tamano(1,i)) & (Tamano_Cell_Patologico_max>Tamano(1,i)))
    % if ((Tamano_Cell_Patologico_min<stats(i).Area) & (Tamano_Cell_Patologico_max>stats(i).Area))
        numero_de_celulas_patologicas_f=numero_de_celulas_patologicas_f+1;
        posi_patologica_f=[posi_patologica_f i];
        BW_rescued(L_failed == i) = 1;
    end
end

figure; imshow(BW_rescued);
title('Rescued cells after size testing');

BW_failed_eroded = imerode(BW_failed, strel('disk', 4));
figure; imshow(BW_failed_eroded)
title('After eroding')

L_eroded = bwlabel(BW_failed_eroded);
stats_eroded = regionprops(L_eroded, 'Solidity', 'Area');
solidity_eroded = [stats_eroded.Solidity];

for k = 1:length(solidity_eroded)
    if solidity_eroded(k) >= 0.9 * avg_solidity
        BW_rescued(L_eroded == k) = 1;
    end
end

figure; imshow(BW_rescued);
title('Rescued cells after eroding and solidity test');

BW_final = BW | BW_rescued;


figure; imshow(BW_rescued);
title('Rescued cells after erosion and re-testing');

figure; imshow(BW_final);
title('Final image with rescued cells re-added');


numero_de_celulas = max(max(L));

L = bwlabel(BW_final);
numero_de_celulas = max(max(L));

Tamano_fin = [];
stats = regionprops(L, 'Area');
Tamano_fin = [stats.Area];

figure; stem(Tamano_fin);
title(['Size of cells in final pic'])

Tamano_smallsmall = 300;
for i=1:numero_de_celulas,
    if (Tamano_smallsmall>Tamano_fin(1,i))
        BW_final(L == i) = 0;
    end
end

figure; imshow(BW_final);
title('Final image with after removing small small things');

L = bwlabel(BW_final);

I2=I;
for i=1:numero_de_celulas,
    [aux3,aux4] = find(L==i);
    for j=1:size(aux3),
        I2(aux3(j),aux4(j),1) = 255;
        I2(aux3(j),aux4(j),2) = 0;
        I2(aux3(j),aux4(j),3) = 0;
    end
end
figure; imshow(I2)
title(['Cells located eliminating small and irregular ones?']);

%% Detectar cuantas son patológicas

Tamano_Cell_Patologico=2500;
numero_de_celulas_patologicas=0;
posi_patologica=[];

L = bwlabel(BW_final);
numero_de_celulas = max(max(L))

Tamano = [];

stats = regionprops(L, 'Area');
Tamano = [stats.Area];

%Poner a cero el valor si lo cumple
for i=1:numero_de_celulas,
    if (Tamano_Cell_Patologico<Tamano(1,i))
    % if (Tamano_Cell_Patologico<stats(i).Area)
        numero_de_celulas_patologicas=numero_de_celulas_patologicas+1;
        posi_patologica=[posi_patologica i];
    end
end

% Etiquetar y contar
numero_de_celulas_patologicas;

%%Pintar la imagen original e indicar número de células

fprintf('\n El número de glóbulos rojos es %.f, de las que células patológicas hay %.f\n', numero_de_celulas, numero_de_celulas_patologicas);

I2=I;
for i=1:numero_de_celulas,
    [aux3,aux4] = find(L==i);
    for j=1:size(aux3),
        I2(aux3(j),aux4(j),1) = 255;
        I2(aux3(j),aux4(j),2) = 0;
        I2(aux3(j),aux4(j),3) = 0;
    end
end

for i=1:numero_de_celulas_patologicas,
    [aux3,aux4] = find(L==posi_patologica(i));
    for j=1:size(aux3),
        I2(aux3(j),aux4(j),1) = 0;
        I2(aux3(j),aux4(j),2) = 255;
        I2(aux3(j),aux4(j),3) = 0;
    end
end
figure; imshow(I2)
title(['Normal cells and pathological cells located !!! '])

NombreFichero=['capture_1548243990843_Glob_Rojo_',int2str(numero_de_celulas),'_Pat_',int2str(numero_de_celulas_patologicas),'.jpg'];
% imwrite(I2, NombreFichero);
% Guardar la imagen como un archivo JPEG con alta calidad
%imwrite(I2, NombreFichero, 'Quality', 100);
