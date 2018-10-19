display('Inserte el nombre del video con formato incluido y luego el tiempo de lectura.');
name = input('Nombre del video: ','s');
time = input('Tiempo [s]: ');       %Tramo del video a utilizar
display('procesando...')
v = VideoReader(name);              %Importar video
numFrames = get(v,'NumberOfFrames');    %Ajustar tiempo a largo de video
t = round(120*time); %Frames a trabajar
t = min(t,numFrames); %No pasarse del largo del video
VidFrame = read(v,[1 t]);           %Leer video y dejarlo en tama�o deseado
VidFrame = imresize(VidFrame,[540 960]);
Im = zeros(540,960,t,'uint8');
figure('Name','Marcadores a utilizar')      %Clickear n marcadores y apretar Enter
imshow(VidFrame(:,:,:,1));
title('Seleccione los marcadores en el orden correspondiente')
click = ginput;
if isempty(click);
    return
end
close figure 1
figure('Name','Horizontal a utilizar')      %Pedir horizontal en caso de ser necesaria
imshow(VidFrame(:,:,:,1));
title('Marque el plano horizontal que se desea utilizar ("enter" en caso de no necesitarlo)')
horizontal = ginput(2);
close figure 1
%Crear matriz de centroides
c = zeros(length(click(:,1)),2,t);
c(:,:,1) = click;
Pred = click;
nummar = length(c(:,1,1));
%Dimensiones de pantalla
scrsz = get(groot,'ScreenSize');
%PROCESAMIENTO DEL VIDEO
%Tama�o caja para prediccion
l = 25;
%Procesar video
for i = 1:t
%Transformar video a escala de grises (3 matrices a 1 matriz)
VidFrame2 = rgb2gray(VidFrame(:,:,:,i));
%Colocar solo la zona del click en la matriz de zeros
for j = 1:nummar
    c = round(c);
    if i == 1;
        Im(c(j,2,i)-l:c(j,2,i)+l,c(j,1,i)-l:c(j,1,i)+l,i) = VidFrame2(c(j,2,i)-l:c(j,2,i)+l,c(j,1,i)-l:c(j,1,i)+l);
    else    
        Im(c(j,2,i-1)-l:c(j,2,i-1)+l,c(j,1,i-1)-l:c(j,1,i-1)+l,i) = VidFrame2(c(j,2,i-1)-l:c(j,2,i-1)+l,c(j,1,i-1)-l:c(j,1,i-1)+l);
    end           
end
%B�sicamente filtrar imperfecciones
VidFrame3 = medfilt2(Im(:,:,i), [1 1]);
%Convertir escala de grises a binario
VidFrame4 = im2bw(VidFrame3,0.85); %95
%Quita objetos de menos de x pixeles
VidFrame5 = bwareaopen(VidFrame4,10); %35
%Detectar objetos
bw = logical(VidFrame5);
%Obtener los centros
stats = regionprops(bw,'BoundingBox','Centroid');
centro = cat(1,stats.Centroid);
if i > 2 ;
    Pred = 2*c(:,:,i-1)-c(:,:,i-2);
end
for k = 1:nummar
    d = sqrt((centro(:,1)-Pred(k,1)).^2 + (centro(:,2)-Pred(k,2)).^2);
    num = find(d==min(d));
    c(k,:,i) = centro(num,:);
end
end
%Obtener �ngulos
%--------------------------
display('Escriba el n�mero correspondiente al �ngulo deseado.')
display('A continuaci�n presione el marcador asociado.')
display('Una vez seleccionados todos los �ngulos presione "Enter".')
display('Plano lateral:')
display([' 1. Tronco ';' 2. Rodilla';' 3. Tobillo'])
display('Plano frontal:')
display(' 4. Valgo din�mico')
display('Plano trasero:')
display(' 5. Tilt p�lvico')
display('Si no requiere cotas:')
display(' 6. Otro')
display('Ingrese vertice del �ngulo ("Enter" si ya termin�)')
%Valores de cotas
Nombreang = {'Tronco','Rodilla','Tobillo','Valgo din�mico','Tilt p�lvico',''};
CotasRuta = [45 45 45 45;65 65 145 155;70 80 95 105;5 5 7 7;0 0 0 0];
CotasMB = [47 47 53 53;65 65 145 155;70 80 95 100;5 5 7 7; 0 0 0 0];
%Crear ciclo mientras ingresa �ngulos
ingresando = 1;
nAng = 0;
while ingresando == 1
    anguloin = input('');
    if isempty(anguloin)
        break
    end
    nAng = nAng+1;
    ncotas(nAng) = anguloin; %#ok<SAGROW>
    figure('Name','Selecci�n de �ngulos') 
    imshow(VidFrame(:,:,:,1));
    title('Seleccione el marcador correspondiente al v�rtice del �ngulo')
    ang(nAng,:) = ginput(1); %#ok<SAGROW>
    close all
end
if ~exist('ang','var');
    %Vector de colores
    color = 'ygbmkcwygbmkcw';
    %Mostrar recorrido
    cp = permute(c,[3 2 1]);
    figure('Name','Trayectoria de los marcadores')
    imshow(VidFrame(:,:,:,1))
    set(figure(1),'Position',[20 scrsz(4)/2-100 scrsz(3)/2 scrsz(4)/2])
    hold on
    for j = 0:nummar-1
        plot(cp(:,1,nummar-j),cp(:,2,nummar-j),color(nummar-j))
    end
    hold off
    return
end
Pto = zeros(length(ang(:,1)),1);
%close figure 1
%Asociar centro a marcador
for k = 1:nAng
    d = sqrt((c(:,1,1)-ang(k,1)).^2 + (c(:,2,1)-ang(k,2)).^2);
    num = find(d==min(d));
    Pto(k) = min(num) ;
end
%Definir funcion para insertar valor
insert = @(a, x, n)cat(1,  x(1:n,:), a, x(n+1:end,:));
%CASOS PARA CALCULO DE ANGULOS --------------
%Caso Torso Lateral
if any(ncotas==1) || any(ncotas==5)
    if isempty(horizontal)
        display('ERROR: Horizontal no entregada')
        return
    end
    if any(ncotas==1)
        ind = find(ncotas==1);
    else
        ind = find(ncotas==5);
    end    
    %Calcular vector horizontal
    mhor = (horizontal(4)-horizontal(3)) / (horizontal(2)-horizontal(1));
    chor1 = [c(Pto(ind),1,:)+50 c(Pto(ind),2,:)+50*mhor];
    chor2 = [c(Pto(ind),1,:)-50 c(Pto(ind),2,:)-50*mhor];
    tam = size(c);
    moly1 = zeros(tam(1)+1,tam(2),tam(3));
    moly2 = zeros(tam(1)+2,tam(2),tam(3));
    moly3 = zeros(tam(1)+3,tam(2),tam(3));
    for ti = 1:t
        moly1(:,:,ti) = insert(chor1(:,:,ti),c(:,:,ti),Pto(ind));
        moly2(:,:,ti) = insert(chor2(:,:,ti),moly1(:,:,ti),Pto(ind)+1);
        moly3(:,:,ti) = insert(c(Pto(ind),:,ti),moly2(:,:,ti),Pto(ind)+2);
    end
    c = moly3;
    nummar = nummar + 3; 
    %Calcular modulo de los vectores
    a = [c(Pto(1:ind)-1,1,:)-c(Pto(1:ind),1,:) c(Pto(1:ind)-1,2,:)-c(Pto(1:ind),2,:);
         c(Pto(ind+1:end)+2,1,:)-c(Pto(ind+1:end)+3,1,:) c(Pto(ind+1:end)+2,2,:)-c(Pto(ind+1:end)+3,2,:)];
    b = [c(Pto(1:ind)+1,1,:)-c(Pto(1:ind),1,:) c(Pto(1:ind)+1,2,:)-c(Pto(1:ind),2,:);
         c(Pto(ind+1:end)+4,1,:)-c(Pto(ind+1:end)+3,1,:) c(Pto(ind+1:end)+4,2,:)-c(Pto(ind+1:end)+3,2,:)];
    distA = sqrt(a(:,1,:).^2 + a(:,2,:).^2);
    distB = sqrt(b(:,1,:).^2 + b(:,2,:).^2);
else
    %Calcular modulo de los vectores
    a = [c(Pto-1,1,:)-c(Pto,1,:) c(Pto-1,2,:)-c(Pto,2,:)];
    b = [c(Pto+1,1,:)-c(Pto,1,:) c(Pto+1,2,:)-c(Pto,2,:)];
    distA = sqrt((c(Pto-1,1,:)-c(Pto,1,:)).^2 + (c(Pto-1,2,:)-c(Pto,2,:)).^2);
    distB = sqrt((c(Pto+1,1,:)-c(Pto,1,:)).^2 + (c(Pto+1,2,:)-c(Pto,2,:)).^2);
end
%Calcular angulo entre vectores 
dots = dot(a,b,2);
ppunto = permute(dots,[3 1 2]);
dist = permute((distA(:,1,:).*distB(:,1,:)),[3 1 2]);
angulo = acos(ppunto./dist);
angulogr = angulo*360/(2*pi);
if any(ncotas==4)
    ind = find(ncotas==4);
    angulogr(:,ind) = 180-angulogr(:,ind);
elseif any(ncotas==5)
    ind = find(ncotas==5);
    if angulogr(:,ind) > 90
        angulogr(:,ind) = 180-angulogr(:,ind);
    end
else
    if any(ncotas==1)
        ind = find(ncotas==1);
        for g = 1:t
            if angulogr(g,ind) > 90
                angulogr(g,ind) = angulogr(g,ind)-90;
            end
        end
    end
end
% --------------------------------------------
%MOSTRAR IMAGENES
%Vector de colores
color = 'ygbmkcwygbmkcw';
%Mostrar recorrido
cp = permute(c,[3 2 1]);
figure('Name','Trayectoria de los marcadores')
imshow(VidFrame(:,:,:,1))
set(figure(1),'Position',[20 scrsz(4)/2-100 scrsz(3)/2 scrsz(4)/2])
hold on
for j = 0:nummar-1
    plot(cp(:,1,nummar-j),cp(:,2,nummar-j),color(nummar-j))
end
hold off
%Mostrar angulos
figure('Name','Gr�ficos de �ngulos vs tiempo')
for k = 1:nAng
    minang = min(angulogr(:,k));
    maxang = max(angulogr(:,k));
    xmin = find(angulogr(:,k)==minang,1);
    xmax = find(angulogr(:,k)==maxang,1);
    subplot(nAng,1,k)
    hold on
    plot(0:1/120:time-1/120,angulogr(:,k),color(Pto(k)));
    plot((xmin-1)/120,minang,'*k')
    plot((xmax-1)/120,maxang,'*k')
    if ncotas(k) ~= 6
        pl1 = plot([0 time-1/120],[CotasRuta(ncotas(k),:)' CotasRuta(ncotas(k),:)'],'r');
        pl2 = plot([0 time-1/120],[CotasMB(ncotas(k),:)' CotasMB(ncotas(k),:)'],'c--');
    end
    hold off
    title(Nombreang(ncotas(k)))
    text((xmin + 5)/120,minang,['M�n = ' num2str(minang) '�']);
    text((xmax + 5)/120,maxang,['M�x = ' num2str(maxang) '�']);
    grid on
    xlabel('Tiempo [s]')
    ylabel('\theta [grados]')
end
if any(ncotas~=6)
    legend([pl1(1) pl2(1)],'Ruta','Mountain Bike','Location','SouthEast')
end
set(figure(2),'Position',[scrsz(3)/2+30 60 scrsz(3)/2-50 scrsz(4)-200])
%Mostrar video con marcadores
figure('Name','Cinem�tica de los marcadores')
%set(figure(3),'Position',[20 60 scrsz(3)/2-350 scrsz(4)/2-250])
for i = 1:t;
    imshow(VidFrame(:,:,:,i))
    hold on
    plot(c(:,1,i),c(:,2,i),'r-o');
    hold off
    F(i) = getframe(gcf) ;
    drawnow
end
display('Listo!')
return
