%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       Script para extrair e interpolar dados do MSWEP 1.2 para o MGB       %%%
%%%%%%%%%%%%             Sly Wongchuig Correa - IPH/UFRGS             %%%%%%%%%%%%
%%%%                                Junho 2018                                %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% versão 1.1 de JP
%%%%%%%%%%%% é posterior e menos elegante
%%%%%%% essa aqui eh pro MERGE

% 0. Precipitation
% 1. Wheight total
% 2. time
% 3. lat
% 4. lon

function chuvabin=IDW_MERGE(ultimo_dia,link,dist)

run 'D:\JP\ANA\extensao_chuva\matlab\MERGE_toolbox\nctoolbox-master\setup_nctoolbox.m'

date1 = datevec(ultimo_dia+1);
date2 = datevec(datenum(date)-1);  %MERGE

if datenum(date1)>datenum(date2)
    disp('A precipitacao ja esta atualizada!')
    chuvabin=[];
end
    
periodo=datenum(date1):1:datenum(date2); periodo=periodo';
fecha=datevec(periodo);
ano=fecha(:,1);mes=fecha(:,2);dia=fecha(:,3);


file=dir([link '*.grib2']);
xx = isempty(file);
if xx==1
    disp('There is no available data')
else
    if length(file)==1
    disp(['Available data just for ' num2str(file.name)]);
    elseif length(file)>1
    disp(['Available data from ' num2str(file(1,1).name) ' to ' num2str(file(end,1).name)]);
    end
end

%% Leer datos de estaciones a ser interpoladas
disp('Select mini.gtp or interpolated points');
temp=importdata('Mini.gtp');     % Cambiar
X=temp.data(:,3);
Y=temp.data(:,4);

% flag=input('Write 0 for fixed radius or 1 for number of neighbours = ');
% 
% if flag == 0
    
xflag=0;

% r2=input('Write factor of min distance = ');
r2=3;
e=2;

%% coordenadas do MERGE

nx=1001;
ny=924;
dxdy=0.1;
xini=-120.05;
yini=-60.05;

lat=yini+dxdy/2:dxdy:yini+ny*dxdy-dxdy/2;
lon=xini+dxdy/2:dxdy:xini+nx*dxdy-dxdy/2;

lat=lat';
lon=lon';
lat2=lat; lon2=lon;


%% Find extraction area
disp(['Write latitud coordinate between ' num2str(min(lat2)) ' to ' num2str(max(lat2))]);
latitudi=min(Y)-1;%input('Write minimum latitud (?) botton = ');
latitudf=max(Y)+1;%input('Write maximum latitud (?) top = ');

disp('-----------------------------------------------');

disp(['Write longitud coordinate between ' num2str(min(lon2)) ' to ' num2str(max(lon2))]);
longitudi=min(X)-1;%input('Write minimum longitud (?) left = ');
longitudf=max(X)+1;%input('Write maximum longitud (?) right = ');

%% Position within matrix
% Latitud
errolat=abs(latitudi+0.06-lat2);
xx=find(errolat<0.05);
if isempty(xx)
    xx=1;
end
lati=xx;

errolat=abs(latitudf-0.06-lat2);
xx=find(errolat<0.05);
if isempty(xx)
    xx=length(lat2);
end
latf=xx;

% Longitud
errolon=abs(longitudi-0.06-lon2);
xx=find(errolon<0.05);
if isempty(xx)
    xx=1;
end
loni=xx;

errolon=abs(longitudf+0.06-lon2);
xx=find(errolon<0.05);
if isempty(xx)
    xx=length(lon2);
end
lonf=xx;

% Extract data for study area
Longitud=lon2;%(loni:lonf,1);
Latitud=lat2;%(lati:latf,1);

% Creating longitud and latitud vectors
for j=1:length(Latitud)            
Lo(:,j)=Longitud;       
end
for k=1:length(Longitud)
La(:,k)=Latitud;
end
N(:,1)=reshape(Lo',[],1);
N(:,2)=reshape(La,[],1);


% %% Creamos los vectores de distancia
% mini_dist=zeros(length(X),1);
% r3=mini_dist;
% for i=1:length(X)
%         D=[]; 
%         D= sqrt((X(i)-N(:,1)).^2 +(Y(i)-N(:,2)).^2);
%         if min(D)==0
%             disp('Error: One or more stations have the coordinates of an interpolation point')
%             return
%         end
%         mini_dist(i)= min(D);
% end
% 
% % Para una minima distancia de radio
% if xflag==0
%     for i=1:length(X)
%         r3(i)=r2.*mini_dist(i);
%     end
% end
% 
% 
% %% JP- determinar as linhas de N que são menores que r3(i)
% 
% disp('tá acabando')
% for i=1:length(X)
%         D=[]; 
%         D= sqrt((X(i)-N(:,1)).^2 +(Y(i)-N(:,2)).^2);
%         if min(D)==0
%             disp('Error: One or more stations have the coordinates of an interpolation point')
%             return
%         end
%         dist(i).zzz= find(D<r3(i)); %posicao das celulas no vetor
%         dist(i).D=D(dist(i).zzz);%distancia das celulas ao centro da MB
% end

%% Leer datos MSWEP
% iT=1;
chuvabin=NaN(length(X),length(periodo));
for i=1:length(periodo)
    disp('------------------------------------------------');
    disp(['Processing day ' num2str(dia(i)) '/' num2str(mes(i)) '/' num2str(ano(i))]);

% if i>anoi % Ya fue leido anteriormente
    % Read data:
    fecha=datestr(datenum(ano(i),mes(i),dia(i)),'yyyymmdd');
    url=[link 'MERGE_CPTEC_' fecha '.grib2'] ;
    geo = ncgeodataset(url);
    nm=geo.geovariable('Precipitation_surface');
    k=nm(1,1:ny,1:nx);
    P=reshape(k,[924,1001]);
% 	subplot(1,2,1)
%     imagesc(P)
%     set(gca,'ydir','normal')
    prec=P(lati:latf,loni:lonf);
    P=prec(:);
    P(P<0)=0; % !! mudei pq tinha dias miseraveis NaN

    Pint=NaN(length(X),1);
    for k=1:length(X)
        wV =[];
        zzz=dist(k).zzz;
        Pc=P(zzz,:); D=dist(k).D;
        Ptemp = Pc.*(1./D.^e);
        wV = 1./D.^e;
        if isempty(D)
            Ptemp=NaN;
        else
            Ptemp=sum(Ptemp)./sum(wV);
        end
        Pint(k,:)=Ptemp;
    end
%     subplot(1,2,2)
%     imagesc(prec)
%     set(gca,'ydir','normal')
%     set(gcf,'position',[50 50 1400 700])
%     pause(0.4)
    chuvabin(:,i)=Pint;
    clear prec
    clear P

end

end
