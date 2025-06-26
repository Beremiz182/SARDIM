% SARDIM - South America River Discharge Monitor
% MGB-SA - Modelo Hidrológico de Grandes Bacias (South America)
% HGE - Hidrologia de Grande Escala
% IPH - Instituto de Pesquisas Hidráulicas (UFRGS)
% Versão Matlab
% Processar automaticamente o arquivo de saída de vazão do MGB em variáveis 
% utilizadas nos atributos do shapefile da plataforma SARDIM.

% Limpar todas as variáveis
clear

% Número de minibacias
nMini = 33749;


%% JUNTANDO QTUDO

fid2=fopen('..\MGB_SA\Output_MSWEP\QTUDO.mgb');
Q2=fread(fid2,'single');
nT=length(Q2)/nMini;
Q=reshape(Q2,[nMini,nT]);

fid2=fopen('..\MGB_SA\Output\QTUDO.mgb');
Q2=fread(fid2,'single');
nT2=length(Q2)/nMini;
Q2=reshape(Q2,[nMini,nT2]);

Q3=[Q Q2(:,nT+1:end)];
nT2=size(Q3,2);
fclose all;

fid2=fopen('..\MGB_SA\Output\QTUDO.mgb','w+');
fwrite(fid2,Q3,'single');
fclose all;
disp('Dados de vazão lidos')


%% ESTATISTICA RODRIGO

% Data inicial dos dados de chuva
InitialDate = datenum('1-1-1979','mm-dd-yyyy');

% Número de intervalos de tempo
NT = load('NT.sim');
%%%NT = 15493;

% Data do último dia de atualização do MGB-SA
EndDate = InitialDate+NT-1;

% Vetor com datas em formato numérico
dates = (InitialDate:EndDate)';

% Gravar dados de vazão em uma variável matricial
flows = Q3;

% Calcular vazão atual (Qhoje)
Qhoje = flows(:,end);

% Calcular vazão média semanal dos últimos 7 dias (Semana_1)
Semana_1 = mean(flows(:,(end-6):end),2);

% Calcular vazão média semanal da semana anterior (Semana_2)
Semana_2 = mean(flows(:,(end-13):(end-7)),2);

% Calcular vazão média semanal 3 (Semana_3)
Semana_3 = mean(flows(:,(end-20):(end-14)),2);

% Calcular vazão média semanal 4 (Semana_4)
Semana_4 = mean(flows(:,(end-27):(end-21)),2);

% Calcular vazão média semanal 5 (Semana_5)
Semana_5 = mean(flows(:,(end-34):(end-28)),2);

% Calcular vazão média semanal 6 (Semana_6)
Semana_6 = mean(flows(:,(end-41):(end-35)),2);

% Calcular vazão média semanal 7 (Semana_7)
Semana_7 = mean(flows(:,(end-48):(end-42)),2);

% Calcular vazão média semanal 8 (Semana_8)
Semana_8 = mean(flows(:,(end-55):(end-49)),2);

% Vazoes dos últimos 3 meses
Qrecente = flows(:,(end-89):end);

% Vazoes dos últimos 30 dias
Qrecente30 = flows(:,(end-29):end);

flows = flows';
Qrecente = Qrecente';   

% Remover o primeiro ano da analise:
Date = datenum('1-1-1980','mm-dd-yyyy');
ii = find(Date == dates);
flows = flows(ii:end,:);
dates = dates(ii:end);

% Calcular a curva de duracao para todas as minibacias
quantiles = [0:100]';	
fdc = ones(length(quantiles),nMini);
fdc2 = ones(nMini,length(quantiles),12);

dates_vec=datevec(dates);
dates_vec=dates_vec(:,1:3);

for j=1:nMini
  
  %for i=1:length(quantiles)
      fdc(:,j) = prctile(flows(:,j),quantiles);
      for im=1:12
          g=find(dates_vec(:,2)==im);
          x= prctile(flows(g,j),quantiles);
          x= flipud(x);
          fdc2(j,:,im) =x;
      end
  %end
end

fdc = flipud(fdc);
fdc = fdc';

% Calcular Q90 
Q90 = fdc(:,91);
Q902 = fdc2(:,91,:);
Q902 = reshape(Q902,[nMini,12]);

Perm=zeros(nMini,30);
PermanenciaMax=zeros(nMini,1);
PermanenciaMin=zeros(nMini,1);
Permanencia=zeros(nMini,1); 

% Interpolar valores de permanencia para as vazões dos últimos 30 dias
disp('Interpolar valores de permanencia para as vazões dos últimos 30 dias')
for iMini = 1:nMini
  [cl,ia1]=unique(fdc(iMini,:));
  Perm(iMini,:) = interp1(fdc(iMini,ia1)',quantiles(ia1),Qrecente30(iMini,:));
  
  % Permanencia máxima dos últimos 30 dias
  PermanenciaMax(iMini) = max(Perm(iMini,:));
  
  % Permanência mínima nos últimos 30 dias
  PermanenciaMin(iMini) = min(Perm(iMini,:));
  
  % Permanencia atual
  Permanencia(iMini) = Perm(iMini,end);
  Permanencia(iMini) = min([max([Permanencia(iMini),0.1]),99.9]);
end

PermanenciaMax = PermanenciaMax';
PermanenciaMin = PermanenciaMin';

% Permanência máxima ou mínima nos últimos 30 dias
disp('Permanência máxima ou mínima nos últimos 30 dias')
Permanencia_30 = zeros(nMini,1);
for iMini = 1:nMini
  
  if (PermanenciaMax(iMini) >= 50 && PermanenciaMin(iMini) <= 50) 
      if (PermanenciaMax(iMini)-50) > (50-PermanenciaMin(iMini))
          Permanencia_30(iMini) = PermanenciaMax(iMini);
      else
          Permanencia_30(iMini) = PermanenciaMin(iMini);
      end
  elseif (PermanenciaMax(iMini) >= 50 && PermanenciaMin(iMini) >= 50)
      if (PermanenciaMax(iMini)-50) > (PermanenciaMin(iMini)-50)
          Permanencia_30(iMini) = PermanenciaMax(iMini);
      else
          Permanencia_30(iMini) = PermanenciaMin(iMini);
      end
  elseif (PermanenciaMax(iMini) <= 50 && PermanenciaMin(iMini) <= 50)
      if (50-PermanenciaMax(iMini)) > (50-PermanenciaMin(iMini))
          Permanencia_30(iMini) = PermanenciaMax(iMini);
      else
          Permanencia_30(iMini) = PermanenciaMin(iMini);
      end
  else
      Permanencia_30(iMini) = 50;
  end
  Permanencia_30(iMini) = min(max(Permanencia_30(iMini),0.1),99.9);
end

% Verificar qual curva de permanencia será utilizada na interpolação
DataAtual = datevec(datestr(now));
MesAtual = DataAtual(2);

for i=1:12
    
    if (MesAtual == i)
        fdc2 = fdc2(:,:,i);
        Q902 = Q902(:,i);
    end    
end

Perm2=zeros(nMini,30);
PermanenciaMax2=zeros(nMini,1);
PermanenciaMin2=zeros(nMini,1);
Permanencia2=zeros(nMini,1); 

% Interpolar valores de permanencia para as vazões dos últimos 30 dias
disp('Interpolar valores de permanencia sazonal para as vazões dos últimos 30 dias')
for iMini = 1:nMini
  [cl,ia1]=unique(fdc2(iMini,:));
  % mudanca JP 
  % pq fdc2 tem uma vazao constante para toda curva de permanencia
  if ia1==1
      Perm2(iMini,:)=50;
  else       
      Perm2(iMini,:) = interp1(fdc2(iMini,ia1)',quantiles(ia1),Qrecente30(iMini,:));
  end
  
  % Permanencia máxima dos últimos 30 dias
  PermanenciaMax2(iMini) = max(Perm2(iMini,:));
  
  % Permanência mínima nos últimos 30 dias
  PermanenciaMin2(iMini) = min(Perm2(iMini,:));
  
  % Permanencia atual
  Permanencia2(iMini) = Perm2(iMini,end);
  Permanencia2(iMini) = min([max([Permanencia2(iMini),0.1]),99.9]);
end

PermanenciaMax2 = PermanenciaMax2';
PermanenciaMin2 = PermanenciaMin2';

% Permanência máxima ou mínima nos últimos 30 dias
disp('Permanência sazonal máxima ou mínima nos últimos 30 dias')
Permanencia2_30 = zeros(iMini,1);
for iMini = 1:nMini
  
  if (PermanenciaMax2(iMini) >= 50 && PermanenciaMin2(iMini) <= 50) 
      if (PermanenciaMax2(iMini)-50) > (50-PermanenciaMin2(iMini))
          Permanencia2_30(iMini) = PermanenciaMax2(iMini);
      else
          Permanencia2_30(iMini) = PermanenciaMin2(iMini);
      end
  elseif (PermanenciaMax2(iMini) >= 50 && PermanenciaMin2(iMini) >= 50)
      if (PermanenciaMax2(iMini)-50) > (PermanenciaMin2(iMini)-50)
          Permanencia2_30(iMini) = PermanenciaMax2(iMini);
      else
          Permanencia2_30(iMini) = PermanenciaMin2(iMini);
      end
  elseif (PermanenciaMax2(iMini) <= 50 && PermanenciaMin2(iMini) <= 50)
      if (50-PermanenciaMax2(iMini)) > (50-PermanenciaMin2(iMini))
          Permanencia2_30(iMini) = PermanenciaMax2(iMini);
      else
          Permanencia2_30(iMini) = PermanenciaMin2(iMini);
      end
  else
      Permanencia2_30(iMini) = 50;
  end
  Permanencia2_30(iMini) = min(max(Permanencia2_30(iMini),0.1),99.9);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculos TR:

% Seleciona máximos anuais:
datesVEC = datevec(dates);
ano = datesVEC(:,1);
anos = (min(ano):max(ano))';
n = max(anos)-min(anos)+1;  

% Pensar em como usar o ano hidrologico no lugar do calendário
% Nao tem nenhuma verificacao se o ano está completo

% JP - incluindo as maiores vazões dos anos
Qmax2=nan(nMini,n-1);
Qmin2=Qmax2;

for iMini = 1:nMini
    for i = 1:n-1 % Excluir o ano atual
        STAT{iMini}.Qmax(i, 1) = max(flows(ano == anos(i), iMini));
        STAT{iMini}.Qmin(i, 1) = min(flows(ano == anos(i), iMini));
    end
    % JP - incluindo as maiores vazões dos anos
    Qmax2(iMini,:)=STAT{iMini}.Qmax(:);
    Qmin2(iMini,:)=STAT{iMini}.Qmin(:);
end

    
% Ajusta distribuições estatísticas:
% Máximas:
disp('Ajusta distribuições estatísticas: Máximas')
for iMini = 1:nMini
    
    STAT{iMini}.Qmax_emp = sort(STAT{iMini}.Qmax ,'descend');
    TR_emp = 1./((1:n-1)'/n);
    STAT{iMini}.TR_emp=TR_emp;
    
    % Ajuste parametros gumbel:
    % Parâmetros:
    x = mean(STAT{iMini}.Qmax);
    s = std(STAT{iMini}.Qmax);
    alpha = 0.7797*s;
    beta = x-0.5772*alpha;
    STAT{iMini}.x = x;
    STAT{iMini}.s = s;
    STAT{iMini}.alpha = alpha;
    STAT{iMini}.beta = beta;
    
    % Vazoes máximas:
    STAT{iMini}.QmaxGumbel = beta-alpha*log(-log(1-1./TR_emp)); 
    
    % Cálcula erros do ajuste:
    erro = STAT{iMini}.QmaxGumbel-STAT{iMini}.Qmax_emp; 
    
    % Erro médio quadrático percentual    
    STAT{iMini}.erroMax = 100*mean(erro.^2)^0.5 / mean(STAT{iMini}.Qmax_emp); 
    
    % Se ajuste é ruim, nao calcula TRs:
    erroCriticoCheia = 50;
    if STAT{iMini}.erroMax > erroCriticoCheia
      STAT{iMini}.TR_max = -999999;
    else
      % Série temporal de TRs:
      STAT{iMini}.TR_max = 1./(1-exp(-exp(-(Qrecente(:,iMini)-beta)/alpha)));
    end
end

% TR dos últimos 3 meses.
% Colocar tudo na mesma matriz:
for iMini = 1:nMini
    
    TRrecenteMax(:,iMini) = STAT{iMini}.TR_max;
    erroMAX(iMini) = STAT{iMini}.erroMax;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mínimas:
disp('Ajusta distribuições estatísticas: Mínimas')
for iMini = 1:nMini
  
  STAT{iMini}.Qmin_emp = sort(STAT{iMini}.Qmin);
  TR_emp = STAT{iMini}.TR_emp;
  STAT{iMini}.erroMin = -1;
  Qmin = STAT{iMini}.Qmin;
    
  if min(Qmin)>0    
    % Ajuste parametros weibull:
    % Parâmetros:
    xmin = mean(Qmin);
    smin = std(Qmin);
    G = skewness(Qmin);
    
    % A equacao de lambda é valida para G entre -1 e 2. Vamos limitar o
    % G.
    G=min(max(-1,G),2);
    
    H0 = 0.2777757913;
    H1 = 0.3132617714;
    H2 = 0.0575670910;
    H3 = -0.0013038566;
    H4 = -0.0081523408 ;

    lambda=1/(H0+H1*G+H2*G^2+H3*G^3+H4*G^4);
      
    if gamma(1+1/lambda)^2 >= gamma(1+2/lambda)
      STAT{iMini}.TR_min = -999999;
      STAT{iMini}.QminWeibull = TR_emp*0.0;
    else
      B=(gamma(1+2/lambda)-gamma(1+1/lambda)^2 )^-0.5;
      A=(1-gamma(1+1/lambda))*B;  
      
      STAT{iMini}.xmin = xmin;
      STAT{iMini}.smin = smin;
      STAT{iMini}.G = G;
      STAT{iMini}.lambda = lambda;
      STAT{iMini}.A = A;
      STAT{iMini}.B = B;

      % Vazoes minimas:
      K = A+B*((-log( 1-1./TR_emp)).^(1/lambda)-1);
      STAT{iMini}.QminWeibull = xmin+K*smin;     
      
      % Cálcula erros do ajuste:
      erro = STAT{iMini}.QminWeibull-STAT{iMini}.Qmin_emp;
      % Erro médio quadrático percentual  
      STAT{iMini}.erroMin = 100*mean(erro.^2)^0.5 / mean(STAT{iMini}.Qmin_emp);
      % Se ajuste é ruim, nao calcula TRs:
      erroCriticoSeca = 50;
      if STAT{iMini}.erroMin > erroCriticoSeca 
        STAT{iMini}.TR_min = -999999;
        STAT{iMini}.QminWeibull = TR_emp*0.0;
      else
        % Série temporal de TRs:
        Kmin = (Qrecente(:,iMini)-xmin)/smin;
        STAT{iMini}.TR_min =  1./(1- exp(-(max((Kmin-A)/B + 1, 0.0)).^lambda));
      end
    end
    else
    STAT{iMini}.TR_min = -999999;
    STAT{iMini}.QminWeibull = TR_emp*0.0;
    end
end

% TR dos últimos 3 meses.
% Colocar tudo na mesma matriz:
for iMini = 1:nMini
  
  TRrecenteMin(:,iMini) = STAT{iMini}.TR_min;
  QminWeibull(:,iMini) = STAT{iMini}.QminWeibull;
  erroMIN(iMini) = STAT{iMini}.erroMin;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Guarda o maior tempo de retorno. Positivo para cheias e negativo para
% estiagem.
TRrecenteMax = TRrecenteMax';
TRrecenteMin = TRrecenteMin';

for iMini = 1:nMini
  
  for iT = 1:90
    if TRrecenteMin(iMini, iT) == -999999
      TRrecenteMin(iMini, iT) = -999999;
    else
      TRrecenteMin(iMini, iT) = -TRrecenteMin(iMini, iT);
     end
  end
end

for iMini = 1:nMini
  
  for iT = 1:90
    if TRrecenteMax(iMini, iT) > abs(TRrecenteMin(iMini, iT))
      TR(iMini, iT) = TRrecenteMax(iMini, iT);
    elseif TRrecenteMax(iMini, iT) < abs(TRrecenteMin(iMini, iT))
      TR(iMini, iT) = TRrecenteMin(iMini, iT);
    else 
      TR(iMini, iT) = 1;
    end
  end
end

TR30 = TR(:,(end-29):end);
Tr = TR30(:,end);

% Calcular tempo de retorno máximo ou mínimo do último mês para cada minibacia
for iMini = 1:nMini
  TRMAX(iMini) = max(TR30(iMini,:));
end
TRMAX = TRMAX';
for iMini = 1:nMini
  TRMIN(iMini) = min(TR30(iMini,:));
end
TRMIN = TRMIN';
for iMini = 1:nMini
    aux=sign(TRMAX(iMini)+TRMIN(iMini));
    if aux==0
        aux=1;
    end
  Tr_30(iMini) = max(TRMAX(iMini),-TRMIN(iMini))*sign(TRMAX(iMini)+TRMIN(iMini));
end
Tr_30 = Tr_30';



%% Escrevendo a Permanencia para todos os meses - JP

anoss=unique(dates_vec(:,1));
Qpermhist=zeros(nMini,1);
str_meses=[];
Permanencia_300=Permanencia_30;
% loop anos
for yy=1:length(anoss)
    mfim=12;
    if yy==length(anoss)
        mfim=dates_vec(end,2)-1;
    end

    % loop meses
    for mm=1:mfim
        
        xx=find(dates_vec(:,1)==anoss(yy) & dates_vec(:,2)==mm);
        if mm>9
            str_meses=[str_meses ', Y' num2str(anoss(yy)) '-' num2str(mm)];
        else
            str_meses=[str_meses ', Y' num2str(anoss(yy)) '-0' num2str(mm)];
        end
        disp([num2str(anoss(yy)) '-' num2str(mm)])
        
        Perm=zeros(nMini,length(xx));
        PermanenciaMax=zeros(nMini,1);
        PermanenciaMin=zeros(nMini,1);
             
        Qrecente30=flows(xx,:);
        Qrecente30=Qrecente30';
        
        % Interpolar valores de permanencia para as vazões dos últimos 30 dias
        %disp('Interpolar valores de permanencia para as vazões dos últimos 30 dias')
        for iMini = 1:nMini
            
            [cl,ia1]=unique(fdc(iMini,:));
            Perm(iMini,:) = interp1(fdc(iMini,ia1)',quantiles(ia1),Qrecente30(iMini,:));
            
            % Permanencia máxima dos últimos 30 dias
            PermanenciaMax(iMini) = max(Perm(iMini,:));
            
            % Permanência mínima nos últimos 30 dias
            PermanenciaMin(iMini) = min(Perm(iMini,:));
            
        end
        
        PermanenciaMax = PermanenciaMax';
        PermanenciaMin = PermanenciaMin';      
        Permanencia_30 = zeros(iMini,1);
        
        for iMini = 1:nMini
            
            if (PermanenciaMax(iMini) >= 50 && PermanenciaMin(iMini) <= 50)
                if (PermanenciaMax(iMini)-50) > (50-PermanenciaMin(iMini))
                    Permanencia_30(iMini) = PermanenciaMax(iMini);
                else
                    Permanencia_30(iMini) = PermanenciaMin(iMini);
                end
            elseif (PermanenciaMax(iMini) >= 50 && PermanenciaMin(iMini) >= 50)
                if (PermanenciaMax(iMini)-50) > (PermanenciaMin(iMini)-50)
                    Permanencia_30(iMini) = PermanenciaMax(iMini);
                else
                    Permanencia_30(iMini) = PermanenciaMin(iMini);
                end
            elseif (PermanenciaMax(iMini) <= 50 && PermanenciaMin(iMini) <= 50)
                if (50-PermanenciaMax(iMini)) > (50-PermanenciaMin(iMini))
                    Permanencia_30(iMini) = PermanenciaMax(iMini);
                else
                    Permanencia_30(iMini) = PermanenciaMin(iMini);
                end
            else
                Permanencia_30(iMini) = 50;
            end
            Permanencia_30(iMini) = min(max(Permanencia_30(iMini),0.1),99.9);
        end
        Qpermhist=[Qpermhist Permanencia_30];
    end
end
Qpermhist(:,1)=[];
Permanencia_30=Permanencia_300;

%% Escrevendo a nova saida

Mini=1:nMini;
Mini=Mini';

% Semana_1
varsM=[Mini, Qhoje, Semana_1,Semana_2, Semana_3, Semana_4, Semana_5, ...
    Semana_6, Semana_7,  Semana_8, Permanencia, Permanencia_30, ...
    Permanencia2, Permanencia2_30, Tr, Tr_30, Qmax2, Qmin2, Qpermhist];


anomm=1980:1980+n-2;
str=[];
for i=1:n-1
    str=[str ', ' num2str(anomm(i)) '_max'];
end

str_min=[];
for i=1:n-1
    str_min=[str_min ', ' num2str(anomm(i)) '_min'];
end


varsS={['Mini, Qhoje, Semana_1, Semana_2, Semana_3, Semana_4, Semana_5, Semana_6, Semana_7,  Semana_8, Permanencia, Permanencia_30, Permanencia2, Permanencia2_30, Tr, Tr_30' ...
    str str_min str_meses]}; 
% Tirei porque nao tem no site SARDIM 'Q90', 'Q902', Tr_intervalo
% e TR_30_intervalo dá pra fazer dentro do HTML

dlmwrite('teste.csv',varsS,'delimiter','');
dlmwrite('teste.csv',varsM,'delimiter',',','precision',6,'-append');

%% Escrevendo o dia de ultima atualizacao
date=dates(end);
time = datestr(date,'YYYY-mm-dd');
fileID = fopen('date.txt','w+');
fprintf(fileID,'%10s\n',time);
fclose(fileID);