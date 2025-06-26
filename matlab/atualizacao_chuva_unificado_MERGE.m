%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%        Script para atualizar a chuva do MGB AS com dados do MERGE          %%%
%%%%%%%%%%%%             Joao Paulo L F Breda - IPH/UFRGS             %%%%%%%%%%%%
%%%%                                Maio 2020                                 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Eh necessario os seguintes arquivos e rotinas:
%
% downloadData.m: responsavel por baixar os dados do MERGE
% IDW_MERGE.m: responsavel por ler e interpolar os dados do MERGE
%
% chuvabin_MSWEP.pbi: conjunto de dados de precipitacao pronto desde 1980
%                    que faz a mescla de MSWEP e MERGE
% gamma_dist.mat: arquivo com os valores da distribuicao gamma para a
%                 correcao de vies
% mini.gtp: arquivo necessario para fazer a interpolacao espacial dos dados
%           que informa a posicao geografica das minibacias

clear 
clc
close all

run '.\nctoolbox-master\setup_nctoolbox.m'

% indica a pasta onde estao (ou serao) localizados os dados baixados do
% MERGE
foldername='..\files\';
nMini=33749; % numero de minibacias do modelo

%***************************************************************************
%% baixa os dados ate os dias mais atuais

downloadData2
disp('Dados do MERGE baixados');

%***************************************************************************
%% chuvabin anterior mais chuvabin atual

% a secao abaixo junta os chuvabins do MSWEP e do MERGE apos 2015
% se os dados do MERGE e do MSWEP nao tiverem sido mesclados ainda, dai
% entra no if (primera rodada) 
% se jah tiver um chuvabin com os dados em conjunto, dai entra no else

% verificar se jah existe um chuvabin_tudo.pbi, 
if exist('chuvabin_tudo.pbi','file')==2

    fid=fopen('chuvabin_tudo.pbi');% mswep versao 1-1990/2009
    a=fread(fid,'single');
    nT=length(a)/nMini;
    P=reshape(a,nMini,nT);
  
else
    %% abrindo os dados de chuva do MSWEP
    fid=fopen('chuvabin_MSWEP.pbi');% mswep versao 1-1990/2009
    a=fread(fid,'single');
    nT=length(a)/nMini;
    P=reshape(a,nMini,nT);

end

iniP=datenum(1979,1,1);
ultimo_dia=iniP+nT-1;
datevec(ultimo_dia)
    
disp('Dados de chuva lidos')

%**************************************************************************
%% abrindo e interpolando os novos dados de chuva

load dist_matfile
% ultimo dia - ultimo dia com dados
% foldername - pasta com os arquivos de chuva
% dist - vetor de pesos da grade do MGB pra grade do MERGE
MERGE=IDW_MERGE(ultimo_dia,foldername,dist);


%**************************************************************************
%% Tirando o vies dos dados

load gamma_dist_MERGE_climatology

iniM=ultimo_dia+1;
nT=size(MERGE,2);
datas=datevec(iniM:(iniM+nT-1));
mes=datas(:,2);
mes(mes<=2)=mes(mes<=2)+12;

MERGE_bc=0.*MERGE;

for j=1:4
    g=find(mes<=j*3+2);
    if isempty(g)
        continue
    end
    mes(g)=mes(g)+12;
    for im=1:nMini
        disp([im j])
        %for it=1:nT
        % probabilidade da chuva do MERGE não ser superada
        prob=gamcdf(MERGE(im,g),coeffMERGE(im,1,j),coeffMERGE(im,2,j));
        if prob~=1
        % chuva com o vies corrigido utilizando a probabilidade anterior na
        % curva de distribuição do MSWEP
            MERGE_bc(im,g)=gaminv(prob,coeffMSWEP(im,1,j),coeffMSWEP(im,2,j));
        else
            MERGE_bc(im,g)=MERGE(im,g)*2;
        end
        %end
    end
end


%**************************************************************************
%% Juntando e gravando os dados

Ptudo=[P MERGE_bc];

fid=fopen('chuvabin_tudo.pbi','w+');
fwrite(fid,Ptudo,'single');


nT=size(Ptudo,2);
% escrevendo arquivo com a ultima data para simulacao
fid=fopen('nT.sim','w+');
fprintf(fid,'%8i',nT);

fclose all