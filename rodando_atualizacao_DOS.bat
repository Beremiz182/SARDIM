rem !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
rem
rem     SARDIM
rem
rem	Hidrologia de Grande Escala (HGE)
rem     IPH/UFRGS
rem     Agosto - 2020
rem
rem           
rem rem -> É o comando utilizado para comentar o código
rem Esse arquivo consiste numa tentativa de criar um código .bat
rem Que atualize automaticamente e diariamente o SARDIM - Informação da Vazão até o dia anterior
rem link do website: https://sardim.herokuapp.com/
rem
rem !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

echo off
cls
echo INICIANDO


rem !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
rem Chamamos a unidade D:/ onde estão os arquivos

D:


rem !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
rem Mudamos o endereço e colocamos para rodar uma rotina no Matlab
rem Rotina que baixa e prepara os dados de chuva para o modelo hidrologico

cd D:\JP\ANA\extensao_chuva\SARDIMv2\matlab
echo Baixando, Interpolando e Removendo Vies de Chuva (MERGE) pelo Matlab
"C:\Program Files\Polyspace\R2021a\bin\matlab.exe" -wait -nodisplay -nosplash -nodesktop -r "atualizacao_chuva_unificado_MERGE; exit"


rem !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
rem Agora vamos rodar o MGB

echo Rodando o MGB
cd D:\JP\ANA\extensao_chuva\SARDIMv2\MGB_SA
MGBSouthAmerica.exe


rem !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
rem Acabando o MGB, vamos juntar a vazao dos dois períodos

cls
echo O MGB acabou de rodar


rem !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
rem Aqui fazemos o calculo das estatísticas em cima das vazoes do MGB

echo Fazendo o calculo das estatisticas do MGB
cd D:\JP\ANA\extensao_chuva\SARDIMv2\matlab
"C:\Program Files\Polyspace\R2021a\bin\matlab.exe" -wait -nodisplay -nosplash -r "SARDIM3; exit"


rem !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
rem Agora nos vamos fazer o upload do SARDIM no site (o mais demorado)

"C:\Program Files (x86)\WinSCP\winscp.com" /command "open sftp://sardim:ENEg2WnSl56E@kepler.ufrgs.br:9022 -hostkey=""ssh-ed25519 255 tYI4jgr0yleRPzrZdh37/m462/pMKzVwAqMCVwFJ9vY""" "put -latest D:\JP\ANA\extensao_chuva\SARDIMv2\matlab\teste.csv /matlab/" "exit"
"C:\Program Files (x86)\WinSCP\winscp.com" /command "open sftp://sardim:ENEg2WnSl56E@kepler.ufrgs.br:9022 -hostkey=""ssh-ed25519 255 tYI4jgr0yleRPzrZdh37/m462/pMKzVwAqMCVwFJ9vY""" "put -latest D:\JP\ANA\extensao_chuva\SARDIMv2\matlab\date.txt /matlab/" "exit"
start https://www.ufrgs.br/sardim/
echo Entrando no Sardim

echo on
rem !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
