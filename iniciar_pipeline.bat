@echo off
title Ampullaceana balthica 16S Pipeline
color 0A

echo ================================================================
echo    Ampullaceana balthica - 16S Microbiome Pipeline GUI
echo    Centro de Ecologia, Evolucao e Alteracoes Ambientais (cE3c)
echo ================================================================
echo.

:: Verificar se o Python esta instalado
where python >nul 2>nul
if %errorlevel% neq 0 (
    echo [ERRO] Python nao encontrado no sistema!
    echo Por favor, instale o Python 3.9+ ou adicione o Python ao PATH do Windows.
    echo.
    pause
    exit /b 1
)

echo [OK] Python detectado com sucesso.
echo [INFO] A verificar dependencias (Streamlit, Pandas, PyYAML)...
python -m pip install --quiet streamlit pandas pyyaml

echo.
echo ================================================================
echo [INFO] A iniciar a Interface Grafica no navegador predefinido...
echo O servidor estara acessivel em: http://localhost:8501
echo Para fechar a aplicacao, encerre esta janela.
echo ================================================================
echo.

python -m streamlit run app.py --browser.serverAddress localhost --server.headless false

pause
