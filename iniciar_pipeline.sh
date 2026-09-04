#!/usr/bin/env bash
# ================================================================
# Ampullaceana balthica - 16S Microbiome Pipeline GUI Launcher (Linux/Ubuntu)
# Centro de Ecologia, Evolução e Alterações Ambientais (cE3c)
# ================================================================

set -e

echo "================================================================"
echo "   Ampullaceana balthica - 16S Microbiome Pipeline GUI"
echo "   Centro de Ecologia, Evolucao e Alteracoes Ambientais (cE3c)"
echo "================================================================"
echo ""

# Navegar para a pasta do script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd "$DIR"

# Ativar Conda caso esteja instalado
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
elif [ -f "/opt/conda/etc/profile.d/conda.sh" ]; then
    source "/opt/conda/etc/profile.d/conda.sh"
fi

# Tentar ativar o ambiente QIIME 2 se existir
if command -v conda >/dev/null 2>&1; then
    echo "[INFO] A ativar ambiente Conda (qiime2-amplicon-2026.1)..."
    conda activate qiime2-amplicon-2026.1 2>/dev/null || true
fi

# Verificar Python
if ! command -v python3 >/dev/null 2>&1 && ! command -v python >/dev/null 2>&1; then
    echo "[ERRO] Python nao encontrado no sistema!"
    echo "Por favor, instale o Python 3.9+ ou ative o ambiente Conda com o QIIME 2."
    exit 1
fi

PYTHON_CMD="python3"
if ! command -v python3 >/dev/null 2>&1; then
    PYTHON_CMD="python"
fi

echo "[OK] Python detectado: $($PYTHON_CMD --version)"
echo "[INFO] A verificar dependencias (Streamlit, Pandas, PyYAML, Paramiko)..."
$PYTHON_CMD -m pip install --quiet streamlit pandas pyyaml paramiko

# Descobrir o IP local para informar o utilizador
LOCAL_IP=$(hostname -I 2>/dev/null | awk '{print $1}' || echo "localhost")

echo ""
echo "================================================================"
echo "[INFO] A iniciar a Interface Grafica Streamlit..."
echo "Acesso local no servidor:       http://localhost:8501"
echo "Acesso a partir de outro PC:    http://${LOCAL_IP}:8501"
echo "Para encerrar a aplicacao, prima Ctrl + C nesta janela."
echo "================================================================"
echo ""

$PYTHON_CMD -m streamlit run app.py --server.port 8501 --server.address 0.0.0.0 --server.headless true
