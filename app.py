"""
app.py - Streamlit Web Interface for the Ampullaceana balthica 16S Microbiome Pipeline.
Provides an interactive, reproducible dashboard for configuring, running, and visualizing the pipeline.
"""

import os
import sys
import time
import platform
import pandas as pd
import streamlit as st
from pathlib import Path

import pipeline_runner

# Page configuration
st.set_page_config(
    page_title="Ampullaceana balthica Microbiome Pipeline",
    page_icon="🐌",
    layout="wide",
    initial_sidebar_state="expanded",
)

BASE_DIR = Path(__file__).resolve().parent

# Custom CSS styling for a modern scientific look
st.markdown("""
<style>
    .main-header {
        font-size: 2.2rem;
        font-weight: 700;
        color: #1E3A8A;
        margin-bottom: 0.2rem;
    }
    .sub-header {
        font-size: 1.1rem;
        color: #4B5563;
        margin-bottom: 1.5rem;
    }
    .status-badge-ok {
        background-color: #DEF7EC;
        color: #03543F;
        padding: 4px 10px;
        border-radius: 6px;
        font-weight: 600;
        font-size: 0.85rem;
    }
    .status-badge-warn {
        background-color: #FEF3C7;
        color: #92400E;
        padding: 4px 10px;
        border-radius: 6px;
        font-weight: 600;
        font-size: 0.85rem;
    }
    .metric-card {
        background-color: #F9FAFB;
        border: 1px solid #E5E7EB;
        padding: 15px;
        border-radius: 8px;
        text-align: center;
    }
</style>
""", unsafe_allow_html=True)


# --- SIDEBAR: AMBIENTE E CONTROLO ---
with st.sidebar:
    st.image("https://img.icons8.com/color/96/snail.png", width=64)
    st.markdown("### 🐌 Pipeline Controller")
    st.caption("Ecologia Microbiana & Bioinformática — cE3c")
    
    st.markdown("---")
    st.markdown("#### 🖥️ Ambiente de Execução")
    
    envs = pipeline_runner.check_environment()
    
    col_e1, col_e2 = st.columns(2)
    with col_e1:
        st.write("WSL2 Ubuntu:")
        st.write("SSH Remoto:")
    with col_e2:
        st.markdown("<span class='status-badge-ok'>Ativo</span>" if envs["wsl"] else "<span class='status-badge-warn'>Inativo</span>", unsafe_allow_html=True)
        st.markdown("<span class='status-badge-ok'>Ativo</span>" if envs["ssh"] else "<span class='status-badge-warn'>Inativo</span>", unsafe_allow_html=True)
    
    env_defaults = pipeline_runner.load_env_settings()
    is_linux = platform.system() == "Linux"

    if is_linux:
        mode_options = ["local", "ssh"]
        mode_labels = {
            "local": "💻 Servidor Local (Conda Nativo)",
            "ssh": "🌐 Outro Servidor Remoto (SSH)",
        }
        mode_idx = 0
    else:
        mode_options = ["wsl", "ssh", "local"]
        mode_labels = {
            "wsl": "🐧 WSL2 Ubuntu (Local)",
            "ssh": "🌐 Servidor Remoto (SSH)",
            "local": "💻 Conda Direto (Local)",
        }
        mode_idx = 0 if envs["wsl"] else (1 if envs["ssh"] else 2)

    exec_mode = st.selectbox(
        "Modo de Execução:",
        options=mode_options,
        format_func=lambda x: mode_labels.get(x, x),
        index=mode_idx
    )

    if exec_mode == "wsl":
        wsl_distro = st.text_input("Distribuição WSL:", value="Ubuntu")
        ssh_host = ""
        remote_path = ""
        ssh_pass = ""
    elif exec_mode == "ssh":
        wsl_distro = "Ubuntu"
        ssh_host = st.text_input("Host SSH:", value=env_defaults.get("SSH_HOST", ""), placeholder="utilizador@servidor.instituicao.pt")
        remote_path = st.text_input("Caminho Remoto:", value=env_defaults.get("SSH_PATH", ""), placeholder="~/Ampullaceana-balthica")
        ssh_pass = st.text_input("Palavra-passe SSH (opcional):", type="password", help="Necessária se não tiver chaves SSH configuradas. Nunca é guardada no disco.")
        
        if ssh_host and remote_path:
            if st.button("📥 Importar Outputs do Servidor", use_container_width=True, key="sidebar_sync_remote"):
                with st.spinner("A importar ficheiros do servidor via SFTP..."):
                    ok, sync_log = pipeline_runner.sync_results_from_remote(ssh_host, remote_path, password=ssh_pass)
                    if ok:
                        st.success("Ficheiros importados com sucesso!")
                        time.sleep(1)
                        st.rerun()
                    else:
                        st.error(f"Erro:\n{sync_log}")
    else:
        wsl_distro = "Ubuntu"
        ssh_host = ""
        remote_path = ""
        ssh_pass = ""

    default_conda = env_defaults.get("CONDA_ENV", "qiime2-amplicon-2024.10")
    conda_env = st.text_input("Ambiente Conda QIIME 2:", value=default_conda)
    
    st.markdown("---")
    st.info("💡 **Dica de Reproducibilidade:**\nTodos os parâmetros e execuções são registados no ficheiro `config.yaml` para rastreabilidade científica total.")


# --- MAIN CONTENT TABS ---
st.markdown("<div class='main-header'>Ampullaceana balthica — 16S rRNA Pipeline</div>", unsafe_allow_html=True)
st.markdown("<div class='sub-header'>Interface Interativa para Orquestração, Controlo de Qualidade e Análise Ecológica do Microbioma</div>", unsafe_allow_html=True)

tab_overview, tab_data, tab_params, tab_exec, tab_results = st.tabs([
    "📋 Visão Geral",
    "📁 Dados & Metadados",
    "⚙️ Parâmetros do Pipeline",
    "🚀 Execução do Pipeline",
    "📊 Resultados & Galeria"
])


# --- TAB 1: VISÃO GERAL ---
with tab_overview:
    st.markdown("### 🧬 Arquitetura do Pipeline")
    st.markdown("""
    Este pipeline automatizado foi concebido para processar sequências amplicon do gene **16S rRNA** (região V3-V4) de gastrópodes de água doce (*Ampullaceana balthica*), avaliando o impacto sinérgico de **Temperatura**, **Dieta**, **Fósforo** e **Origem Populacional**.
    """)
    
    col_ov1, col_ov2 = st.columns([1, 1])
    with col_ov1:
        st.markdown("#### 🔵 Fase 1: Bioinformática & Filtragem")
        st.markdown("""
        1. **Controlo de Qualidade Bruto:** `fastp` (remoção de adaptadores Illumina e poly-G).
        2. **Remoção de Primers:** `cutadapt` (primers 341F e 806R).
        3. **Denoising e Geração de ASVs:** `DADA2` (correção de erros, fusão paired-end e quimeras).
        4. **Classificação Taxonómica:** `SILVA 138.2` (classificador Naive Bayes).
        5. **Filtragem Taxonómica:** Remoção rigorosa de mitocôndrias, cloroplastos e não atribuídos.
        """)
        
    with col_ov2:
        st.markdown("#### 🟢 Fase 2: Ecologia & Estatística")
        st.markdown("""
        1. **Alinhamento e Filogenia:** `FastTree` (construção da árvore filogenética).
        2. **Rarefação Normalizada:** Sub-amostragem uniforme a 50.000 reads/amostra.
        3. **Diversidade Alfa & Beta:** Shannon, Faith's PD, Observed Features, UniFrac e Bray-Curtis (PCoA).
        4. **Modelos Lineares Mistos (GLMMs):** R `glmmTMB` com controlo morfológico (`Size_PC1`) e `DHARMa`.
        5. **Abundância Diferencial:** `ANCOM-BC2` e cladogramas hierárquicos (*Heat Trees* com `metacoder`).
        """)

    st.markdown("---")
    st.markdown("#### 🚀 Como Utilizar:")
    st.markdown("""
    1. Aceda à aba **📁 Dados & Metadados** para verificar as amostras FASTQ e a tabela de metadados.
    2. Ajuste parâmetros se necessário na aba **⚙️ Parâmetros do Pipeline**.
    3. Na aba **🚀 Execução**, faça um teste prévio (*Dry-Run*) e inicie a Fase 1 ou Fase 2.
    4. Visualize gráficos e faça o download do pacote completo na aba **📊 Resultados & Galeria**.
    """)


# --- TAB 2: DADOS & METADADOS ---
with tab_data:
    st.markdown("### 📁 Gestão de Ficheiros de Entrada")
    
    config = pipeline_runner.load_config()
    
    col_d1, col_d2 = st.columns(2)
    with col_d1:
        raw_dir = st.text_input("Diretório das Leituras FASTQ:", value=config["raw"]["data_directory"])
        meta_file = st.text_input("Caminho do Ficheiro de Metadados (.tsv):", value=config["raw"]["metadata"])
        
        # Check FASTQs
        raw_path = BASE_DIR / raw_dir
        if raw_path.exists():
            fastqs = list(raw_path.glob("*.fastq.gz"))
            st.success(f"✅ Encontrados **{len(fastqs)} ficheiros FASTQ.gz** ({len(fastqs)//2} amostras estimadas)")
        else:
            st.warning(f"⚠️ Diretório não encontrado no caminho local: `{raw_path}`")
            
    with col_d2:
        tag_r1 = st.text_input("Tag R1 (Forward):", value=config["illumina"]["paired_end_tag_1"])
        tag_r2 = st.text_input("Tag R2 (Reverse):", value=config["illumina"]["paired_end_tag_2"])
        primer_f = st.text_input("Primer Forward (341F):", value=config["illumina"]["primer_f"])
        primer_r = st.text_input("Primer Reverse (806R):", value=config["illumina"]["primer_r"])

    st.markdown("---")
    st.markdown("#### 📋 Pré-visualização e Validação dos Metadados")
    
    meta_path = BASE_DIR / meta_file
    if meta_path.exists():
        try:
            df_meta = pd.read_csv(meta_path, sep="\t")
            st.dataframe(df_meta.head(10), use_container_width=True)
            
            # Validation check
            cols = df_meta.columns.tolist()
            first_col = cols[0]
            st.markdown(f"**Total de Amostras nos Metadados:** {len(df_meta)} linhas | **Coluna de Amostra:** `{first_col}`")
            
            required_cols = ["Temp", "Diet", "Phosphorus", "Pop", "Box"]
            missing_cols = [c for c in required_cols if c not in cols]
            
            if missing_cols:
                st.warning(f"⚠️ Colunas experimentais em falta nos metadados: {missing_cols}")
            else:
                st.success("✅ Todas as variáveis experimentais (`Temp`, `Diet`, `Phosphorus`, `Pop`, `Box`) estão presentes!")
                
        except Exception as e:
            st.error(f"Erro ao ler metadados: {e}")
    else:
        st.error(f"Ficheiro de metadados não encontrado em: `{meta_path}`")


# --- TAB 3: PARÂMETROS DO PIPELINE ---
with tab_params:
    st.markdown("### ⚙️ Ajuste Fino dos Parâmetros")
    st.caption("Ajuste os limiares de corte, filtragem e profundidade de rarefação.")
    
    config = pipeline_runner.load_config()
    
    col_p1, col_p2 = st.columns(2)
    with col_p1:
        st.markdown("#### 🔬 Parâmetros DADA2 (Fase 1)")
        trunc_f = st.number_input("Tamanho de Truncagem Forward (trunc_f):", min_value=0, max_value=350, value=int(config["denoise"]["trunc_f"]), help="0 = sem truncagem (corte já feito por qualidade)")
        trunc_r = st.number_input("Tamanho de Truncagem Reverse (trunc_r):", min_value=0, max_value=350, value=int(config["denoise"]["trunc_r"]), help="0 = sem truncagem")
        max_ee_f = st.slider("Erros Esperados Máximos Forward (max_ee_f):", 0.5, 5.0, float(config["denoise"]["max_ee_f"]), 0.5)
        max_ee_r = st.slider("Erros Esperados Máximos Reverse (max_ee_r):", 0.5, 5.0, float(config["denoise"]["max_ee_r"]), 0.5)
        min_seq = st.number_input("Mínimo de Leituras por Amostra (num_min_seq):", min_value=10, max_value=5000, value=int(config["denoise"]["num_min_seq"]))
        
    with col_p2:
        st.markdown("#### 📊 Parâmetros de Diversidade (Fase 2)")
        raref_depth = st.number_input("Profundidade de Rarefação (raref_depth):", min_value=1000, max_value=200000, value=int(config["diversity"]["raref_depth"]), step=5000, help="Valor onde a curva de Shannon estabiliza")
        max_depth = st.number_input("Profundidade Máxima para Curva (max_depth):", min_value=10000, max_value=500000, value=int(config["diversity"]["max_depth"]), step=10000)
        cpu_threads = st.slider("Número de Threads / Núcleos de CPU:", min_value=1, max_value=32, value=int(config["raw"]["threads"]))
        fixed_effects = st.text_input("Efeitos Fixos (GLMM):", value=config["metadata_columns"]["fixed_effects"])
        random_effect = st.text_input("Efeito Aleatório (GLMM):", value=config["metadata_columns"]["random_effect"])

    col_b1, col_b2 = st.columns(2)
    with col_b1:
        if st.button("💾 Guardar Alterações no config.yaml", use_container_width=True):
            config["denoise"]["trunc_f"] = trunc_f
            config["denoise"]["trunc_r"] = trunc_r
            config["denoise"]["max_ee_f"] = max_ee_f
            config["denoise"]["max_ee_r"] = max_ee_r
            config["denoise"]["num_min_seq"] = min_seq
            config["diversity"]["raref_depth"] = raref_depth
            config["diversity"]["max_depth"] = max_depth
            config["raw"]["threads"] = cpu_threads
            config["metadata_columns"]["fixed_effects"] = fixed_effects
            config["metadata_columns"]["random_effect"] = random_effect
            pipeline_runner.save_config(config)
            st.success("Configurações atualizadas com sucesso em config.yaml!")
            
    with col_b2:
        if st.button("🔄 Restaurar Predefinições Validadas da Tese", use_container_width=True):
            config["denoise"]["trunc_f"] = 0
            config["denoise"]["trunc_r"] = 0
            config["denoise"]["max_ee_f"] = 2.0
            config["denoise"]["max_ee_r"] = 2.0
            config["denoise"]["num_min_seq"] = 100
            config["diversity"]["raref_depth"] = 50000
            config["diversity"]["max_depth"] = 242655
            config["raw"]["threads"] = 20
            pipeline_runner.save_config(config)
            st.info("Valores da tese repostos! (Rarefação: 50.000 reads, Max Depth: 242.655)")
            st.rerun()


# --- TAB 4: EXECUÇÃO DO PIPELINE ---
with tab_exec:
    st.markdown("### 🚀 Orquestrador de Execução Snakemake")
    
    cores = st.number_input("Cores a utilizar nesta execução:", min_value=1, max_value=32, value=20)
    dry_run = st.checkbox("Modo de Simulação (Dry-Run: valida o grafo de regras sem executar tarefas reais)", value=False)
    
    col_ex1, col_ex2 = st.columns(2)
    
    with col_ex1:
        st.markdown("#### 📦 Fase 1: Processamento Inicial & ASVs")
        st.caption("Fastp ➔ Cutadapt ➔ DADA2 ➔ Classificador SILVA ➔ Filtro Taxonómico")
        btn_phase1 = st.button("▶️ Executar Fase 1", type="primary", use_container_width=True)
        
    with col_ex2:
        st.markdown("#### 🌳 Fase 2: Diversidade & Modelos GLMM")
        st.caption("Filogenia ➔ Rarefação ➔ Métricas Core ➔ GLMMs ➔ Heat Trees")
        btn_phase2 = st.button("▶️ Executar Fase 2", type="primary", use_container_width=True)

    # Action handler
    if btn_phase1 or btn_phase2:
        phase = 1 if btn_phase1 else 2
        cmd = pipeline_runner.build_command(
            phase=phase,
            cores=cores,
            dry_run=dry_run,
            exec_mode=exec_mode,
            wsl_distro=wsl_distro,
            ssh_host=ssh_host,
            remote_path=remote_path,
            conda_env=conda_env
        )
        
        st.markdown(f"**Comando Gerado:** `{cmd}`")
        st.markdown("---")
        st.markdown("#### 📟 Consola de Execução ao Vivo:")
        
        log_placeholder = st.empty()
        log_lines = []
        
        with st.spinner(f"A executar Fase {phase}{' (Dry-Run)' if dry_run else ''}..."):
            runner_gen = pipeline_runner.run_pipeline_stream(cmd)
            returncode = -1
            try:
                while True:
                    line = next(runner_gen)
                    log_lines.append(line)
                    # Keep last 50 lines to keep UI responsive
                    display_text = "".join(log_lines[-50:])
                    log_placeholder.code(display_text, language="bash")
            except StopIteration as e:
                returncode = e.value if e.value is not None else 0
                
            if returncode == 0:
                st.success(f"🎉 Fase {phase} concluída com sucesso!")
                if exec_mode == "ssh" and not dry_run:
                    st.info("🔄 A transferir e sincronizar ficheiros gerados no servidor para o computador local...")
                    ok, sync_log = pipeline_runner.sync_results_from_remote(ssh_host, remote_path, password=ssh_pass)
                    if ok:
                        st.success("✅ Resultados sincronizados com sucesso! Aceda à aba 'Resultados & Galeria'.")
                    else:
                        st.warning(f"Aviso na transferência:\n{sync_log}")
            else:
                st.error(f"❌ Execução falhou com código de erro {returncode}. Verifique a consola acima.")


# --- TAB 5: RESULTADOS & GALERIA ---
with tab_results:
    st.markdown("### 📊 Galeria de Resultados & Exportação")
    
    if exec_mode == "ssh":
        col_sync1, col_sync2 = st.columns([3, 1])
        with col_sync1:
            st.info(f"🌐 **Modo Servidor Ativo:** `{ssh_host}` — Os ficheiros gerados remotamente são sincronizados para o teu PC.")
        with col_sync2:
            if st.button("🔄 Sincronizar do Servidor", use_container_width=True, key="btn_manual_sync_ssh"):
                with st.spinner("A transferir gráficos e tabelas do servidor..."):
                    ok, sync_log = pipeline_runner.sync_results_from_remote(ssh_host, remote_path, password=ssh_pass)
                    if ok:
                        st.success("Resultados atualizados com sucesso!")
                        time.sleep(1)
                        st.rerun()
                    else:
                        st.error(f"Erro na sincronização:\n{sync_log}")
    
    res = pipeline_runner.list_output_results()
    
    # Download ZIP Button
    col_z1, col_z2 = st.columns([2, 1])
    with col_z1:
        st.write("Exporte todos os gráficos, ficheiros estatísticos e PDFs com um só clique:")
    with col_z2:
        if st.button("📦 Criar Pacote de Resultados (.ZIP)", use_container_width=True):
            with st.spinner("A compactar ficheiros..."):
                zip_f = pipeline_runner.create_results_zip()
                with open(zip_f, "rb") as fp:
                    st.download_button(
                        label="⬇️ Descarregar Ampullaceana_Resultados.zip",
                        data=fp,
                        file_name="Ampullaceana_Resultados.zip",
                        mime="application/zip",
                        use_container_width=True,
                        key="btn_dl_zip_results"
                    )

    st.markdown("---")
    res_subtab1, res_subtab2, res_subtab3, res_subtab4 = st.tabs([
        "🌳 Árvores de Calor (Heat Trees)",
        "📈 Diversidade Alfa",
        "🧭 Morfologia (PCA Biplots)",
        "📑 Estatística & P-Values"
    ])
    
    with res_subtab1:
        st.markdown("#### 🌳 Cladogramas de Abundância Diferencial (Metacoder)")
        st.caption("Mapeamento filogenético do Log2 Fold Change e testes Wilcoxon (p < 0.05).")
        
        # Display Heat Trees
        ht_files = res.get("heattrees", [])
        if ht_files:
            for idx, ht in enumerate(ht_files):
                if ht.suffix == ".pdf":
                    # Offer PDF download
                    with open(ht, "rb") as f:
                        st.download_button(
                            label=f"📄 Descarregar PDF Vetorial: {ht.name}",
                            data=f,
                            file_name=ht.name,
                            mime="application/pdf",
                            key=f"btn_dl_ht_pdf_{idx}_{ht.name}"
                        )
                elif ht.suffix == ".png":
                    st.image(str(ht), caption=ht.name, use_container_width=True)
        else:
            st.info("Nenhuma Heat Tree encontrada na pasta Outputs/.")
            
    with res_subtab2:
        st.markdown("#### 📈 Boxplots de Diversidade Alfa")
        st.caption("Comparação das métricas de Shannon, Faith's PD e Observed Features entre tratamentos.")
        
        alpha_files = [f for f in res.get("alpha_plots", []) if f.suffix == ".png"]
        if alpha_files:
            for af in alpha_files:
                st.image(str(af), caption=af.name, use_container_width=True)
        else:
            st.info("Nenhum gráfico de diversidade alfa encontrado.")
            
    with res_subtab3:
        st.markdown("#### 🧭 PCA Biplots dos Traços Morfológicos")
        st.caption("Redução dimensional dos traços morfológicos do hospedeiro para gerar Size_PC1.")
        
        pca_files = [f for f in res.get("pca_plots", []) if f.suffix == ".png"]
        if pca_files:
            cols_pca = st.columns(2)
            for i, pf in enumerate(pca_files):
                with cols_pca[i % 2]:
                    st.image(str(pf), caption=pf.name, use_container_width=True)
        else:
            st.info("Nenhum PCA biplot encontrado na pasta Outputs/.")
            
    with res_subtab4:
        st.markdown("#### 📑 Resultados Estatísticos dos Modelos Lineares Mistos (GLMMs)")
        
        pvals_path = BASE_DIR / "glmm_pvalues.tsv"
        if pvals_path.exists() and pvals_path.stat().st_size > 0:
            try:
                df_pval = pd.read_csv(pvals_path, sep="\t")
                if not df_pval.empty:
                    # Format p-values
                    st.dataframe(
                        df_pval.style.map(
                            lambda v: "background-color: #DEF7EC; font-weight: bold;" if isinstance(v, (int, float)) and v < 0.05 else "",
                            subset=["Pr(>|t|)"] if "Pr(>|t|)" in df_pval.columns else []
                        ),
                        use_container_width=True
                    )
                    st.caption("💡 Valores a verde indicam efeitos com significância estatística (p < 0.05).")
            except Exception:
                pass

        st.markdown("---")
        st.markdown("##### 🔬 Relatórios Completos dos Modelos GLMM (ANOVA Tipo II, Coeficientes e Resíduos)")
        metric_opt = st.selectbox(
            "Selecione a Métrica de Diversidade para Inspecionar:",
            options=["Shannon Entropy", "Observed Features", "Faith's Phylogenetic Diversity"],
            index=0
        )
        file_map = {
            "Shannon Entropy": "shannon_entropy_summary_results.txt",
            "Observed Features": "observed_features_summary_results.txt",
            "Faith's Phylogenetic Diversity": "faith_pd_summary_results.txt",
        }
        
        chosen_file = BASE_DIR / file_map[metric_opt]
        if chosen_file.exists() and chosen_file.stat().st_size > 0:
            with open(chosen_file, "r", encoding="utf-8", errors="ignore") as fp:
                st.code(fp.read(), language="text")
        else:
            st.info(f"O ficheiro {file_map[metric_opt]} ainda não foi gerado ou está vazio.")
