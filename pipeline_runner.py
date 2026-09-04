"""
pipeline_runner.py - Backend executor for the Ampullaceana balthica Microbiome Pipeline.
Handles config.yaml management, environment detection (WSL2, SSH, Local),
Snakemake subprocess orchestration, and output packaging.
"""

import os
import sys
import yaml
import shutil
import zipfile
import subprocess
from typing import Generator, Dict, Any, List, Optional
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
CONFIG_PATH = BASE_DIR / "config.yaml"


def load_config(path: Optional[Path] = None) -> Dict[str, Any]:
    """Load configuration from config.yaml."""
    target_path = path or CONFIG_PATH
    if not target_path.exists():
        raise FileNotFoundError(f"Config file not found at: {target_path}")
    with open(target_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def save_config(config_data: Dict[str, Any], path: Optional[Path] = None) -> None:
    """Save configuration data back to config.yaml."""
    target_path = path or CONFIG_PATH
    with open(target_path, "w", encoding="utf-8") as f:
        yaml.dump(config_data, f, default_flow_style=False, sort_keys=False, allow_unicode=True)


def check_environment() -> Dict[str, bool]:
    """Detect available execution environments on the system."""
    env_status = {
        "wsl": False,
        "conda": False,
        "ssh": False,
    }
    # Check WSL
    try:
        res = subprocess.run(["wsl", "-l", "-q"], capture_output=True, text=True, timeout=5)
        if res.returncode == 0 and len(res.stdout.strip()) > 0:
            env_status["wsl"] = True
    except Exception:
        env_status["wsl"] = False

    # Check local conda
    conda_exe = shutil.which("conda")
    if conda_exe:
        env_status["conda"] = True

    # Check SSH client
    ssh_exe = shutil.which("ssh")
    if ssh_exe:
        env_status["ssh"] = True

    return env_status


def load_env_settings() -> Dict[str, str]:
    """Read private settings from local .env if present, or OS environment."""
    env_file = BASE_DIR / ".env"
    settings = {
        "SSH_HOST": os.environ.get("SSH_HOST", ""),
        "SSH_PATH": os.environ.get("SSH_PATH", ""),
        "CONDA_ENV": os.environ.get("CONDA_ENV", "qiime2-amplicon-2026.1"),
    }
    if env_file.exists():
        try:
            with open(env_file, "r", encoding="utf-8") as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith("#") and "=" in line:
                        k, v = line.split("=", 1)
                        settings[k.strip()] = v.strip()
        except Exception:
            pass
    return settings


def build_command(
    phase: int,
    cores: int = 20,
    dry_run: bool = False,
    exec_mode: str = "wsl",
    wsl_distro: str = "Ubuntu",
    ssh_host: str = "",
    remote_path: str = "",
    conda_env: str = "qiime2-amplicon-2026.1",
) -> str:
    """Build the command to execute Snakemake in the specified environment."""
    snakefile = f"Snakefile_phase{phase}"
    dry_flag = "-n" if dry_run else ""
    cores_flag = f"--cores {cores}"
    
    base_cmd = f"snakemake -s {snakefile} --configfile config.yaml {cores_flag} {dry_flag}".strip()

    if exec_mode == "wsl":
        # Run inside WSL with conda activated
        win_path = str(BASE_DIR).replace("\\", "/")
        if ":" in win_path:
            drive, rest = win_path.split(":", 1)
            wsl_path = f"/mnt/{drive.lower()}{rest}"
        else:
            wsl_path = win_path

        wsl_cmd = (
            f"cd '{wsl_path}' && "
            f"source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null || true; "
            f"conda activate {conda_env} 2>/dev/null || true; "
            f"{base_cmd}"
        )
        return f'wsl -d {wsl_distro} -e bash -c "{wsl_cmd}"'

    elif exec_mode == "ssh":
        # Run on remote server via SSH
        remote_cmd = (
            f"cd {remote_path} && "
            f"source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null || true; "
            f"conda activate {conda_env} 2>/dev/null || true; "
            f"{base_cmd}"
        )
        return f'ssh {ssh_host} "{remote_cmd}"'

    else:
        # Local direct run
        return f"conda run -n {conda_env} {base_cmd}"


def run_pipeline_stream(cmd: str) -> Generator[str, None, int]:
    """Execute command and yield stdout/stderr lines in real-time."""
    process = subprocess.Popen(
        cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
        bufsize=1,
    )

    if process.stdout:
        for line in process.stdout:
            yield line
        process.stdout.close()

    returncode = process.wait()
    return returncode


def list_output_results() -> Dict[str, List[Path]]:
    """Scan and list all generated figures, PDFs, tables, and reports."""
    results: Dict[str, List[Path]] = {
        "heattrees": [],
        "alpha_plots": [],
        "pca_plots": [],
        "stats_tables": [],
        "reports": [],
    }

    search_dirs = [BASE_DIR / "Outputs", BASE_DIR / "diversity_results_viz_Ampullaceana", BASE_DIR / "glmm_outputs", BASE_DIR]

    for d in search_dirs:
        if not d.exists():
            continue
        for f in d.iterdir():
            if not f.is_file():
                continue
            fname_lower = f.name.lower()
            if "heattree" in fname_lower:
                results["heattrees"].append(f)
            elif any(k in fname_lower for k in ["shannon", "faith", "observed"]) and f.suffix in [".png", ".pdf", ".txt"]:
                results["alpha_plots"].append(f)
            elif "pca" in fname_lower and f.suffix in [".png", ".pdf", ".txt"]:
                results["pca_plots"].append(f)
            elif f.suffix in [".tsv", ".csv", ".txt"] and ("pvalue" in fname_lower or "summary" in fname_lower or "results" in fname_lower):
                results["stats_tables"].append(f)
            elif f.suffix in [".html", ".qzv"]:
                results["reports"].append(f)

    for k in results:
        seen_names = set()
        unique_list = []
        for p in results[k]:
            if p.name not in seen_names:
                seen_names.add(p.name)
                unique_list.append(p)
        results[k] = sorted(unique_list, key=lambda p: p.name)

    return results


def create_results_zip(zip_filename: str = "Ampullaceana_Resultados.zip") -> Path:
    """Package key result files into a single downloadable ZIP archive."""
    zip_path = BASE_DIR / zip_filename
    res = list_output_results()
    added_names = set()
    
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for category, file_list in res.items():
            for fpath in file_list:
                arcname = f"{category}/{fpath.name}"
                if arcname not in added_names:
                    zf.write(fpath, arcname=arcname)
                    added_names.add(arcname)
        
        if CONFIG_PATH.exists():
            zf.write(CONFIG_PATH, arcname="config_provenance.yaml")

    return zip_path


def sync_results_from_remote(
    ssh_host: str = "",
    remote_path: str = "",
    password: Optional[str] = None,
) -> tuple:
    """Sync output files, figures, and tables from the remote server back to local folder via SFTP."""
    env_cfg = load_env_settings()
    host_str = ssh_host or env_cfg.get("SSH_HOST", "")
    r_path = remote_path or env_cfg.get("SSH_PATH", "")

    if not host_str or not r_path:
        return False, "Host SSH ou caminho remoto não definidos."

    username = None
    port = 22
    hostname = host_str
    if "@" in hostname:
        username, hostname = hostname.split("@", 1)
    if ":" in hostname:
        hostname, port_str = hostname.split(":", 1)
        port = int(port_str)

    local_outputs = BASE_DIR / "Outputs"
    local_outputs.mkdir(exist_ok=True)

    try:
        import paramiko
        client = paramiko.SSHClient()
        client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        
        client.connect(
            hostname=hostname,
            port=port,
            username=username,
            password=password if password else None,
            timeout=10,
            look_for_keys=True,
            allow_agent=True,
        )

        stdin, stdout, stderr = client.exec_command(f'echo "{r_path}"')
        resolved_remote = stdout.read().decode().strip() or r_path

        sftp = client.open_sftp()
        downloaded = []

        # 1. Download files from remote Outputs/
        remote_outputs = f"{resolved_remote}/Outputs"
        try:
            for item in sftp.listdir(remote_outputs):
                r_file = f"{remote_outputs}/{item}"
                l_file = local_outputs / item
                try:
                    sftp.get(r_file, str(l_file))
                    downloaded.append(f"Outputs/{item}")
                except Exception:
                    pass
        except Exception:
            pass

        # 2. Download glmm_pvalues.tsv
        for root_file in ["glmm_pvalues.tsv"]:
            try:
                sftp.get(f"{resolved_remote}/{root_file}", str(BASE_DIR / root_file))
                downloaded.append(root_file)
            except Exception:
                pass

        # 3. Download any HeatTree_*.pdf in root
        try:
            for item in sftp.listdir(resolved_remote):
                if item.startswith("HeatTree_") and item.endswith(".pdf"):
                    sftp.get(f"{resolved_remote}/{item}", str(local_outputs / item))
                    downloaded.append(item)
        except Exception:
            pass

        sftp.close()
        client.close()

        if downloaded:
            return True, f"✓ {len(downloaded)} ficheiros descarregados com sucesso do servidor:\n" + ", ".join(downloaded)
        else:
            return True, "Ligação bem-sucedida, mas nenhum ficheiro novo encontrado na pasta do servidor."

    except Exception as e:
        return False, f"Erro de ligação SSH/SFTP: {e}"
