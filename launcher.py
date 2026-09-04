"""
launcher.py - Desktop Entrypoint and Executable Wrapper for Streamlit.
Can be compiled into a single .exe with PyInstaller.
"""

import os
import sys
import time
import webbrowser
import subprocess
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent

def main():
    print("=" * 60)
    print("   Ampullaceana balthica 16S Microbiome Pipeline")
    print("   Iniciando Interface Grafica...")
    print("=" * 60)
    
    app_path = BASE_DIR / "app.py"
    if not app_path.exists():
        print(f"[ERRO] app.py nao encontrado em: {app_path}")
        input("Pressione Enter para sair...")
        sys.exit(1)

    # Launch Streamlit in subprocess
    cmd = [
        sys.executable,
        "-m",
        "streamlit",
        "run",
        str(app_path),
        "--browser.serverAddress",
        "localhost",
        "--server.headless",
        "false",
        "--server.port",
        "8501",
    ]

    print("[INFO] Servidor a iniciar em http://localhost:8501 ...")
    
    # Open browser after a short delay
    def open_browser():
        time.sleep(2)
        webbrowser.open("http://localhost:8501")
        
    import threading
    t = threading.Thread(target=open_browser, daemon=True)
    t.start()

    try:
        subprocess.run(cmd)
    except KeyboardInterrupt:
        print("\nA encerrar o servidor...")

if __name__ == "__main__":
    main()
