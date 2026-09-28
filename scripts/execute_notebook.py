"""Execute the walkthrough with the current Python, without a global kernel install."""
import os
from pathlib import Path
import sys
import nbformat
from nbclient import NotebookClient
from jupyter_client.kernelspec import KernelSpecManager

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("JUPYTER_RUNTIME_DIR", str(ROOT / ".tools/jupyter-runtime"))
os.environ.setdefault("IPYTHONDIR", str(ROOT / ".tools/ipython"))
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".tools/matplotlib"))


class CurrentPythonKernel(KernelSpecManager):
    def get_kernel_spec(self, kernel_name):
        spec = super().get_kernel_spec(kernel_name)
        spec.argv = [sys.executable, "-m", "ipykernel_launcher", "-f", "{connection_file}"]
        return spec


path = ROOT / "notebooks/01_oscc_analysis.ipynb"
notebook = nbformat.read(path, as_version=4)
client = NotebookClient(notebook, timeout=180, kernel_name="python3", resources={"metadata": {"path": str(ROOT)}},
                        kernel_manager_class="jupyter_client.manager.KernelManager")
client.create_kernel_manager()
client.km.kernel_spec_manager = CurrentPythonKernel()
client.execute()
nbformat.write(notebook, path)
print(f"Executed {sum(c.cell_type == 'code' for c in notebook.cells)} code cells without errors")
