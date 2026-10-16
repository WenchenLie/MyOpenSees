import shutil
import subprocess
import sys
import sysconfig
from pathlib import Path

from Cython.Build import cythonize
from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


ROOT = Path(__file__).resolve().parent
if sys.version_info < (3, 11):
    raise RuntimeError("Stable ABI build requires Python 3.11 or newer")


class BuildExtToExt(build_ext):
    def initialize_options(self):
        super().initialize_options()
        self.build_lib = str(ROOT / "ext")

    def get_ext_filename(self, ext_name):
        filename = Path(super().get_ext_filename(ext_name))
        extension = next((ext for ext in self.extensions if ext.name == ext_name), None)
        if extension is not None and extension.py_limited_api:
            # Windows importlib recognizes .pyd, but not .abi3.pyd.
            suffix = ".pyd" if sys.platform == "win32" else f".abi3{filename.suffix}"
            return str(filename.with_name(
                f"{ext_name.rsplit('.', 1)[-1]}{suffix}"
            ))
        return str(filename)

    def run(self):
        # On Windows, use an available MinGW compiler without requiring an
        # elevated Visual Studio environment. An explicit MSVC choice still
        # uses setuptools' normal build path.
        gcc = shutil.which("gcc") if sys.platform == "win32" else None
        if gcc is None or self.compiler not in (None, "mingw32"):
            super().run()
        else:
            include = sysconfig.get_paths()["include"]
            for extension in self.extensions:
                output = Path(self.get_ext_fullpath(extension.name))
                output.parent.mkdir(parents=True, exist_ok=True)
                command = [
                    gcc, "-shared", "-O2", "-DPy_LIMITED_API=0x030B0000",
                    f"-I{include}", *extension.sources,
                    f"-L{sys.base_prefix}", "-lpython3", "-o", str(output),
                ]
                subprocess.run(command, check=True)

        if sys.platform == "win32":
            platform_tag = sysconfig.get_config_var("EXT_SUFFIX").split("-", 1)[1].removesuffix(".pyd")
            for extension in self.extensions:
                output = Path(self.get_ext_fullpath(extension.name))
                for alias in output.parent.glob(f"{extension.name}.cp*-{platform_tag}.pyd"):
                    shutil.copyfile(output, alias)


extension = Extension(
    name="TSSCB",
    sources=[str(ROOT / "src" / "TSSCB" / "TSSCB.pyx")],
    define_macros=[("Py_LIMITED_API", "0x030B0000")],
    py_limited_api=True,
    libraries=["python3"] if sys.platform == "win32" else [],
)

setup(
    name="TSSCB",
    ext_modules=cythonize(
        [extension],
        build_dir=str(ROOT / "build" / "cython"),
        language_level="3",
    ),
    cmdclass={"build_ext": BuildExtToExt},
    zip_safe=False,
)

# Build: python setup.py build_ext
# Output: ext/xxx.pyd (Windows, stable ABI) or ext/xxx.abi3.so (Linux/macOS)
