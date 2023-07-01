from setuptools import setup, find_packages

setup(
    name='build_complex',
    version='0.1',
    packages= find_packages(),
    author=["qyfu","qcxia"],
    author_email="xiaqiancheng@nibs.ac.cn",
    entry_points={
        "console_scripts": [
            "build_complex = build_complex.build_complex_system:main"
        ]
    }
    # scripts=["db2_converter/build_ligand.py"]
)