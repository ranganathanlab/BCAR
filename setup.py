from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize

setup(
    name="BCAR",
    version="1.14.0",
    description="Barcode Collapse by Aligning Reads",
    author="Bryan Andrews",
    author_email="andrewsb@uchicago.edu",
    packages=find_packages(),
    ext_modules=cythonize(
        [Extension("BCAR.needleman_wunsch", ["BCAR/needleman_wunsch.pyx"])],
        compiler_directives={"language_level": "3"}
    ),
    python_requires=">=3.6",
    install_requires=[
        "Cython",
    ],
    entry_points={
        "console_scripts": [
            "bcar=BCAR.BCAR_v1_14:main",  # Use the fully qualified package name
        ]
    },
)


