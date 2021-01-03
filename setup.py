import setuptools
import os


readme_dir = os.path.dirname(__file__)
readme_path = os.path.join(readme_dir, 'README.md')
with open(readme_path, "r") as f:
    long_description = f.read()


required_packages = [
    "scikit-learn==0.23.2",
    "numpy",
    "pandas",
    "joblib",
    'argparse',
    'pandarallel'
]


setuptools.setup(
    name="mhclovac",
    version="4.0",
    author="Stefan Stojanovic",
    author_email="stefs304@gmail.com",
    description="MHC binding prediction based on modeled physicochemical properties of peptides",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/stefs304/mhclovac",
    packages=[
        'mhclovac',
    ],
    package_data={
        'mhclovac': ['index/*', 'models/*']
    },
    install_requires=required_packages,
    entry_points={
          'console_scripts': [
              "mhclovac=mhclovac.run:run",
          ]
    },
    classifiers=(
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: MIT License",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Operating System :: POSIX :: Linux"
    )
)
