import setuptools
from src import __version__


short_description = """Binary classification for predicting whether proteins are secreted
by the bacterial Type III secretion system based on a model 
trained on sequence- and amino acid property (patterns)-based feature groups."""

with open("README.md", "r", encoding="utf-8") as ifile:
    long_description = ifile.read()


setuptools.setup(
    name="EffectiveT3",
    version=__version__,
    author="Nicolas Robert Nemeth",
    author_email="nicolas.r.nemeth@gmail.com",
    description=short_description,
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/nicolasrnemeth/EffectiveT3",
    download_url="https://github.com/nicolasrnemeth/EffectiveT3/archive/refs/tags/beta.tar.gz",
    project_urls={
        "Bug Tracker": "https://github.com/nicolasrnemeth/EffectiveT3/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    license_files=["LICENSE.txt"],
    packages=setuptools.find_packages(),
    package_data={
        "": [
            "training_config.yaml",
            "protein_sequences/*",
            "src/models/*", 
            "src/training/hyperparameter_space.json",
            "src/training/models_and_parameters"
        ],
    },
    install_requires=[
        "numpy>=1.21.5",
        "lightgbm>=3.3.2",
        "sklearn-genetic>=0.5.1",
        "scikit-learn>=1.1.2",
        "pyyaml>=6.0"
    ],
    include_package_data=True,
    python_requires=">=3.6",
    entry_points={
        "console_scripts": [
            "effectivet3 = src.__predict__:main",
            "effectiveTrain = src.__train__:main",
        ],
    }
)
