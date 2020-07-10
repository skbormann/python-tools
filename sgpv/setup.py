import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="sgpv",
    version="1.0.3",
    author="Sven-Kristjan Bormann",
    author_email="sven-kristjan@gmx.de",
    description="Tools to calculate SGPVs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/skbormann/python-tools/sgpv",
    keywords='"second generation p-values" sgpv "false discovery risk" \
        "false confirmatory risk" power-function',
    packages=setuptools.find_packages(),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
    install_requires=['pandas>=1.0.4', 'matplotlib>=3.2.1',
                      'numpy>=1.18.0', 'scipy>=1.3.2']
)
