import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="sgpv-skbormann", # Replace with your own username
    version="1.0.0",
    author="Sven-Kristjan Bormann",
    author_email="sven-kristjan@gmx.de",
    description="Tools to calculate SGPVs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/skbormann/python-tools/sgpv",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)