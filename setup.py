import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="recount-methylation-server-metamaden",
    version="0.0.1",
    author="Sean Maden",
    author_email="maden.sean@gmail.com",
    description="Backend management of DNAm experiment file downloads.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/metamaden/recount-methylation-server",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=setuptools.find_packages(),
    python_requires='>=3.6',
)