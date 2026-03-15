from setuptools import setup, find_packages

setup(
    name="doseprojection",
    version="0.1.0",
    description="Preclinical dose projection tool for animal studies and human dose estimation",
    author="jmabbott40",
    url="https://github.com/jmabbott40/doseprojection",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20",
        "pandas>=1.3",
        "scipy>=1.7",
        "openpyxl>=3.0",
        "matplotlib>=3.4",
    ],
    extras_require={
        "dev": [
            "pytest>=7.0",
            "jupyter>=1.0",
        ],
    },
)
