from setuptools import setup, find_packages


setup(
    name="moldfoam",
    version="0.1",
    packages=find_packages(),
    py_modules=['cli'],
    entry_points={
        'console_scripts': [
            'moldfoam=moldfoam.cli:cli',  # For Click version
        ],
    },
    install_requires=[
        'click',  # Only needed for the Click version
    ],
)