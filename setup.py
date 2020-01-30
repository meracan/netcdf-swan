from setuptools import setup, find_packages

setup(
    name="netcdf-swan",
    version="0.1",
    packages=find_packages('src', exclude=["tests", "tests.*"]),
    package_dir={'': 'src'},
    include_package_data=True,
    install_requires=[
        'pandas',
        'regex',
        'numpy',
    ],
)
