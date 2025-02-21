from setuptools import setup, find_packages

setup(
    name="anna16",
    version="1.1",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        "mypackage": ["model_files/trimmed/*.zip", "model_files/untrimmed/*.zip"],
    },
)

