import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()
requirements = ['requests==2.27.1', 'numpy==1.21.2', 'scikit-learn==1.0.1']
setuptools.setup(
    name='protFeat',
    version='1.0',
    author='Gokhan Ozsari',
    author_email='gozsari@ceng.metu.edu.tr',
    description='Installation of Package',
    long_description="readme",
    long_description_content_type="text/markdown",
    url='https://github.com/gozsari/protFeat',
    project_urls = {
        "Bug Tracker": "https://github.com/gozsari/protFeat/issues"
    },
    install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.7",
)
